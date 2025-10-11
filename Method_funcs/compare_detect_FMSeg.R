# NotSBS

cusum.fts <- function(x, r = NULL,
                      V.diag = TRUE, lrv = TRUE,
                      n = NULL, trim, lbd = NULL,
                      svd_engine = c("auto","irlba","xtx","svd"),
                      svd_tol = 1e-6, svd_maxit = 1000, svd_seed = 1L) {
  svd_engine <- match.arg(svd_engine)
  p <- nrow(x); Time <- ncol(x)
  if (is.null(r))    r    <- median(abc.factor.number(x)$r[4:6])
  if (is.null(lbd))  lbd  <- round(Time^(max(2/5, 1 - min(1, log(p)/log(Time)))) * log(Time)^1.1)
  if (is.null(n))    n    <- floor(Time^0.25)
  
  # === FAST FACTORS ===
  use_irlba <- (svd_engine %in% c("auto","irlba")) && requireNamespace("irlba", quietly = TRUE)
  if (use_irlba && (svd_engine == "irlba" || p >= Time || svd_engine == "auto")) {
    set.seed(as.integer(svd_seed))
    sv <- irlba::irlba(x, nv = r, nu = 0, tol = svd_tol, maxit = svd_maxit)
    f_hat <- t(sv$v) * sqrt(Time)                 # r x T
  } else if (svd_engine %in% c("auto","xtx") || p >= Time) {
    XtX <- crossprod(x)                            # T x T
    eg  <- eigen(XtX, symmetric = TRUE)
    f_hat <- t(eg$vectors[, 1:r, drop = FALSE]) * sqrt(Time)
  } else {
    sv <- svd(x, nu = r, nv = 0)
    f_hat <- (t(sv$u) %*% x) / sqrt(p)            # r x T
  }
  
  ind      <- lower.tri(diag(r), diag = TRUE)
  covF     <- f_hat %*% t(f_hat) / Time
  mean_vec <- covF[ind]
  d        <- r * (r + 1) / 2
  ff       <- matrix(0, nrow = d, ncol = Time)
  
  for (t in seq_len(Time)) {
    ft <- f_hat[, t, drop = FALSE]
    ff[, t] <- tcrossprod(ft)[ind] - mean_vec
  }
  
  V <- ff %*% t(ff) / Time
  if (lrv && n >= 1) {
    for (ell in 1:n) {
      tmp <- ff[, 1:(Time-ell), drop = FALSE] %*% t(ff[, (1:(Time-ell)) + ell, drop = FALSE]) / Time
      w   <- 1 - ell/(n + 1)
      V   <- V + w * (tmp + t(tmp))
    }
  }
  
  # whitening (diag or full)
  if (V.diag) {
    dd <- diag(V)
    posmin <- min(dd[dd > 0])
    dd <- pmax(dd, posmin)
    Vhalf_inv <- diag(1/sqrt(dd))
  } else {
    eig <- eigen(V, symmetric = TRUE)
    dd  <- pmax(eig$values, min(eig$values[eig$values > 0]))
    Vhalf_inv <- diag(1/sqrt(dd), d) %*% t(eig$vectors)
  }
  
  D <- Vhalf_inv %*% ff
  
  intervals <- seeded_intervals(Time, minl = lbd)
  results <- data.frame()
  for (i in seq_len(nrow(intervals))) {
    st <- intervals$st[i]; ed <- intervals$ed[i]
    if ((ed - st) > 2 * trim) {
      res <- find_single_cp(D, st, ed, trim)
      results <- rbind(results, res)
    }
  }
  results
}




find_single_cp <- function(D, st, ed, trim) {
  
  Ttot <- ncol(D)
  
  st_i <- max(1, as.integer(st))
  ed_i <- min(Ttot, as.integer(ed))
  len  <- ed_i - st_i + 1
  stopifnot(len > 2 * trim + 1)
  
  
  stat <- numeric(len)
  total <- rowSums(D[, st_i:ed_i, drop = FALSE])
  lsum  <- numeric(nrow(D))
  
  for (k in seq_len(len - 1)) {
    lsum <- lsum + D[, st_i + k - 1]
    rsum <- total - lsum
    nL   <- k
    nR   <- len - k
    diff <- rsum/nR - lsum/nL
    stat[k] <- sqrt((nL * nR) / len) * sqrt(sum(diff^2))
  }
  
  stat[seq_len(trim)]         <- 0
  stat[(len - trim + 1):len] <- 0
  
  best_rel <- which.max(stat)
  best_abs_i <- st_i + best_rel - 1
  
  data.frame(
    est.cp = best_abs_i,     
    val    = stat[best_rel],
    st     = st,
    ed     = ed,
    trim   = trim
  )
}

NotSBS <- function(x, r = NULL, type = NULL, m = 3, trim, threshold = NULL, 
                   method = c("fixed", "oracle"), 
                   V.diag = TRUE, lrv = TRUE, lbd = NULL) {
  
  method <- match.arg(method)
  Time <- ncol(x)
  
  if (is.null(trim)) {
    stop("Must supply 'trim'.")
  }
  
  results <- cusum.fts(x, r, V.diag = V.diag, lrv = lrv, trim = trim, lbd = lbd)
  results <- results[order(results$ed - results$st, decreasing = FALSE), ]
  
  # Method!
  if (method == "fixed") {
    
    if (is.null(threshold)) {
      stop("Must supply 'threshold'.")
    }
    
    results <- results[results$val >= threshold + 5e-3, ]
    
    cid <- 0
    selected_indices <- integer()
    while (cid < nrow(results)) {
      cid <- cid + 1
      selected_indices <- c(selected_indices, cid)
      
      if (cid < nrow(results)) {
        tmp <- results[(cid + 1):nrow(results), ]
        tmp <- tmp[!(tmp$st < results$est.cp[cid] & tmp$ed >= results$est.cp[cid]), ]
        results <- rbind(results[1:cid, ], tmp)
      }
    }
    
    results <- results[selected_indices, ]
    results <- results[order(results$est.cp), ]
    
    
  } else if (method == "oracle") {
    if (is.null(m)){
      if (identical(type, 0)) { m <- 0 } else { m <- length(type) }}
    
    thds <- sort(results$val, decreasing = TRUE)
    selected_threshold <- NA
    next_highest_threshold <- NULL  
    no_change_threshold <- max(thds, na.rm = TRUE)
    
    selected_solution <- NULL
    
    for (i in seq_along(thds)) {
      thd <- thds[i]
      aux <- results[results$val >= thd, ]
      cid <- 0
      
      while (cid < m && cid < nrow(aux)) {
        cid <- cid + 1
        if (cid < nrow(aux)) {
          tmp <- aux[(cid + 1):nrow(aux), ]
          tmp <- tmp[!(tmp$st < aux$est.cp[cid] & tmp$ed >= aux$est.cp[cid]), ]
          aux <- rbind(aux[1:cid, ], tmp)
        }
      }
      
      if (cid == m) {
        selected_threshold <- thd
        selected_solution <- aux
        S <- selected_solution$est.cp
        
        # === Updated logic for next_highest_threshold ===
        for (j in (i + 1):length(thds)) {
          candidate_thd <- thds[j]
          candidate_intervals <- results[results$val == candidate_thd, ]
          
          for (k in seq_len(nrow(candidate_intervals))) {
            intvl <- candidate_intervals[k, ]
            st <- intvl$st
            ed <- intvl$ed
            if (all(S - st < log(Time) | ed - S < log(Time))) {
              next_highest_threshold <- candidate_thd
              break
            }
          }
          
          if (!is.null(next_highest_threshold)) break
        }
        
        
        results <- selected_solution
        break
      }
    }
    
    results$selected_threshold <- selected_threshold
    results$next_highest_threshold <- next_highest_threshold
    results$no_change_threshold <- no_change_threshold
    results <- results[order(results$est.cp), ]
  }
  
  
  # results <- post_process_cusum(results, nrow(results), F_hat, V = V, post = TRUE)
  
  rownames(results) <- NULL
  return(results)
}


abc.factor.number <- function(x, r.max = NULL, center = TRUE, 
                              p.seq = NULL, n.seq = NULL, do.plot = FALSE) {
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
  xx <- x - mean.x
  
  if(is.null(r.max)) r.max <- min(50, floor(sqrt(min(n - 1, p))))
  
  if(is.null(p.seq)) p.seq <- floor(4 * p / 5 + (1:10) * p / 50)
  if(is.null(n.seq)) n.seq <- floor(4 * n / 5 + (1:10) * n / 50)
  const.seq <- seq(.01, 3, by = 0.01)
  IC <- array(Inf, dim = c(r.max + 1, length(const.seq), 10, 6))
  
  for(kk in 1:min(length(n.seq), length(p.seq))) {
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    
    int <- sort(sample(n, nn, replace = FALSE))
    
    pen <- c((nn + pp) / (nn * pp) * log(nn * pp / (nn + pp)),
             (nn + pp) / (nn * pp) * log(min(nn, pp)),
             log(min(nn, pp)) / min(nn, pp))
    
    covx <- xx[, int] %*% t(xx[, int]) / nn
    
    sv <- svd(covx[1:pp, 1:pp], nu = 0, nv = 0)
    tmp <- rev(cumsum(rev(sv$d))) / pp
    if(pp > r.max) tmp <- tmp[1:(r.max + 1)]
    for (jj in 1:length(const.seq)) {
      for (ic.op in 1:3) {
        IC[1:length(tmp), jj, kk, ic.op] <-
          tmp + (1:length(tmp) - 1) * const.seq[jj] * pen[ic.op]
        IC[1:length(tmp), jj, kk, 3 * 1 + ic.op] <-
          log(tmp) + (1:length(tmp) - 1) * const.seq[jj] * pen[ic.op]
      }
    }
  }
  
  r.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(r.mat, c(1, 3), var)
  r.hat <- rep(0, 6)
  for(ii in 1:6){
    ss <- Sc[, ii]
    if(min(ss) > 0){
      r.hat[ii] <- min(r.mat[max(which(ss == min(ss))),, ii]) - 1
    } else{
      if(sum(ss[-length(const.seq)] != 0 & ss[-1] == 0)) {
        r.hat[ii] <-
          r.mat[which(ss[-length(const.seq)] != 0 &
                        ss[-1] == 0)[1] + 1, dim(r.mat)[2], ii] - 1
      }else{
        r.hat[ii] <- min(r.mat[max(which(ss == 0)),, ii]) - 1
      }
    }
  }
  
  out <- list(r.hat = r.hat)
  attr(out, "data") <- list(Sc = Sc, const.seq = const.seq, r.mat = r.mat)
  
  if(do.plot){
    data <- attr(out, "data")
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    
    par(mfrow = c(2, 3))
    Sc <- data$Sc
    const.seq <- data$const.seq
    q.mat <- data$q.mat
    for(ii in 1:6){
      plot(const.seq, r.mat[, dim(r.mat)[2], ii] - 1, type = "b", pch = 1, col = 2, 
           bty = "n", axes = FALSE, xlab = "constant", ylab = "", main = paste("IC ", ii))
      box()
      axis(1, at = pretty(range(const.seq)))
      axis(2, at = pretty(range(r.mat[, dim(r.mat)[2], ii] - 1)),
           col = 2, col.ticks = 2, col.axis = 2)
      par(new = TRUE)
      plot(const.seq, Sc[, ii], col = 4, pch = 2, type = "b", 
           bty = "n", axes = FALSE, xlab = "", ylab = "")
      axis(4, at = pretty(range(Sc[, ii])), col = 4, col.ticks = 4, col.axis = 4)
      legend("topright", legend = c("r", "Sc"), col = c(2, 4), lty = c(1, 1), pch = c(1, 2), bty = "n")
    }    
  }
  
  return(out)
  
}



