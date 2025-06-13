# lower-triangular vectorisation
vech <- function(M) M[lower.tri(M, diag = TRUE)]

find_single_cp <- function(D, st, ed, trim) {
  Time <- ncol(D)
  st_i <- max(1, as.integer(st) + 1)
  ed_i <- min(Time, as.integer(ed))
  len <- ed_i - st_i + 1
  #stopifnot(len > 2 * trim + 1)
  
  stat <- numeric(len)
  total <- rowSums(D[, st_i:ed_i, drop = FALSE])
  lsum <- numeric(nrow(D))
  
  for (k in seq_len(len - 1)) {
    lsum <- lsum + D[, st_i + k - 1]
    rsum <- total - lsum
    nL <- k; nR <- len - k
    diff <- rsum / nR - lsum / nL
    stat[k] <- sqrt((nL * nR)/len) * base::norm(matrix(diff, 1), "2")
  }
  #stat[seq_len(trim)] <- 0
  #stat[(len - trim + 1L):len] <- 0
  
  best_rel <- which.max(stat)
  data.frame(est.cp = st_i + best_rel - 1,
             val    = stat[best_rel],
             stringsAsFactors = FALSE)
}

# candidates under SBS, using stacked cusum
cand_sbs <- function(G, G_dim,
                     trim = NULL,
                     V.diag = TRUE,
                     lrv = TRUE,
                     n = NULL,
                     lbd = NULL,
                     single = FALSE) {
  
  K <- length(G_dim) # how many mode
  Time <- dim(G)[K + 1]
  
  if (is.null(trim)) trim <- 2 * round(log(Time))
  trim <- as.integer(trim);  stopifnot(trim >= 0)
  if (is.null(n)) n <- floor(Time^0.25)
  
  if (single) {
    intervals <- data.frame(st = 0, ed = Time)
  } else {
    if (is.null(lbd)) lbd <- round(2 * log(Time))
    intervals <- seeded_intervals(Time, minl = lbd)
  }
  
  D_list <- vector("list", K)
  for (i in seq_len(K)) {
    p_i <- G_dim[i]
    d_i <- p_i * (p_i + 1) / 2
    FF <- matrix(0, d_i, Time)
    
    for (t in seq_len(Time)) {
      if (K == 1) {
        g_t <- G[, t]
        S <- g_t %o% g_t
      } else {
        idx <- c(rep(list(TRUE), K), list(t))
        g_t <- do.call(`[`, c(list(G), idx))
        perm <- c(i, setdiff(seq_len(K), i))
        X_i <- matrix(aperm(g_t, perm), nrow = p_i)
        S <- X_i %*% t(X_i)
      }
      FF[, t] <- vech(S)
    }
    
    Z <- FF - rowMeans(FF)
    
    ## long-run variance HAC
    V <- Z %*% t(Z) / Time
    if (lrv && n >= 1) {
      for (ell in seq_len(n)) {
        C <- Z[, 1:(Time - ell)] %*% t(Z[, (1:(Time - ell)) + ell]) / Time
        w <- 1 - ell / (n + 1)
        V <- V + w * (C + t(C))
      }
    }
    if (V.diag) {
      d <- diag(V); d[d <= 0] <- min(d[d > 0])
      Vinvhalf <- diag(1 / sqrt(d))
    } else {
      eg <- eigen(V, symmetric = TRUE)
      vals <- pmax(eg$values, .Machine$double.eps)
      Vinvhalf <- diag(1 / sqrt(vals)) %*% t(eg$vectors)
    }
    D_list[[i]] <- Vinvhalf %*% Z  
  }
  D_stack <- do.call(rbind, D_list)
  
  out_rows <- lapply(seq_len(nrow(intervals)), function(j) {
    iv  <- intervals[j, ]
    res <- find_single_cp(D_stack, iv$st, iv$ed, trim)
    data.frame(st = iv$st,
               ed = iv$ed,
               est.cp = res$est.cp,
               val = res$val,
               trim = trim)
  })
  
  cands <- do.call(rbind, out_rows)
  rownames(cands) <- NULL
  cands
}

# With NOT Principle
TNotSBS <- function(G, G_dim,
                    trim = round(2 * log(dim(G)[length(G_dim) + 1])),
                    method = c("fixed", "oracle"),
                    threshold = NULL,
                    m = NULL,
                    V.diag = TRUE,
                    lrv = TRUE,
                    n = NULL,
                    lbd = NULL,
                    single = FALSE) {
  
  method <- match.arg(method)
  
  cands <- cand_sbs(G, G_dim,
                    trim = trim,
                    V.diag = V.diag,
                    lrv = lrv,
                    n = n,
                    lbd = lbd,
                    single = single)
  
  cands <- cands[order(cands$ed - cands$st), ]
  st_vec <- cands$st
  ed_vec <- cands$ed
  cp_vec <- cands$est.cp
  val_vec <- cands$val
  
  if (method == "fixed") {
    
    if (is.null(threshold))
      stop("Must supply 'threshold' in fixed mode.")
    
    avail <- which(val_vec >= threshold + 5e-3)
    sel <- integer(0)
    
    while (length(avail) > 0) {
      k <- avail[1]
      sel <- c(sel, k)
      
      avail <- avail[-1]
      overlap <- st_vec[avail] < cp_vec[k] & ed_vec[avail] >= cp_vec[k]
      if (any(overlap)) avail <- avail[!overlap]
    }
    
    res <- cands[sel, , drop = FALSE]
    attr(res, "threshold") <- threshold
    
  } else {  
    ## oracle
    if (is.null(m))
      stop("Must supply 'm' in oracle method.")
    
    thds <- sort(unique(val_vec), decreasing = TRUE)
    sel <- NULL; chosen_th <- NA
    
    for (thd in thds) {
      avail <- which(val_vec >= thd)
      sel_tmp <- integer(0)
      
      while (length(sel_tmp) < m && length(avail) > 0) {
        k <- avail[1]
        sel_tmp <- c(sel_tmp, k)
        
        avail <- avail[-1]
        overlap <- st_vec[avail] < cp_vec[k] & ed_vec[avail] >= cp_vec[k]
        if (any(overlap)) avail <- avail[!overlap]
      }
      
      if (length(sel_tmp) == m) { sel <- sel_tmp; chosen_th <- thd; break }
    }
    
    if (is.null(sel))
      stop(sprintf("Cannot find threshold yielding m = %d", m))
    
    res <- cands[sel, , drop = FALSE]
    res <- res[order(res$est.cp), , drop = FALSE]
    
    attr(res, "selected_threshold")  <- chosen_th
    attr(res, "no_change_threshold") <- max(val_vec)
  }
  
  rownames(res) <- NULL
  res
}
