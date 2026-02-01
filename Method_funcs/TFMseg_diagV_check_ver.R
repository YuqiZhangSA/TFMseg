# lower-triangular vectorisation
vech <- function(M) M[lower.tri(M, diag = TRUE)]

find_single_cp_t <- function(D, st, ed, trim) {
  Time <- ncol(D)
  st_i <- max(1, as.integer(st))
  ed_i <- min(Time, as.integer(ed))
  len <- ed_i - st_i + 1
  stopifnot(len > 2 * trim + 1)
  
  stat <- numeric(len)
  total <- rowSums(D[, st_i:ed_i, drop = FALSE])
  lsum <- numeric(nrow(D))
  
  for (k in seq_len(len - 1)) {
    lsum <- lsum + D[, st_i + k - 1]
    rsum <- total - lsum
    nL <- k
    nR <- len - k
    diff <- rsum / nR - lsum / nL
    stat[k] <- sqrt((nL * nR) / len) * sqrt(sum(diff^2))
  }
  
  stat[seq_len(trim)] <- 0
  stat[(len - trim + 1):len] <- 0
  
  best_rel <- which.max(stat)
  best_abs_i <- st_i + best_rel - 1
  
  data.frame(
    est.cp = best_abs_i,
    val = stat[best_rel],
    st = st,
    ed = ed,
    trim = trim
  )
}


# candidates under SBS, using stacked cusum
cand_sbs <- function(G, G_dim,
                     trim = NULL,
                     V.diag = TRUE,
                     lrv = TRUE,
                     n = NULL,
                     #lbd = NULL,
                     single = FALSE,
                     return_diagV = FALSE) {
  
  K <- length(G_dim) # how many mode
  Time <- dim(G)[K + 1]
  
  # To check
  diagV_list <- vector("list", K)
  
  if (is.null(trim)) trim <- 3 * round(log(Time))
  trim <- as.integer(trim);  stopifnot(trim >= 0)
  if (is.null(n)) n <- floor(Time^0.25)
  
  if (single) {
    intervals <- data.frame(st = 0, ed = Time)
  } else {
    intervals <- seeded_intervals(Time = Time, minl = floor(0.5*Time/log(Time)))
  }
  
  # for each unfolding
  D_list <- vector("list", K)
  for (i in seq_len(K)) {
    r_i <- G_dim[i]
    d_i <- r_i * (r_i + 1) / 2
    GG <- matrix(0, d_i, Time)
    
    for (t in seq_len(Time)) {
      if (K == 1) {
        g_t <- G[, t]
        S <- g_t %o% g_t
      } else {
        G_pure <- slice_time(G, t)
        G_mi <- unfold_mode_k(G_pure, i) #r_k\times\rmk (G minus i)
        S <- G_mi %*% t(G_mi)
      }
      GG[, t] <- vech(S)
    }
    
    #GG_0 <- GG - rowMeans(GG)
    mean_vec <- rowMeans(GG)
    GG_0 <- sweep(GG, 1, mean_vec, "-")
    
    ## long-run variance HAC
    V <- GG_0 %*% t(GG_0) / Time
    if (lrv && n >= 1) {
      for (ell in seq_len(n)) {
        temp <- GG_0[, 1:(Time - ell), drop=FALSE] %*% t(GG_0[, (1:(Time - ell)) + ell, drop=FALSE]) / Time
        w <- 1 - ell / (n + 1)
        V <- V + w * (temp + t(temp))
      }
    }
    
    # TO check
    diagV <- if (is.matrix(V)) diag(V) else as.numeric(V)
    diagV[!is.finite(diagV)] <- NA_real_
    if (all(is.na(diagV))) warning("All diag(V) are NA - check V construction.")
    diagV_list[[i]] <- diagV
    
    
    if (V.diag) {
      if (dim(V)[1] == 1) {
        d <- V; Vinvhalf <- as.matrix(1 / sqrt(d))
      } else {
        d <- diag(V); d[d <= 0] <- min(d[d > 0])
        Vinvhalf <- diag(1 / sqrt(d))
      }
    } else {
      if (dim(V)[1] == 1) {
        d <- V; Vinvhalf <- as.matrix(1 / sqrt(d))
      } else {
        eg <- eigen(V, symmetric = TRUE)
        vals <- pmax(eg$values, .Machine$double.eps)
        Vinvhalf <- diag(1 / sqrt(vals)) %*% t(eg$vectors)
      }
    }
    D_list[[i]] <- Vinvhalf %*% GG_0  
  }
  D_stack <- do.call(rbind, D_list)
  
  out_rows <- lapply(seq_len(nrow(intervals)), function(j) {
    iv  <- intervals[j, ]
    res <- find_single_cp_t(D_stack, iv$st, iv$ed, trim)
    data.frame(st = res$st,
               ed = res$ed,
               est.cp = res$est.cp,
               val = res$val,
               trim = res$trim)
  })
  
  cands <- do.call(rbind, out_rows)
  rownames(cands) <- NULL
  #cands
  if (return_diagV) {
    return(list(cands = cands, diagV_list = diagV_list))
  } else {
    return(cands)
  }
}





unfold_mode_k <- function(tensor, mode) {
  dims <- dim(tensor)
  perm <- c(mode, setdiff(seq_along(dims), mode))
  unfolded <- aperm(tensor, perm)
  matrix(unfolded, nrow = dims[mode])
}

compute_next_threshold <- function(current_threshold, S, thds, all_results, Time) {
  cand_candidates <- sort(thds[thds < current_threshold], decreasing = TRUE)
  
  for (cand in cand_candidates) {
    candidate_intervals <- all_results[all_results$val == cand, , drop = FALSE]
    if (nrow(candidate_intervals) == 0) next
    
    for (j in seq_len(nrow(candidate_intervals))) {
      st <- candidate_intervals$st[j]
      ed <- candidate_intervals$ed[j]
      if (all(S - st < log(Time) | ed - S < log(Time))) {
        return(cand)
      }
    }
  }
  return(NA_real_)
}







TFMseg <- function(G, G_dim,
                   trim = round(0.25 * dim(G)[length(G_dim) + 1] / log(dim(G)[length(G_dim) + 1])),
                   method = c("fixed", "oracle"),
                   threshold = NULL,
                   m = NULL,
                   V.diag = TRUE,
                   lrv = TRUE,
                   n = NULL,
                   single = FALSE,
                   return_diagV = FALSE) {
  
  method <- match.arg(method)
  
  tmp <- cand_sbs(G, G_dim,
                  trim = trim,
                  V.diag = V.diag,
                  lrv = lrv,
                  n = n,
                  single = single,
                  return_diagV = return_diagV)
  
  if (return_diagV) {
    cands <- tmp$cands
    diagV_list <- tmp$diagV_list
  } else {
    cands <- tmp
    diagV_list <- NULL
  }
  
  cands <- cands[order(cands$ed - cands$st), ]
  st_vec <- cands$st
  ed_vec <- cands$ed
  cp_vec <- cands$est.cp
  val_vec <- cands$val
  
  if (method == "fixed") {
    if (is.null(threshold)) stop("Must supply 'threshold' in fixed mode.")
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
    res <- res[order(res$est.cp), , drop = FALSE]
    rownames(res) <- NULL
    
  } else {
    if (is.null(m)) stop("Must supply 'm' in oracle method.")
    thds <- sort(unique(val_vec), decreasing = TRUE)
    sel <- NULL
    chosen_th <- NA
    
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
      if (length(sel_tmp) == m) {
        sel <- sel_tmp
        chosen_th <- thd
        break
      }
    }
    if (is.null(sel)) stop(sprintf("Cannot find threshold yielding m = %d", m))
    
    res <- cands[sel, , drop = FALSE]
    res <- res[order(res$est.cp), , drop = FALSE]
    res$selected_threshold  <- chosen_th
    res$no_change_threshold <- max(val_vec, na.rm = TRUE)
    
    Time <- dim(G)[length(G_dim) + 1]
    S <- res$est.cp
    all_results <- cands[, c("st","ed","val")]
    res$next_highest_threshold <- compute_next_threshold(
      current_threshold = chosen_th,
      S = S,
      thds = thds,
      all_results = all_results,
      Time = Time
    )
    rownames(res) <- NULL
  }
  
  if (return_diagV) return(list(res = res, diagV_list = diagV_list))
  res
}

