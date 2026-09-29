vech <- function(M) M[lower.tri(M, diag = TRUE)]


unfold_mode_k <- function(tensor, mode) {
  dims <- dim(tensor)
  perm <- c(mode, setdiff(seq_along(dims), mode))
  unfolded <- aperm(tensor, perm)
  matrix(unfolded, nrow = dims[mode])
}


# interval-wise HAC/LRV standardisation
int_lrv <- function(Z, V.diag = TRUE, lrv = TRUE,
                    n = NULL, eps = 1e-8, centre = TRUE) {
  
  d <- nrow(Z)
  L <- ncol(Z)
  
  if (centre) {
    Zc <- sweep(Z, 1, rowMeans(Z), "-")
  } else {
    Zc <- Z
  }
  
  if (is.null(n)) {
    n_use <- floor(L^0.25)
  } else {
    n_use <- min(n, L - 1)
  }
  
  n_use <- max(0, n_use)
  
  V <- Zc %*% t(Zc) / L
  
  if (lrv && n_use >= 1) {
    for (ell in seq_len(n_use)) {
      Z1 <- Zc[, 1:(L - ell), drop = FALSE]
      Z2 <- Zc[, (1 + ell):L, drop = FALSE]
      
      tmp <- Z1 %*% t(Z2) / L
      w <- 1 - ell / (n_use + 1)
      
      V <- V + w * (tmp + t(tmp))
    }
  }
  
  V <- (V + t(V)) / 2
  
  flag <- FALSE
  
  if (V.diag) {
    
    dd <- diag(V)
    
    pos <- dd[is.finite(dd) & dd > eps]
    floor_val <- if (length(pos) > 0) min(pos) else eps
    
    dd_new <- dd
    bad <- !is.finite(dd_new) | dd_new <= floor_val
    
    if (any(bad)) {
      dd_new[bad] <- floor_val
      flag <- TRUE
    }
    
    Vhalf_inv <- diag(1 / sqrt(dd_new), nrow = d, ncol = d)
    
  } else {
    
    eig <- eigen(V, symmetric = TRUE)
    dd <- eig$values
    
    pos <- dd[is.finite(dd) & dd > eps]
    floor_val <- if (length(pos) > 0) min(pos) else eps
    
    dd_new <- dd
    bad <- !is.finite(dd_new) | dd_new <= floor_val
    
    if (any(bad)) {
      dd_new[bad] <- floor_val
      flag <- TRUE
    }
    
    Vhalf_inv <- diag(1 / sqrt(dd_new), nrow = d, ncol = d) %*% t(eig$vectors)
  }
  
  list(
    Vhalf_inv = Vhalf_inv,
    flag = flag,
    bandwidth = n_use
  )
}


make_GG_list <- function(G, G_dim) {
  
  K <- length(G_dim)
  Time <- dim(G)[K + 1]
  
  GG_list <- vector("list", K)
  
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
        G_mi <- unfold_mode_k(G_pure, i)
        S <- G_mi %*% t(G_mi)
      }
      
      GG[, t] <- vech(S)
    }
    
    GG_list[[i]] <- GG
  }
  
  GG_list
}



# stacked interval-wise standardisation over all tensor modes
int_std_stack <- function(GG_list,
                          st,
                          ed,
                          Time,
                          V.diag = TRUE,
                          lrv = TRUE,
                          n = NULL) {
  
  st_i <- max(1, as.integer(st))
  ed_i <- min(Time, as.integer(ed))
  idx <- st_i:ed_i
  len_i <- length(idx)
  
  D_list <- lapply(GG_list, function(GG) {
    
    Z <- GG[, idx, drop = FALSE]
    
    lrv_obj <- int_lrv(
      Z = Z,
      V.diag = V.diag,
      lrv = lrv,
      n = if (is.null(n)) NULL else min(n, len_i - 1),
      centre = TRUE
    )
    
    lrv_obj$Vhalf_inv %*% Z
  })
  
  do.call(rbind, D_list)
}



# find single cp on a given interval
find_single_cp_local_lrv <- function(D, st, ed, trim) {
  
  len <- ncol(D)
  stopifnot(len >= 2 * trim + 1)
  
  st_i <- max(1, as.integer(st))
  
  stat <- numeric(len)
  total <- rowSums(D)
  lsum <- numeric(nrow(D))
  
  for (k in seq_len(len - 1)) {
    lsum <- lsum + D[, k]
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



# next-highest threshold: largest cusum value among non-overlapping intervals
compute_next_threshold <- function(current_threshold, S, thds, all_results, Time) {
  
  cand_candidates <- sort(
    unique(thds[thds < current_threshold]),
    decreasing = TRUE
  )
  
  if (length(S) == 0 || length(cand_candidates) == 0) {
    return(list(
      thd = NA_real_,
      st = NA_integer_,
      ed = NA_integer_,
      cp = NA_integer_,
      len = NA_integer_,
      a = NA_real_
    ))
  }
  
  for (cand in cand_candidates) {
    candidate_intervals <- all_results[all_results$val == cand, , drop = FALSE]
    
    if (nrow(candidate_intervals) == 0) next
    
    for (j in seq_len(nrow(candidate_intervals))) {
      st <- candidate_intervals$st[j]
      ed <- candidate_intervals$ed[j]
      cp <- candidate_intervals$est.cp[j]
      
      if (all(S - st < log(Time) | ed - S < log(Time))) {
        return(list(
          thd = cand,
          st = st,
          ed = ed,
          cp = cp,
          len = ed - st,
          a = min(ed - cp, cp - st)
        ))
      }
    }
  }
  
  return(list(
    thd = NA_real_,
    st = NA_integer_,
    ed = NA_integer_,
    cp = NA_integer_,
    len = NA_integer_,
    a = NA_real_
  ))
}



# candidates under SBS, using stacked interval-wise HAC standardisation
cand_sbs_local_lrv <- function(G, 
                     G_dim,
                     trim = NULL,
                     V.diag = TRUE,
                     lrv = TRUE,
                     n = NULL,
                     single = FALSE) {
  
  K <- length(G_dim) 
  Time <- dim(G)[K + 1]
  
  if (is.null(trim)) {
    trim <- floor(0.25 * Time / log(Time))
  }
  
  trim <- as.integer(trim)
  stopifnot(trim >= 0)
  
  if (is.null(n)) {
    n <- floor(Time^0.15)
  }
  
  if (single) {
    intervals <- data.frame(st = 0, ed = Time)
  } else {
    intervals <- seeded_intervals(
      Time = Time,
      minl = 2 * trim + 1
    )
  }
  
  GG_list <- make_GG_list(G, G_dim)
  
  out_rows <- lapply(seq_len(nrow(intervals)), function(j) {
    
    iv <- intervals[j, ]
    
    st_i <- max(1, as.integer(iv$st))
    ed_i <- min(Time, as.integer(iv$ed))
    len_i <- ed_i - st_i + 1
    
    if (len_i < 2 * trim + 1) {
      return(NULL)
    }
    
    D_stack <- int_std_stack(
      GG_list = GG_list,
      st = iv$st,
      ed = iv$ed,
      Time = Time,
      V.diag = V.diag,
      lrv = lrv,
      n = n
    )
    
    res <- find_single_cp_local_lrv(
      D = D_stack,
      st = iv$st,
      ed = iv$ed,
      trim = trim
    )
    
    data.frame(
      st = res$st,
      ed = res$ed,
      est.cp = res$est.cp,
      val = res$val,
      trim = res$trim
    )
  })
  
  out_rows <- out_rows[!vapply(out_rows, is.null, logical(1))]
  
  if (length(out_rows) == 0) {
    return(data.frame(
      st = integer(0),
      ed = integer(0),
      est.cp = integer(0),
      val = numeric(0),
      trim = integer(0)
    ))
  }
  
  cands <- do.call(rbind, out_rows)
  rownames(cands) <- NULL
  
  cands
}


# TFMseg under local HAC
TFMseg_local_lrv <- function(G, G_dim,
                   trim = floor(0.25 * dim(G)[length(G_dim) + 1] /
                                  log(dim(G)[length(G_dim) + 1])),
                   method = c("fixed", "oracle"),
                   threshold = NULL,
                   m = NULL,
                   V.diag = TRUE,
                   lrv = TRUE,
                   n = NULL,
                   single = FALSE,
                   tie_breaking = TRUE) {
  
  method <- match.arg(method)
  
  cands <- cand_sbs_local_lrv(
    G = G,
    G_dim = G_dim,
    trim = trim,
    V.diag = V.diag,
    lrv = lrv,
    n = n,
    single = single
  )
  
  Time <- dim(G)[length(G_dim) + 1]
  
  if (nrow(cands) == 0) {
    res <- data.frame(
      est.cp = NA_integer_,
      val = NA_real_,
      st = NA_integer_,
      ed = NA_integer_,
      trim = as.integer(trim),
      stringsAsFactors = FALSE
    )
    rownames(res) <- NULL
    return(res)
  }
  
  # sort by interval length, i.e. NOT rule
  if(tie_breaking){
    cands <- cands[order(cands$ed - cands$st, -cands$val, cands$st, cands$ed), , drop = FALSE]
  } else {
    cands <- cands[order(cands$ed - cands$st), ]
  }

  st_vec <- cands$st
  ed_vec <- cands$ed
  cp_vec <- cands$est.cp
  val_vec <- cands$val
  
  length_sbs <- ed_vec - st_vec
  bound_vec <- pmin(ed_vec - cp_vec, cp_vec - st_vec)
  
  is_full <- (st_vec == 0 & ed_vec == Time) | 
    (st_vec == 1 & ed_vec == Time)
  
  idx_max <- which.max(val_vec)
  
  null_max_meta <- list(
    null_max = val_vec[idx_max],
    null_max_st = st_vec[idx_max],
    null_max_ed = ed_vec[idx_max],
    null_max_cp = cp_vec[idx_max],
    null_max_len = length_sbs[idx_max],
    null_max_bound = bound_vec[idx_max]
  )
  
  val_excl <- val_vec[!is_full]
  null_max_excl <- if (length(val_excl) > 0) {
    max(val_excl, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  attach_null_diagnostics <- function(res) {
    res$null_max <- null_max_meta$null_max
    res$null_max_excl <- null_max_excl
    res$null_max_len <- null_max_meta$null_max_len
    res$null_max_bound <- null_max_meta$null_max_bound
    res$null_max_st <- null_max_meta$null_max_st
    res$null_max_ed <- null_max_meta$null_max_ed
    res$null_max_cp <- null_max_meta$null_max_cp
    res
  }
  
  if (method == "fixed") {
    
    if (is.null(threshold)) {
      dr <- sum(G_dim * (G_dim + 1L) / 2L)
      threshold <- 747.0283085 * sqrt(log(Time)) +
                    0.9132869 * sqrt(dr) +
                    5538.8887436 * sqrt(1 / log(Time)) +
                    5329.8863270 * log(log(Time)) / sqrt(log(Time)) -
                    7984.4027860
      }
    
    
    avail <- which(val_vec >= threshold + 5e-3)
    sel <- integer(0)
    
    while (length(avail) > 0) {
      k <- avail[1]
      sel <- c(sel, k)
      avail <- avail[-1]
      
      overlap <- st_vec[avail] < cp_vec[k] & 
        ed_vec[avail] >= cp_vec[k]
      
      if (any(overlap)) {
        avail <- avail[!overlap]
      }
    }
    
    if (length(sel) == 0) {
      res <- data.frame(
        est.cp = NA_integer_,
        val = NA_real_,
        st = NA_integer_,
        ed = NA_integer_,
        trim = as.integer(trim),
        stringsAsFactors = FALSE
      )
    } else {
      res <- cands[sel, , drop = FALSE]
      res <- res[order(res$est.cp), , drop = FALSE]
      rownames(res) <- NULL
    }
    
    res <- attach_null_diagnostics(res)
    rownames(res) <- NULL
    return(res)
  }
  
  # oracle mode
  if (is.null(m)) {
    stop("Must supply 'm' in oracle method.")
  }
  
  thds <- sort(unique(val_vec), decreasing = TRUE)
  
  sel <- NULL
  chosen_th <- NA_real_
  
  for (thd in thds) {
    
    avail <- which(val_vec >= thd)
    sel_tmp <- integer(0)
    
    while (length(sel_tmp) < m && length(avail) > 0) {
      k <- avail[1]
      sel_tmp <- c(sel_tmp, k)
      avail <- avail[-1]
      
      overlap <- st_vec[avail] < cp_vec[k] & 
        ed_vec[avail] >= cp_vec[k]
      
      if (any(overlap)) {
        avail <- avail[!overlap]
      }
    }
    
    if (length(sel_tmp) == m) {
      sel <- sel_tmp
      chosen_th <- thd
      break
    }
  }
  
  if (is.null(sel)) {
    stop(sprintf("Cannot find threshold yielding m = %d", m))
  }
  
  res <- cands[sel, , drop = FALSE]
  res <- res[order(res$est.cp), , drop = FALSE]
  
  res$selected_threshold <- chosen_th
  
  res <- attach_null_diagnostics(res)
  
  S <- res$est.cp
  
  all_results <- cands[, c("st", "ed", "est.cp", "val")]
  
  next_obj <- compute_next_threshold(
    current_threshold = chosen_th,
    S = S,
    thds = thds,
    all_results = all_results,
    Time = Time
  )
  
  res$next_highest_threshold <- next_obj$thd
  res$next_len <- next_obj$len
  res$next_bound <- next_obj$a
  res$next_st <- next_obj$st
  res$next_ed <- next_obj$ed
  res$next_cp <- next_obj$cp
  
  rownames(res) <- NULL
  
  res
}