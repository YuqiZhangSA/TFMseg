# lower-triangular vectorisation
vech <- function(M) M[lower.tri(M, diag = TRUE)]

# mode unfolding
unfold_mode_k <- function(tensor, mode) {
  dims <- dim(tensor)
  perm <- c(mode, setdiff(seq_along(dims), mode))
  unfolded <- aperm(tensor, perm)
  matrix(unfolded, nrow = dims[mode])
}

# "next-highest" threshold: largest cusum value among non-overlapping intervals
compute_next_threshold <- function(current_threshold, S, thds, all_results, Time) {
  cand_candidates <- sort(thds[thds < current_threshold], decreasing = TRUE)
  
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
    thd = NA_real_, st = NA_integer_, ed = NA_integer_,
    cp  = NA_integer_, len = NA_integer_, a = NA_real_
  ))
}




# find single cp on given interval
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
cand_sbs <- function(G, 
                     G_dim,
                     trim = NULL,
                     V.diag = TRUE,
                     lrv = TRUE,
                     n = NULL,
                     single = FALSE) {
  
  K <- length(G_dim) 
  Time <- dim(G)[K + 1]
  
  if (is.null(trim)) trim <- round(0.25 * dim(G)[length(G_dim) + 1] / log(dim(G)[length(G_dim) + 1]))
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
    
    ## long-run covariance HAC
    V <- GG_0 %*% t(GG_0) / Time
    if (lrv && n >= 1) {
      for (ell in seq_len(n)) {
        temp <- GG_0[, 1:(Time - ell), drop=FALSE] %*% t(GG_0[, (1:(Time - ell)) + ell, drop=FALSE]) / Time
        w <- 1 - ell / (n + 1)
        V <- V + w * (temp + t(temp))
      }
    }
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
  cands
}



# -------------------------------------------------------------------
# TFMseg
# -------------------------------------------------------------------
TFMseg <- function(G, G_dim,
                   trim = round(0.25 * dim(G)[length(G_dim) + 1] / log(dim(G)[length(G_dim) + 1])),
                   method = c("fixed", "oracle"),
                   threshold = NULL,
                   m = NULL,
                   V.diag = TRUE,
                   lrv = TRUE,
                   n = NULL,
                   single = FALSE) {
  
  method <- match.arg(method)
  
  cands <- cand_sbs(G, G_dim,
                    trim = trim,
                    V.diag = V.diag,
                    lrv = lrv,
                    n = n,
                    single = single)
  
  # sort by interval length (NOT)
  cands <- cands[order(cands$ed - cands$st), ]
  st_vec <- cands$st
  ed_vec <- cands$ed
  cp_vec <- cands$est.cp
  val_vec <- cands$val
  
  Time <- dim(G)[length(G_dim) + 1]
  

  length_sbs <- ed_vec - st_vec
  bound_vec <- pmin(ed_vec - cp_vec, cp_vec - st_vec) 
  
  is_full <- (st_vec == 0 & ed_vec == Time) | (st_vec == 1 & ed_vec == Time)
  
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
  null_max_excl <- if (length(val_excl) > 0) max(val_excl, na.rm = TRUE) else NA_real_
  
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
    return(res)
  }
  
  # oracle
  if (is.null(m))
    stop("Must supply 'm' in oracle method.")
  
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
      overlap <- st_vec[avail] < cp_vec[k] & ed_vec[avail] >= cp_vec[k]
      if (any(overlap)) avail <- avail[!overlap]
    }
    
    if (length(sel_tmp) == m) {
      sel <- sel_tmp
      chosen_th <- thd
      break
    }
  }
  
  if (is.null(sel))
    stop(sprintf("Cannot find threshold yielding m = %d", m))
  
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
