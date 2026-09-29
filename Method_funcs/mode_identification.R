# covariance for specific unfolding
mode_cov_from_G <- function(G, mode_k, a, b) {
  dims <- dim(G)
  KK <- length(dims) - 1L
  rk <- dims[mode_k]
  C <- matrix(0, rk, rk)
  
  for (t in (a + 1L):b) {
    idx <- c(rep(list(TRUE), KK), list(t))
    Gt <- do.call(`[`, c(list(G), idx, list(drop = FALSE)))
    Gt <- array(Gt, dim = dims[1:KK])
    perm <- c(mode_k, setdiff(seq_len(KK), mode_k))
    Mk <- matrix(aperm(Gt, perm), nrow = rk)
    C <- C + Mk %*% t(Mk)
  }
  
  C / (b - a)
}


# eq(3.5) in the main text
zeta_pair <- function(G_pre, G_post) {
  s0 <- sum(diag(G_pre))
  s1 <- sum(diag(G_post))
  
  if (!is.finite(s0) || !is.finite(s1) || s0 <= 0 || s1 <= 0) {
    return(NA_real_)
  }
  
  Xi <- (G_post / s1) - (G_pre / s0)
  norm(Xi, type = "2")
}


id_modes <- function(G,
                    G_dim,
                    dim_obs,
                    cp_vec,
                    st = 1L,
                    ed = dim(G)[length(G_dim) + 1L],
                    id_thd_coef = NULL
) {
  stopifnot(is.array(G), length(dim(G)) == length(G_dim) + 1L)
  stopifnot(all(as.integer(dim(G))[seq_along(G_dim)] == as.integer(G_dim)))
  
  if(is.null(id_thd_coef)) id_thd_coef <- 3.8
  
  dim_obs <- as.integer(dim_obs)
  p <- prod(dim_obs)
  K <- length(G_dim)
  Ttot <- dim(G)[K + 1L]
  
  st <- as.integer(st)
  ed <- as.integer(ed)
  if (st < 1L || ed > Ttot || st >= ed) stop("Invalid st/ed.")
  
  cp_vec <- sort(unique(as.integer(cp_vec)))
  cp_vec <- cp_vec[cp_vec > st & cp_vec < ed]
  if (!length(cp_vec)) stop("No interior cps.")
  
  delta_Tp <- (1 / sqrt(Ttot)) + (1 / p)
  
  bnds <- c(st - 1L, cp_vec, ed)
  qhat <- length(cp_vec)
  
  zeta_mat <- matrix(
    NA_real_,
    nrow = K,
    ncol = qhat,
    dimnames = list(paste0("mode", seq_len(K)), paste0("cp", seq_len(qhat)))
  )
  ratio_Tp <- zeta_mat
  changed <- matrix(FALSE, nrow = K, ncol = qhat, dimnames = dimnames(zeta_mat))
  
  for (j in seq_len(qhat)) {
    a_pre <- bnds[j]
    b_pre <- bnds[j + 1L]
    a_post <- bnds[j + 1L]
    b_post <- bnds[j + 2L]
    
    for (k in seq_len(K)) {
      G_pre <- mode_cov_from_G(G, mode_k = k, a = a_pre, b = b_pre)
      G_post <- mode_cov_from_G(G, mode_k = k, a = a_post, b = b_post)
      
      z <- zeta_pair(G_pre, G_post)
      
      zeta_mat[k, j] <- z
      ratio_Tp[k, j] <- z / delta_Tp
      changed[k, j] <- is.finite(ratio_Tp[k, j]) && (ratio_Tp[k, j] > id_thd_coef)
    }
  }
  
  list(
    zeta_mat = zeta_mat,
    ratio_Tp = ratio_Tp,
    changed = changed,
    delta_Tp = delta_Tp,
    cp_vec = cp_vec
  )
}



# evaluator for one setting
eval_modes <- function(
    zeta_by_cp,
    dim_obs = c(10, 10, 10),
    id_thd_coef,
    true_mode_list
) {
  q <- length(true_mode_list)
  K <- length(dim_obs)
  
  if (length(zeta_by_cp) != q) {
    stop("zeta_by_cp and true_mode_list must have same length.")
  }
  
  sims_list <- lapply(zeta_by_cp, function(df) {
    if (is.data.frame(df) && nrow(df)) sort(unique(as.integer(df$sim))) else integer(0)
  })
  sims <- sims_list[[1]]
  if (q > 1L) {
    for (jj in 2:q) sims <- intersect(sims, sims_list[[jj]])
  }
  sims <- sort(unique(as.integer(sims)))
  if (!length(sims)) stop("No common sims found across cp tables.")
  
  eval_rows <- list()
  ie <- 0L
  
  for (sim_id in sims) {
    for (j in seq_len(q)) {
      dfj <- zeta_by_cp[[j]]
      dfjs <- dfj[as.integer(dfj$sim) == sim_id, , drop = FALSE]
      
      if (!nrow(dfjs)) {
        ie <- ie + 1L
        eval_rows[[ie]] <- data.frame(
          sim = sim_id, cp_index = j,
          TPR = NA_real_, FPR = NA_real_, J = NA_real_
        )
        next
      }
      
      if (!("ratio_Tp" %in% names(dfjs))) {
        stop("Missing column 'ratio_Tp' in zeta_by_cp[[", j, "]].")
      }
      
      pred_set <- sort(unique(dfjs$mode[
        is.finite(dfjs$ratio_Tp) & (dfjs$ratio_Tp > id_thd_coef)
      ]))
      
      true_set <- sort(unique(true_mode_list[[j]]))
      true_set <- true_set[true_set >= 1L & true_set <= K]
      
      tp <- length(intersect(true_set, pred_set))
      fp <- length(setdiff(pred_set, true_set))
      
      tpr <- if (length(true_set)) tp / length(true_set) else NA_real_
      fpr <- if ((K - length(true_set)) > 0L) fp / (K - length(true_set)) else NA_real_
      J <- if (is.finite(tpr) && is.finite(fpr)) tpr - fpr else NA_real_
      
      ie <- ie + 1L
      eval_rows[[ie]] <- data.frame(
        sim = sim_id, cp_index = j,
        TPR = tpr, FPR = fpr, J = J
      )
    }
  }
  
  eval_per_rep <- do.call(rbind, eval_rows)
  
  eval_mean <- aggregate(
    cbind(TPR, FPR, J) ~ cp_index,
    data = eval_per_rep,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  n_valid <- aggregate(!is.na(TPR) ~ cp_index, data = eval_per_rep, FUN = sum)
  names(n_valid)[2] <- "n_valid_reps"
  eval_mean <- merge(eval_mean, n_valid, by = "cp_index", all.x = TRUE, sort = TRUE)
  
  list(
    eval_per_rep = eval_per_rep,
    eval_mean = eval_mean
  )
}


# wrapper over many settings
eval_all <- function(
    res_all_true,
    id_thd_coef = 3,
    true_mode_list
) {
  stopifnot(is.list(res_all_true), length(res_all_true) > 0)
  
  out <- list()
  ii <- 0L
  
  for (store in res_all_true) {
    Time <- store$meta$Time
    dobs <- store$meta$dim_obs
    dobs_chr <- paste(dobs, collapse = "x")
    
    ev <- eval_modes(
      zeta_by_cp = store$zeta_by_cp,
      dim_obs = dobs,
      id_thd_coef = id_thd_coef,
      true_mode_list = true_mode_list
    )
    
    tab <- ev$eval_mean
    tab$Time <- Time
    tab$dim_obs <- dobs_chr
    
    ii <- ii + 1L
    out[[ii]] <- tab
  }
  
  do.call(rbind, out)
}