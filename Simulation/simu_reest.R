mode_informed_proj <- function(X, dim_X, cp,
                               k,
                               Lambda1_pre, Lambda1_post,
                               Lambda_comp_pool,
                               r_k, centre = TRUE) {
  stopifnot(length(dim_X) == 3L)
  stopifnot(k %in% c(2L, 3L))
  
  Time <- dim(X)[4]
  stopifnot(!is.na(Time), Time >= 2L)
  stopifnot(cp >= 1L, cp < Time)
  
  if (centre) {
    X_use <- centre_along_axis(X, axis = length(dim(X)))
  } else {
    X_use <- X
  }
  
  p_k <- dim_X[k]
  
  if (k == 2L) {
    p_mk <- dim_X[1] * dim_X[3]
    Lambda_kron_pre  <- kronecker_list(list(Lambda_comp_pool, Lambda1_pre))
    Lambda_kron_post <- kronecker_list(list(Lambda_comp_pool, Lambda1_post))
  } else {
    p_mk <- dim_X[1] * dim_X[2]
    Lambda_kron_pre  <- kronecker_list(list(Lambda_comp_pool, Lambda1_pre))
    Lambda_kron_post <- kronecker_list(list(Lambda_comp_pool, Lambda1_post))
  }
  
  perm <- c(k, setdiff(seq_len(3L), k), 4L)
  Xp <- aperm(X_use, perm)
  
  sum_Y <- matrix(0, nrow = p_k, ncol = p_k)
  
  for (t in seq_len(cp)) {
    Xkt <- matrix(slice_time(Xp, t), nrow = p_k)
    Ykt <- (1 / p_mk) * Xkt %*% Lambda_kron_pre
    sum_Y <- sum_Y + Ykt %*% t(Ykt)
  }
  
  for (t in seq.int(cp + 1L, Time)) {
    Xkt <- matrix(slice_time(Xp, t), nrow = p_k)
    Ykt <- (1 / p_mk) * Xkt %*% Lambda_kron_post
    sum_Y <- sum_Y + Ykt %*% t(Ykt)
  }
  
  Gamma_Y <- sum_Y / (Time * p_k)
  eigY <- eigen(Gamma_Y, symmetric = TRUE)
  
  eigY$vectors[, seq_len(r_k), drop = FALSE] * sqrt(p_k)
}

compare_reest <- function(
    dim_obs,
    dim_latent,
    Time,
    dep = TRUE,
    theta_coef = 0.5,
    data_setting = c("s3_A2","s3_3I", "s0", "s3_A1"),
    r_hat_pre = NULL,
    r_hat_post = NULL,
    r_hat = NULL,
    detect_thd = NULL,
    V.diag = TRUE,
    lrv = TRUE,
    id_thd_coef = 3.5,
    seed = 901,
    use_true_cp = FALSE,
    N_exact = c("log", "mid")
) {
  stopifnot(length(dim_obs) == 3L, length(dim_latent) == 3L)
  stopifnot(all(dim_latent >= 3))
  data_setting <- match.arg(data_setting)
  N_exact <- match.arg(N_exact)
  
  theta_true <- floor(Time * theta_coef)
  if (theta_true <= 0 || theta_true >= Time) {
    stop("theta_true must lie in {1, ..., Time - 1}.")
  }
  
  set.seed(seed)
  
  if (data_setting == "s3_A2") {
    A <- matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 0
    ), nrow = 3, byrow = TRUE)
    true_changed <- c(TRUE, FALSE, FALSE)
  } else if (data_setting == "s3_3I") {
    A <- 3 * diag(3)
    true_changed <- c(TRUE, FALSE, FALSE)
  } else if (data_setting == "s0") {
    A <- diag(3)
    true_changed <- c(FALSE, FALSE, FALSE)
  } else if (data_setting == "s3_A1") {
    A <- matrix(0, 3, 3)
    A[1, 1] <- 0.5
    A[2, 1] <- rnorm(1, 0, 1); A[2, 2] <- 1
    A[3, 1] <- rnorm(1, 0, 1); A[3, 2] <- rnorm(1, 0, 1); A[3, 3] <- 1.5
    true_changed <- c(TRUE, FALSE, FALSE)
  }
  
  
  
  data_sim <- dgp_general(
    model = "tensor",
    dim_obs = dim_obs,
    dim_latent = dim_latent,
    Time = Time,
    dep = dep,
    coeff = 0.7,
    true_cp = theta_true,
    type_mode = list(1),
    type_change = list("l"),
    shift_ind = NULL,
    shift_mean = NULL,
    shift_var = NULL,
    transform_list = list(
      list(A, diag(dim_latent[2]), diag(dim_latent[3]))
    ),
    add_factors = NULL,
    add_factors_coeff = NULL,
    dist = "Gaussian"
  )
  
  X <- data_sim$X
  dim_X <- dim(X)[1:3]
  
  #X <- X[, , , Time:1, drop = FALSE]
  
  fit_full0 <- global_pca(
    X = X,
    dim_X = dim_X,
    r_hat = r_hat,
    centre = TRUE,
    proj = TRUE
  )
  
  G <- fit_full0$G_proj
  G_dim <- if (is.null(r_hat)) as.integer(fit_full0$r_hat) else as.integer(r_hat)
  
  dr <- sum(G_dim * (G_dim + 1L) / 2L)
  thd <- if (is.null(detect_thd)) {
    747.0283085 * sqrt(log(Time)) +
      0.9132869 * sqrt(dr) +
      5538.8887436 * sqrt(1 / log(Time)) +
      5329.8863270 * log(log(Time)) / sqrt(log(Time)) -
      7984.4027860
  } else {
    detect_thd
  }
  
  out_det <- TFMseg(
    G, G_dim,
    method = "fixed",
    threshold = thd,
    V.diag = V.diag,
    lrv = lrv
  )
  
  cp_hat_all <- out_det$est.cp %||% integer(0)
  cp_hat_all <- sort(unique(as.integer(cp_hat_all)))
  
  acc_win <- max(1L, as.integer(round(2 * log(Time))))
  
  if (isTRUE(use_true_cp)) {
    valid <- TRUE
    cp_used <- theta_true
    cp_hat_all <- if (data_setting == "s0") integer(0) else cp_hat_all
    
    cp_hat_store <- if (length(cp_hat_all) >= 1L) cp_hat_all[1] else NA_integer_
    cp_err_store <- if (!is.na(cp_hat_store)) cp_hat_store - theta_true else NA_integer_
  } else {
    if (N_exact == "log") {
      valid <- (length(cp_hat_all) == 1L) &&
        (abs(cp_hat_all[1] - theta_true) <= acc_win)
    } else if (N_exact == "mid") {
      left_bd  <- floor((0 + theta_true) / 2)
      right_bd <- floor((theta_true + Time) / 2)
      
      valid <- (length(cp_hat_all) == 1L) &&
        (cp_hat_all[1] > left_bd) &&
        (cp_hat_all[1] <= right_bd)
    }
    
    if (!valid) {
      return(list(
        valid = FALSE,
        theta_true = theta_true,
        cp_hat = NA_integer_,
        cp_used = NA_integer_,
        cp_hat_all = cp_hat_all,
        cp_err = NA_integer_,
        changed_hat = c(NA, NA, NA),
        mode_id_correct = c(NA, NA, NA),
        r_hat_full = as.integer(fit_full0$r_hat),
        r_hat_pre = rep(NA_integer_, 3L),
        r_hat_post = rep(NA_integer_, 3L),
        loss_mode = NULL
      ))
    }
    
    cp_used <- cp_hat_all[1]
    cp_hat_store <- cp_used
    cp_err_store <- cp_used - theta_true
  }
  
  id_out <- id_modes(
    G = G,
    G_dim = G_dim,
    dim_obs = dim_obs,
    cp_vec = cp_used,
    st = 1L,
    ed = Time,
    id_thd_coef = id_thd_coef
  )
  
  changed_hat <- as.logical(id_out$changed[, 1])
  mode_id_correct <- (changed_hat == true_changed)
  
  true_pre <- data_sim$loadings_hist[[1]]
  true_post <- data_sim$loadings_hist[[theta_true + 1L]]
  
  fit_pre <- global_pca(
    X = X, dim_X = dim_X, r_hat = r_hat_pre,
    st = 1L, ed = cp_used,
    centre = TRUE, proj = TRUE
  )
  
  fit_post <- global_pca(
    X = X, dim_X = dim_X, r_hat = r_hat_post,
    st = cp_used + 1L, ed = Time,
    centre = TRUE, proj = TRUE
  )
  
  r_full_used <- as.integer(fit_full0$r_hat)
  r_pre_used  <- as.integer(fit_pre$r_hat)
  r_post_used <- as.integer(fit_post$r_hat)
  
  Lambda_init_pool <- fit_full0$Lambda_init
  
  ## M2
  ## changed mode: same as M1, segment by segment
  Lambda_M2_pre <- vector("list", 3L)
  Lambda_M2_post <- vector("list", 3L)
  Lambda_M2_pool <- vector("list", 3L)
  
  Lambda_M2_pre[[1]]  <- fit_pre$Lambda_proj[[1]] #fit_pre$Lambda_proj[[1]]
  Lambda_M2_post[[1]] <- fit_post$Lambda_proj[[1]] 
  Lambda_M2_pool[[1]] <- NULL
  
  ## unchanged mode 2: one pooled estimator over the whole sample
  Lambda_M2_pool[[2]] <- mode_informed_proj(
    X = X, dim_X = dim_X, cp = cp_used, k = 2L,
    Lambda1_pre = fit_pre$Lambda_proj[[1]],
    Lambda1_post = fit_post$Lambda_proj[[1]],
    Lambda_comp_pool = Lambda_init_pool[[3]],
    r_k = r_full_used[2],
    centre = TRUE
  )
  
  Lambda_M2_pre[[2]]  <- NULL
  Lambda_M2_post[[2]] <- NULL
  
  ## unchanged mode 3: one pooled estimator over the whole sample
  Lambda_M2_pool[[3]] <- mode_informed_proj(
    X = X, dim_X = dim_X, cp = cp_used, k = 3L,
    Lambda1_pre = fit_pre$Lambda_proj[[1]],
    Lambda1_post = fit_post$Lambda_proj[[1]],
    Lambda_comp_pool = Lambda_init_pool[[2]],
    r_k = r_full_used[3],
    centre = TRUE
  )
  Lambda_M2_pre[[3]]  <- NULL
  Lambda_M2_post[[3]] <- NULL
  
  # Match the true loading matrix to the estimated rank by taking
  # the first r_hat true columns, consistent with the DGP ordering.
  true_match <- function(Ltrue, Lhat) {
    r_true <- ncol(Ltrue)
    r_hat  <- ncol(Lhat)
    if (r_hat > r_true) {
      stop("ncol(Lhat) exceeds ncol(Ltrue).")
    }
    Ltrue[, seq_len(r_hat), drop = FALSE]
  }
  
  loss_mode <- do.call(rbind, lapply(1:3, function(k) {
    ## M1
    true_pre_M1  <- true_match(true_pre[[k]],  fit_pre$Lambda_proj[[k]])
    true_post_M1 <- true_match(true_post[[k]], fit_post$Lambda_proj[[k]])
    
    M1_pre  <- tensorMiss::fle(fit_pre$Lambda_proj[[k]],  true_pre_M1)
    M1_post <- tensorMiss::fle(fit_post$Lambda_proj[[k]], true_post_M1)
    
    ## M2
    if (k == 1L) {
      true_pre_M2  <- true_match(true_pre[[k]],  Lambda_M2_pre[[k]])
      true_post_M2 <- true_match(true_post[[k]], Lambda_M2_post[[k]])
      
      M2_pre  <- tensorMiss::fle(Lambda_M2_pre[[k]],  true_pre_M2)
      M2_post <- tensorMiss::fle(Lambda_M2_post[[k]], true_post_M2)
      M2_pool <- NA_real_
    } else {
      true_pool_M2 <- true_match(true_pre[[k]], Lambda_M2_pool[[k]])
      M2_pool <- tensorMiss::fle(Lambda_M2_pool[[k]], true_pool_M2)
      
      M2_pre  <- NA_real_
      M2_post <- NA_real_
    }
    
    ## M3
    if (changed_hat[k]) {
      true_pre_M3  <- true_match(true_pre[[k]],  fit_pre$Lambda_proj[[k]])
      true_post_M3 <- true_match(true_post[[k]], fit_post$Lambda_proj[[k]])
      
      M3_pre  <- tensorMiss::fle(fit_pre$Lambda_proj[[k]],  true_pre_M3)
      M3_post <- tensorMiss::fle(fit_post$Lambda_proj[[k]], true_post_M3)
    } else {
      true_pre_M3  <- true_match(true_pre[[k]],  fit_full0$Lambda_proj[[k]])
      true_post_M3 <- true_match(true_post[[k]], fit_full0$Lambda_proj[[k]])
      
      M3_pre  <- tensorMiss::fle(fit_full0$Lambda_proj[[k]], true_pre_M3)
      M3_post <- tensorMiss::fle(fit_full0$Lambda_proj[[k]], true_post_M3)
    }
    
    data.frame(
      mode = k,
      true_changed = true_changed[k],
      changed_hat = changed_hat[k],
      id_correct = mode_id_correct[k],
      r_full = r_full_used[k],
      r_pre = r_pre_used[k],
      r_post = r_post_used[k],
      M1_pre = M1_pre,
      M1_post = M1_post,
      M2_pre = M2_pre,
      M2_post = M2_post,
      M2_pool = M2_pool,
      M3_pre = M3_pre,
      M3_post = M3_post
    )
  }))
  
  list(
    valid = TRUE,
    theta_true = theta_true,
    cp_hat = cp_hat_store,
    cp_used = cp_used,
    cp_hat_all = cp_hat_all,
    cp_err = cp_err_store,
    changed_hat = changed_hat,
    mode_id_correct = mode_id_correct,
    ratio_Tp = id_out$ratio_Tp[, 1],
    r_hat_full = r_full_used,
    r_hat_pre = r_pre_used,
    r_hat_post = r_post_used,
    loss_mode = loss_mode
  )
}




simu_reest <- function(
    nrep = 100,
    Time_list = NULL,
    dim_obs_list = NULL,
    dim_latent = c(3, 3, 3),
    dep = TRUE,
    theta_coef = 0.5,
    r_hat_pre = NULL,
    r_hat_post = NULL,
    r_hat = NULL,
    detect_thd = NULL,
    V.diag = TRUE,
    lrv = TRUE,
    id_thd_coef = 3.5,
    data_setting = c("s3_A2","s3_3I", "s0", "s3_A1"),
    use_true_cp = FALSE,
    N_exact = c("log", "mid")
) {
  if (is.null(Time_list)) stop("Please provide Time_list.")
  if (is.null(dim_obs_list)) stop("Please provide dim_obs_list.")
  data_setting <- match.arg(data_setting)
  N_exact <- match.arg(N_exact)
  
  settings <- expand.grid(
    Time = Time_list,
    idx_dim = seq_along(dim_obs_list),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  res_store <- vector("list", nrow(settings))
  
  for (i in seq_len(nrow(settings))) {
    Time_i <- settings$Time[i]
    dim_obs_i <- dim_obs_list[[settings$idx_dim[i]]]
    dim_chr <- paste(dim_obs_i, collapse = "x")
    
    message(sprintf(
      "[reest] %d/%d: Time=%d, dim_obs=%s, setting=%s",
      i, nrow(settings), Time_i, dim_chr, data_setting
    ))
    
    out_list <- vector("list", nrep)
    for (sim in seq_len(nrep)) {
      out_list[[sim]] <- compare_reest(
        dim_obs = dim_obs_i,
        dim_latent = dim_latent,
        Time = Time_i,
        dep = dep,
        theta_coef = theta_coef,
        data_setting = data_setting,
        r_hat_pre = r_hat_pre,
        r_hat_post = r_hat_post,
        r_hat = r_hat,
        detect_thd = detect_thd,
        V.diag = V.diag,
        lrv = lrv,
        id_thd_coef = id_thd_coef,
        seed = 900 + sim,
        use_true_cp = use_true_cp,
        N_exact = N_exact
      )
    }
    
    meta_df <- do.call(rbind, lapply(seq_len(nrep), function(sim) {
      x <- out_list[[sim]]
      
      valid_detect <- if (isTRUE(use_true_cp)) {
        !is.null(x$loss_mode)
      } else {
        isTRUE(x$valid)
      }
      
      valid_idt <- if (data_setting == "s0") {
        valid_detect &&
          !isTRUE(x$changed_hat[1]) &&
          !isTRUE(x$changed_hat[2]) &&
          !isTRUE(x$changed_hat[3])
      } else {
        valid_detect &&
          isTRUE(x$changed_hat[1]) &&
          !isTRUE(x$changed_hat[2]) &&
          !isTRUE(x$changed_hat[3])
      }
      
      data.frame(
        sim = sim,
        Time = Time_i,
        dim_obs = dim_chr,
        valid_detect = valid_detect,
        valid_idt = valid_idt,
        theta_true = x$theta_true,
        cp_hat = x$cp_hat,
        cp_used = x$cp_used,
        cp_err = x$cp_err,
        mode1_changed = x$changed_hat[1],
        mode2_changed = x$changed_hat[2],
        mode3_changed = x$changed_hat[3],
        mode1_id_correct = x$mode_id_correct[1],
        mode2_id_correct = x$mode_id_correct[2],
        mode3_id_correct = x$mode_id_correct[3],
        r1_hat = x$r_hat_full[1],
        r2_hat = x$r_hat_full[2],
        r3_hat = x$r_hat_full[3],
        r1_hat_pre = x$r_hat_pre[1],
        r2_hat_pre = x$r_hat_pre[2],
        r3_hat_pre = x$r_hat_pre[3],
        r1_hat_post = x$r_hat_post[1],
        r2_hat_post = x$r_hat_post[2],
        r3_hat_post = x$r_hat_post[3]
      )
    }))
    
    loss_mode_df <- do.call(rbind, lapply(seq_len(nrep), function(sim) {
      x <- out_list[[sim]]
      if (!isTRUE(x$valid) || is.null(x$loss_mode)) return(NULL)
      cbind(
        sim = sim,
        Time = Time_i,
        dim_obs = dim_chr,
        x$loss_mode
      )
    }))
    
    valid_idt_keys <- meta_df %>%
      dplyr::filter(valid_idt) %>%
      dplyr::select(sim, Time, dim_obs)
    
    loss_mode_df_idt <- if (!is.null(loss_mode_df) && nrow(loss_mode_df)) {
      dplyr::semi_join(loss_mode_df, valid_idt_keys, by = c("sim", "Time", "dim_obs"))
    } else {
      loss_mode_df
    }
    
    res_store[[i]] <- list(
      meta = meta_df,
      loss_mode = loss_mode_df,
      loss_mode_idt = loss_mode_df_idt,
      raw = out_list,
      setting = list(Time = Time_i, dim_obs = dim_obs_i)
    )
  }
  
  meta_all <- do.call(rbind, lapply(res_store, `[[`, "meta"))
  loss_mode_all <- do.call(rbind, lapply(res_store, `[[`, "loss_mode"))
  loss_mode_all_idt <- do.call(rbind, lapply(res_store, `[[`, "loss_mode_idt"))
  
  error_table <- do.call(rbind, lapply(res_store, function(obj) {
    meta_df <- obj$meta
    loss_mode_df <- obj$loss_mode_idt
    
    valid_detect_n <- sum(meta_df$valid_detect, na.rm = TRUE)
    valid_idt_n <- sum(meta_df$valid_idt, na.rm = TRUE)
    
    get_mean <- function(mode_k, var_nm) {
      if (!is.null(loss_mode_df) && nrow(loss_mode_df)) {
        mean(loss_mode_df[[var_nm]][loss_mode_df$mode == mode_k], na.rm = TRUE)
      } else {
        NA_real_
      }
    }
    
    data.frame(
      Time = unique(meta_df$Time),
      dim_obs = unique(meta_df$dim_obs),
      valid_detect = valid_detect_n,
      valid_idt = valid_idt_n,
      cp_true = unique(meta_df$theta_true),
      
      mode1_pre_M1  = get_mean(1, "M1_pre"),
      mode1_pre_M2  = get_mean(1, "M2_pre"),
      mode1_pre_M3  = get_mean(1, "M3_pre"),
      mode1_post_M1 = get_mean(1, "M1_post"),
      mode1_post_M2 = get_mean(1, "M2_post"),
      mode1_post_M3 = get_mean(1, "M3_post"),
      
      mode2_pre_M1  = get_mean(2, "M1_pre"),
      mode2_M2      = get_mean(2, "M2_pool"),
      mode2_pre_M3  = get_mean(2, "M3_pre"),
      mode2_post_M1 = get_mean(2, "M1_post"),
      mode2_post_M3 = get_mean(2, "M3_post"),
      
      mode3_pre_M1  = get_mean(3, "M1_pre"),
      mode3_M2      = get_mean(3, "M2_pool"),
      mode3_pre_M3  = get_mean(3, "M3_pre"),
      mode3_post_M1 = get_mean(3, "M1_post"),
      mode3_post_M3 = get_mean(3, "M3_post")
    )
  }))
  
  error_table <- error_table[order(error_table$Time, error_table$dim_obs), ]
  
  id_table <- do.call(rbind, lapply(res_store, function(obj) {
    meta_df <- obj$meta
    meta_detect <- meta_df[meta_df$valid_detect, , drop = FALSE]
    
    valid_detect_n <- sum(meta_df$valid_detect, na.rm = TRUE)
    valid_idt_n <- sum(meta_df$valid_idt, na.rm = TRUE)
    
    if (nrow(meta_detect) == 0L) {
      return(data.frame(
        Time = unique(meta_df$Time),
        dim_obs = unique(meta_df$dim_obs),
        valid_detect = valid_detect_n,
        valid_idt = valid_idt_n,
        mode1_tpr = NA_real_,
        mode2_fpr = NA_real_,
        mode3_fpr = NA_real_
      ))
    }
    
    data.frame(
      Time = unique(meta_df$Time),
      dim_obs = unique(meta_df$dim_obs),
      valid_detect = valid_detect_n,
      valid_idt = valid_idt_n,
      mode1_tpr = mean(meta_detect$mode1_changed, na.rm = TRUE),
      mode2_fpr = mean(meta_detect$mode2_changed, na.rm = TRUE),
      mode3_fpr = mean(meta_detect$mode3_changed, na.rm = TRUE)
    )
  }))
  
  id_table <- id_table[order(id_table$Time, id_table$dim_obs), ]
  
  list(
    error_table = error_table,
    id_table = id_table,
    meta_all = meta_all,
    loss_mode_all = loss_mode_all,
    loss_mode_all_idt = loss_mode_all_idt,
    by_setting = res_store
  )
}

plot_reest_box <- function(reest_obj,
                           Time_select,
                           dim_obs_list,
                           mode = c(1, 2, 3),
                           show_M3 = FALSE,
                           y_lim = NULL,
                           outlier_size = 1.2,
                           box_width = 0.65,
                           title_suffix = NULL,
                           legend_title = NULL,
                           legend_labels = NULL) {
  mode <- as.integer(mode)
  if (!all(mode %in% c(1L, 2L, 3L))) {
    stop("`mode` must be chosen from 1, 2 and 3.")
  }
  
  loss_df <- reest_obj$loss_mode_all_idt
  
  if (is.null(loss_df) || nrow(loss_df) == 0L) {
    stop("No loss data available in `reest_obj$loss_mode_all_idt`.")
  }
  
  if (length(Time_select) != 1L) {
    stop("`Time_select` must be a single value.")
  }
  
  dim_obs_chr_list <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
  
  df_sub <- loss_df %>%
    dplyr::filter(
      Time == !!Time_select,
      dim_obs %in% !!dim_obs_chr_list,
      mode %in% !!mode
    )
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found for the given `Time_select`, `dim_obs_list`, and `mode`.")
  }
  
  df_list <- list()
  
  for (m in mode) {
    df_m <- df_sub %>% dplyr::filter(mode == m)
    
    if (m == 1L) {
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "Pre",
          error = M1_pre,
          x_pos = 1
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "Post",
          error = M1_post,
          x_pos = 2
        )
    }
    
    if (m == 2L) {
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M1 pre",
          error = M1_pre,
          x_pos = 4
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M1 post",
          error = M1_post,
          x_pos = 5
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M2",
          error = M2_pool,
          x_pos = 6
        )
      
      if (isTRUE(show_M3)) {
        df_list[[length(df_list) + 1L]] <- df_m %>%
          dplyr::transmute(
            sim, Time, dim_obs, mode,
            box_type = "M3",
            error = 0.5 * (M3_pre + M3_post),
            x_pos = 7
          )
      }
    }
    
    if (m == 3L) {
      base_pos <- if (isTRUE(show_M3)) 9 else 8
      
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M1 pre",
          error = M1_pre,
          x_pos = base_pos
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M1 post",
          error = M1_post,
          x_pos = base_pos + 1
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M2",
          error = M2_pool,
          x_pos = base_pos + 2
        )
      
      if (isTRUE(show_M3)) {
        df_list[[length(df_list) + 1L]] <- df_m %>%
          dplyr::transmute(
            sim, Time, dim_obs, mode,
            box_type = "M3",
            error = 0.5 * (M3_pre + M3_post),
            x_pos = base_pos + 3
          )
      }
    }
  }
  
  df_plot <- dplyr::bind_rows(df_list) %>%
    dplyr::filter(!is.na(error))
  
  if (nrow(df_plot) == 0L) {
    stop("No non-missing errors available for plotting.")
  }
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      fill_key = dplyr::case_when(
        box_type %in% c("Pre", "M1 pre") ~ "M1 pre",
        box_type %in% c("Post", "M1 post") ~ "M1 post",
        box_type == "M2" ~ "M2",
        box_type == "M3" ~ "M3"
      )
    )
  
  used_fill_levels <- c("M1 pre", "M1 post", "M2")
  if (isTRUE(show_M3)) {
    used_fill_levels <- c(used_fill_levels, "M3")
  }
  
  if (is.null(legend_labels)) {
    legend_labels <- used_fill_levels
  }
  
  if (!is.character(legend_labels) || length(legend_labels) != length(used_fill_levels)) {
    stop("`legend_labels` must be NULL or a character vector with length equal to the number of legend entries.")
  }
  
  names(legend_labels) <- used_fill_levels
  
  dim_obs_labels <- paste0("(", gsub("x", ",", dim_obs_chr_list), ")")
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      fill_key = factor(
        fill_key,
        levels = used_fill_levels,
        labels = unname(legend_labels)
      ),
      Time = factor(
        Time,
        levels = Time_select,
        labels = paste0("T = ", Time_select)
      ),
      dim_obs = factor(
        dim_obs,
        levels = dim_obs_chr_list,
        labels = dim_obs_labels
      )
    )
  
  box_colours <- c(
    "M1 pre"  = "#F4B6C2",
    "M1 post" = "#D97A9A",
    "M2"      = "#A1D99B",
    "M3"      = "#9ECAE1"
  )
  
  # if (!is.null(title_suffix) && nzchar(title_suffix)) {
  #   title_main <- title_suffix
  # } else {
  #   title_main <- paste0("T = ", Time_select)
  # }
  
  vline_pos <- c()
  if (all(c(1L, 2L) %in% mode)) {
    vline_pos <- c(vline_pos, 3)
  }
  if (all(c(2L, 3L) %in% mode)) {
    vline_pos <- c(vline_pos, if (isTRUE(show_M3)) 8 else 7)
  }
  
  x_breaks <- c()
  x_labels <- c()
  if (1L %in% mode) {
    x_breaks <- c(x_breaks, 1.5)
    x_labels <- c(x_labels, "Mode 1")
  }
  if (2L %in% mode) {
    x_breaks <- c(x_breaks, if (isTRUE(show_M3)) 5.5 else 5)
    x_labels <- c(x_labels, "Mode 2")
  }
  if (3L %in% mode) {
    x_breaks <- c(x_breaks, if (isTRUE(show_M3)) 10.5 else 9)
    x_labels <- c(x_labels, "Mode 3")
  }
  
  p <- ggplot(df_plot, aes(x = x_pos, y = error, fill = fill_key, group = x_pos)) +
    {if (length(vline_pos) > 0) geom_vline(xintercept = vline_pos, linetype = "dashed",
                                           colour = "grey50", linewidth = 0.4)} +
    geom_boxplot(
      width = box_width,
      outlier.size = outlier_size,
      outlier.alpha = 0.75
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels
    ) +
    scale_fill_manual(
      values = unname(box_colours[used_fill_levels]),
      breaks = unname(legend_labels),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = "Loading estimation error",
      fill = legend_title
      #title = title_main
    ) +
    facet_grid(Time ~ dim_obs, scales = "fixed") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 12),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y.right = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(y_lim)) {
    if (!is.numeric(y_lim) || length(y_lim) != 2L) {
      stop("`y_lim` must be NULL or a numeric vector of length 2.")
    }
    p <- p + coord_cartesian(ylim = y_lim)
  }
  
  p
}


library(dplyr)
library(tidyr)
library(ggplot2)

plot_reest_box_all <- function(reest_obj,
                               Time_list,
                               dim_obs_list,
                               mode = c(1, 2, 3),
                               show_M3 = FALSE,
                               y_lim = c(0, 0.1),
                               outlier_size = 1.2,
                               box_width = 0.65,
                               title_suffix = NULL,
                               legend_title = NULL,
                               legend_labels = NULL) {
  mode <- as.integer(mode)
  if (!all(mode %in% c(1L, 2L, 3L))) {
    stop("`mode` must be chosen from 1, 2 and 3.")
  }
  
  loss_df <- reest_obj$loss_mode_all_idt
  
  if (is.null(loss_df) || nrow(loss_df) == 0L) {
    stop("No loss data available in `reest_obj$loss_mode_all_idt`.")
  }
  
  if (length(Time_list) == 0L) {
    stop("`Time_list` must contain at least one value.")
  }
  
  dim_obs_chr_list <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
  
  df_sub <- loss_df %>%
    dplyr::filter(
      Time %in% Time_list,
      dim_obs %in% dim_obs_chr_list,
      mode %in% mode
    )
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found for the given `Time_list`, `dim_obs_list`, and `mode`.")
  }
  
  df_list <- list()
  
  for (m in mode) {
    df_m <- df_sub %>% dplyr::filter(mode == m)
    
    if (m == 1L) {
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "Pre",
          error = M1_pre
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "Post",
          error = M1_post
        )
    } else {
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M1 pre",
          error = M1_pre
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M1 post",
          error = M1_post
        )
      df_list[[length(df_list) + 1L]] <- df_m %>%
        dplyr::transmute(
          sim, Time, dim_obs, mode,
          box_type = "M2",
          error = M2_pool
        )
      
      if (isTRUE(show_M3)) {
        df_list[[length(df_list) + 1L]] <- df_m %>%
          dplyr::transmute(
            sim, Time, dim_obs, mode,
            box_type = "M3",
            error = 0.5 * (M3_pre + M3_post)
          )
      }
    }
  }
  
  df_plot <- dplyr::bind_rows(df_list) %>%
    dplyr::filter(!is.na(error))
  
  if (nrow(df_plot) == 0L) {
    stop("No non-missing errors available for plotting.")
  }
  
  all_box_levels <- c("Pre", "Post", "M1 pre", "M1 post", "M2", "M3")
  used_box_levels <- all_box_levels[all_box_levels %in% unique(df_plot$box_type)]
  
  if (is.null(legend_labels)) {
    legend_labels <- used_box_levels
  }
  
  if (!is.character(legend_labels) || length(legend_labels) != length(used_box_levels)) {
    stop("`legend_labels` must be NULL or a character vector with length equal to the number of displayed box types.")
  }
  
  names(legend_labels) <- used_box_levels
  
  dim_obs_labels <- paste0("(", gsub("x", ",", dim_obs_chr_list), ")")
  mode_labels <- c("Mode 1", "Mode 2", "Mode 3")
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      box_type = factor(
        box_type,
        levels = used_box_levels,
        labels = unname(legend_labels)
      ),
      mode = factor(
        mode,
        levels = c(1, 2, 3),
        labels = mode_labels
      ),
      Time = factor(
        Time,
        levels = Time_list,
        labels = paste0("T = ", Time_list)
      ),
      dim_obs = factor(
        dim_obs,
        levels = dim_obs_chr_list,
        labels = dim_obs_labels
      ),
      x_dummy = ""
    )
  
  box_colours <- c(
    "Pre" = "#F4B6C2",
    "Post" = "#D97A9A",
    "M1 pre" = "#F4B6C2",
    "M1 post" = "#D97A9A",
    "M2" = "#A1D99B",
    "M3" = "#9ECAE1"
  )
  
  title_main <- paste(levels(droplevels(df_plot$mode)), collapse = ", ")
  if (!is.null(title_suffix) && nzchar(title_suffix)) {
    title_main <- paste0(title_main, ": ", title_suffix)
  }
  
  facet_scales <- if (is.list(y_lim)) "free_y" else "fixed"
  
  p <- ggplot(df_plot, aes(x = x_dummy, y = error, fill = box_type)) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = box_width,
      outlier.size = outlier_size,
      outlier.alpha = 0.75
    ) +
    scale_fill_manual(
      values = unname(box_colours[used_box_levels]),
      breaks = unname(legend_labels),
      drop = TRUE
    ) +
    labs(
      x = NULL,
      y = "Loading estimation error",
      fill = legend_title,
      title = title_main
    ) +
    facet_grid(dim_obs ~ Time, scales = facet_scales) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold", angle = 270),
      strip.text.y.right = element_text(size = 10, face = "bold", angle = 270),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ## y_lim can be:
  ## 1. NULL -> do nothing
  ## 2. numeric vector length 2 -> same limit for all panels
  ## 3. list of length = number of row panels -> row-specific limits
  if (!is.null(y_lim)) {
    n_row <- length(levels(droplevels(df_plot$mode))) * length(dim_obs_labels)
    
    if (is.numeric(y_lim) && length(y_lim) == 2L) {
      p <- p + coord_cartesian(ylim = y_lim)
    } else if (is.list(y_lim)) {
      if (!requireNamespace("ggh4x", quietly = TRUE)) {
        stop("Package `ggh4x` is required for row-specific y limits.")
      }
      if (length(y_lim) != n_row) {
        stop("When `y_lim` is a list, its length must equal the number of facet rows, i.e. length(selected modes) × length(dim_obs_list).")
      }
      
      y_scales <- lapply(seq_len(n_row), function(i) {
        lim_i <- y_lim[[i]]
        if (is.null(lim_i)) {
          ggplot2::scale_y_continuous()
        } else {
          if (!is.numeric(lim_i) || length(lim_i) != 2L) {
            stop("Each element of `y_lim` must be NULL or a numeric vector of length 2.")
          }
          ggplot2::scale_y_continuous(limits = lim_i)
        }
      })
      
      p <- p + ggh4x::facetted_pos_scales(y = y_scales)
    } else {
      stop("`y_lim` must be NULL, a numeric vector of length 2, or a list.")
    }
  }
  
  p
}


plot_reest_box_all_mode <- function(reest_obj,
                                    Time_list,
                                    dim_obs_list,
                                    y_lim = NULL,
                                    outlier_size = 1.2,
                                    box_width = 0.65,
                                    title_suffix = NULL,
                                    legend_title = NULL,
                                    legend_labels = NULL) {
  loss_df <- reest_obj$loss_mode_all_idt
  
  if (is.null(loss_df) || nrow(loss_df) == 0L) {
    stop("No loss data available in `reest_obj$loss_mode_all_idt`.")
  }
  
  if (length(Time_list) == 0L) {
    stop("`Time_list` must contain at least one value.")
  }
  
  dim_obs_chr_list <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
  
  df_sub <- loss_df %>%
    dplyr::filter(
      Time %in% Time_list,
      dim_obs %in% dim_obs_chr_list,
      mode %in% c(1L, 2L, 3L)
    )
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found for the given `Time_list` and `dim_obs_list`.")
  }
  
  df_list <- list()
  
  ## Mode 1: Pre / Post
  df_m1 <- df_sub %>% dplyr::filter(mode == 1L)
  if (nrow(df_m1) > 0L) {
    df_list[[length(df_list) + 1L]] <- df_m1 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "Pre",
        error = M1_pre,
        x_pos = 1
      )
    df_list[[length(df_list) + 1L]] <- df_m1 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "Post",
        error = M1_post,
        x_pos = 2
      )
  }
  
  ## Mode 2: M1 pre / M1 post / M2
  df_m2 <- df_sub %>% dplyr::filter(mode == 2L)
  if (nrow(df_m2) > 0L) {
    df_list[[length(df_list) + 1L]] <- df_m2 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "M1 pre",
        error = M1_pre,
        x_pos = 4
      )
    df_list[[length(df_list) + 1L]] <- df_m2 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "M1 post",
        error = M1_post,
        x_pos = 5
      )
    df_list[[length(df_list) + 1L]] <- df_m2 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "M2",
        error = M2_pool,
        x_pos = 6
      )
  }
  
  ## Mode 3: M1 pre / M1 post / M2
  df_m3 <- df_sub %>% dplyr::filter(mode == 3L)
  if (nrow(df_m3) > 0L) {
    df_list[[length(df_list) + 1L]] <- df_m3 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "M1 pre",
        error = M1_pre,
        x_pos = 8
      )
    df_list[[length(df_list) + 1L]] <- df_m3 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "M1 post",
        error = M1_post,
        x_pos = 9
      )
    df_list[[length(df_list) + 1L]] <- df_m3 %>%
      dplyr::transmute(
        sim, Time, dim_obs, mode,
        method = "M2",
        error = M2_pool,
        x_pos = 10
      )
  }
  
  df_plot <- dplyr::bind_rows(df_list) %>%
    dplyr::filter(!is.na(error))
  
  if (nrow(df_plot) == 0L) {
    stop("No non-missing errors available for plotting.")
  }
  
  ## Map plotting categories to a smaller legend
  df_plot <- df_plot %>%
    dplyr::mutate(
      fill_key = dplyr::case_when(
        method %in% c("Pre", "M1 pre") ~ "M1 pre",
        method %in% c("Post", "M1 post") ~ "M1 post",
        method == "M2" ~ "M2"
      )
    )
  
  used_fill_levels <- c("M1 pre", "M1 post", "M2")
  
  if (is.null(legend_labels)) {
    legend_labels <- used_fill_levels
  }
  
  if (!is.character(legend_labels) || length(legend_labels) != length(used_fill_levels)) {
    stop("`legend_labels` must be NULL or a character vector of length 3.")
  }
  
  names(legend_labels) <- used_fill_levels
  
  dim_obs_labels <- paste0("(", gsub("x", ",", dim_obs_chr_list), ")")
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      fill_key = factor(
        fill_key,
        levels = used_fill_levels,
        labels = unname(legend_labels)
      ),
      Time = factor(
        Time,
        levels = Time_list,
        labels = paste0("T = ", Time_list)
      ),
      dim_obs = factor(
        dim_obs,
        levels = dim_obs_chr_list,
        labels = dim_obs_labels
      )
    )
  
  box_colours <- c(
    "M1 pre"  = "#F4B6C2",
    "M1 post" = "#D97A9A",
    "M2"      = "#A1D99B"
  )
  
  #title_main <- "Modes 1, 2 and 3"
  # if (!is.null(title_suffix) && nzchar(title_suffix)) {
  #   title_main <- paste0(title_suffix)
  # }
  
  facet_scales <- if (is.list(y_lim)) "free_y" else "fixed"
  
  p <- ggplot(df_plot, aes(x = x_pos, y = error, fill = fill_key, group = x_pos)) +
    geom_vline(xintercept = c(3, 7), linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_boxplot(
      width = box_width,
      outlier.size = outlier_size,
      outlier.alpha = 0.75
    ) +
    scale_x_continuous(
      breaks = c(1.5, 5, 9),
      labels = c("Mode 1", "Mode 2", "Mode 3")
    ) +
    scale_fill_manual(
      values = unname(box_colours[used_fill_levels]),
      breaks = unname(legend_labels),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = "Loading estimation error",
      fill = legend_title
      #title = title_main
    ) +
    facet_grid(dim_obs ~ Time, scales = facet_scales) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold", angle = 270),
      strip.text.y.right = element_text(size = 10, face = "bold", angle = 270),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(y_lim)) {
    if (is.numeric(y_lim) && length(y_lim) == 2L) {
      p <- p + coord_cartesian(ylim = y_lim)
    } else if (is.list(y_lim)) {
      if (!requireNamespace("ggh4x", quietly = TRUE)) {
        stop("Package `ggh4x` is required for row-specific y limits.")
      }
      
      n_row <- length(dim_obs_labels)
      if (length(y_lim) != n_row) {
        stop("When `y_lim` is a list, its length must equal length(dim_obs_list).")
      }
      
      y_scales <- lapply(seq_len(n_row), function(i) {
        lim_i <- y_lim[[i]]
        if (is.null(lim_i)) {
          ggplot2::scale_y_continuous()
        } else {
          if (!is.numeric(lim_i) || length(lim_i) != 2L) {
            stop("Each element of `y_lim` must be NULL or a numeric vector of length 2.")
          }
          ggplot2::scale_y_continuous(limits = lim_i)
        }
      })
      
      p <- p + ggh4x::facetted_pos_scales(y = y_scales)
    } else {
      stop("`y_lim` must be NULL, a numeric vector of length 2, or a list.")
    }
  }
  
  p
}

check_outliers <- function(reest_obj,
                           Time,
                           dim_obs,
                           mode = 2,
                           methods = c("M1", "M2", "M3"),
                           error_lim = c(0.1, Inf),
                           include_M3 = TRUE,
                           sort_by_max_error = TRUE) {
  methods <- unique(methods)
  methods <- match.arg(methods, choices = c("M1", "M2", "M3"), several.ok = TRUE)
  
  mode <- as.integer(mode)
  if (!mode %in% c(2L, 3L)) {
    stop("`mode` must be 2 or 3.")
  }
  
  if (!is.numeric(error_lim) || length(error_lim) != 2L) {
    stop("`error_lim` must be a numeric vector of length 2.")
  }
  if (error_lim[1] > error_lim[2]) {
    stop("`error_lim` must satisfy error_lim[1] <= error_lim[2].")
  }
  
  loss_df <- reest_obj$loss_mode_all_idt
  meta_df <- reest_obj$meta_all
  dim_obs_chr <- paste(dim_obs, collapse = "x")
  
  df_sub <- loss_df %>%
    dplyr::filter(Time == !!Time, dim_obs == !!dim_obs_chr, mode == !!mode)
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found.")
  }
  
  df_list <- list()
  
  if ("M1" %in% methods) {
    df_list[["M1_pre"]] <- df_sub %>%
      dplyr::transmute(sim, box_type = "M1_pre", error = M1_pre)
    
    df_list[["M1_post"]] <- df_sub %>%
      dplyr::transmute(sim, box_type = "M1_post", error = M1_post)
  }
  
  if ("M2" %in% methods) {
    df_list[["M2"]] <- df_sub %>%
      dplyr::transmute(sim, box_type = "M2", error = M2_pool)
  }
  
  if ("M3" %in% methods && isTRUE(include_M3)) {
    df_list[["M3"]] <- df_sub %>%
      dplyr::transmute(sim, box_type = "M3", error = 0.5 * (M3_pre + M3_post))
  }
  
  df_check <- dplyr::bind_rows(df_list)
  
  sim_keep <- df_check %>%
    dplyr::filter(
      !is.na(error),
      error >= error_lim[1],
      error <= error_lim[2]
    ) %>%
    dplyr::group_by(sim) %>%
    dplyr::summarise(
      max_error = max(error),
      which_box = paste(box_type[which.max(error)], collapse = ", "),
      .groups = "drop"
    )
  
  if (nrow(sim_keep) == 0L) {
    message(
      "No outliers found for Time = ", Time,
      ", dim_obs = ", dim_obs_chr,
      ", mode = ", mode,
      " within error_lim = [", error_lim[1], ", ", error_lim[2], "]."
    )
    return(data.frame())
  }
  
  out <- meta_df %>%
    dplyr::filter(Time == !!Time, dim_obs == !!dim_obs_chr) %>%
    dplyr::semi_join(sim_keep, by = "sim") %>%
    dplyr::select(
      sim,
      valid_detect, valid_idt,
      cp_hat, cp_err,
      mode1_changed, mode2_changed, mode3_changed,
      r1_hat, r2_hat, r3_hat,
      r1_hat_pre, r2_hat_pre, r3_hat_pre,
      r1_hat_post, r2_hat_post, r3_hat_post
    )
  
  out <- sim_keep %>%
    {if (sort_by_max_error) dplyr::arrange(., dplyr::desc(max_error)) else dplyr::arrange(., sim)} %>%
    dplyr::left_join(out, by = "sim")
  
  out
}