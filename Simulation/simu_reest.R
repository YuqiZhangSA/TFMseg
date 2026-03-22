mode_informed_proj <- function(X, dim_X, st, ed, k,
                               Lambda1_seg, Lambda_comp_pool,
                               r_k, centre = TRUE) {
  stopifnot(length(dim_X) == 3L)
  stopifnot(k %in% c(2L, 3L))
  
  if (centre) {
    X_use <- centre_along_axis(X, axis = length(dim(X)))
  } else {
    X_use <- X
  }
  
  idx_time <- st:ed
  T_win <- length(idx_time)
  X_win <- slice_time(X_use, idx_time)
  p_k <- dim_X[k]
  
  if (k == 2L) {
    p_mk <- dim_X[1] * dim_X[3]
    Lambda_kron <- kronecker_list(list(Lambda_comp_pool, Lambda1_seg))
  } else {
    p_mk <- dim_X[1] * dim_X[2]
    Lambda_kron <- kronecker_list(list(Lambda_comp_pool, Lambda1_seg))
  }
  
  perm <- c(k, setdiff(seq_len(3L), k), 4L)
  Xp <- aperm(X_win, perm)
  
  sum_Y <- matrix(0, nrow = p_k, ncol = p_k)
  for (t in seq_len(T_win)) {
    Xkt <- matrix(slice_time(Xp, t), nrow = p_k)
    Ykt_t <- (1 / p_mk) * Xkt %*% Lambda_kron
    sum_Y <- sum_Y + Ykt_t %*% t(Ykt_t)
  }
  
  Gamma_Y <- sum_Y / (T_win * p_k)
  eigY <- eigen(Gamma_Y, symmetric = TRUE)
  eigY$vectors[, seq_len(r_k), drop = FALSE] * sqrt(p_k)
}


compare_reest <- function(
    dim_obs,
    dim_latent,
    Time,
    dep = TRUE,
    theta_coef = 0.5,
    r_hat_pre = NULL,
    r_hat_post = NULL,
    r_hat = NULL,
    detect_thd = NULL,
    V.diag = TRUE,
    lrv = TRUE,
    id_thd_coef = 3,
    seed = 901
) {
  stopifnot(length(dim_obs) == 3L, length(dim_latent) == 3L)
  stopifnot(all(dim_latent >= 3))
  
  theta_true <- floor(Time * theta_coef)
  if (theta_true <= 0 || theta_true >= Time) {
    stop("theta_true must lie in {1, ..., Time - 1}.")
  }
  
  set.seed(seed)
  
  A2 <- matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 0
  ), nrow = 3, byrow = TRUE)
  
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
      list(A2, diag(dim_latent[2]), diag(dim_latent[3]))
    ),
    add_factors = NULL,
    add_factors_coeff = NULL,
    dist = "Gaussian"
  )
  
  X <- data_sim$X
  dim_X <- dim(X)[1:3]
  
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
    209.8954613 * sqrt(log(Time)) +
      0.7127491 * sqrt(dr) +
      1566.8875353 * sqrt(1 / log(Time)) +
      1572.7173337 * log(log(Time)) / sqrt(log(Time)) -
      2298.3882769
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
  valid <- (length(cp_hat_all) == 1L) && (abs(cp_hat_all[1] - theta_true) <= acc_win)
  
  if (!valid) {
    return(list(
      valid = FALSE,
      theta_true = theta_true,
      cp_hat = NA_integer_,
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
  
  cp_hat <- cp_hat_all[1]
  
  id_out <- id_modes(
    G = G,
    G_dim = G_dim,
    dim_obs = dim_obs,
    cp_vec = cp_hat,
    st = 1L,
    ed = Time,
    id_thd_coef = id_thd_coef
  )
  
  changed_hat <- as.logical(id_out$changed[, 1])
  true_changed <- c(TRUE, FALSE, FALSE)
  mode_id_correct <- (changed_hat == true_changed)
  
  true_pre <- data_sim$loadings_hist[[1]]
  true_post <- data_sim$loadings_hist[[theta_true + 1L]]
  
  fit_pre <- global_pca(
    X = X, dim_X = dim_X, r_hat = r_hat_pre,
    st = 1L, ed = cp_hat,
    centre = TRUE, proj = TRUE
  )
  
  fit_post <- global_pca(
    X = X, dim_X = dim_X, r_hat = r_hat_post,
    st = cp_hat + 1L, ed = Time,
    centre = TRUE, proj = TRUE
  )
  
  r_full_used <- as.integer(fit_full0$r_hat)
  r_pre_used  <- as.integer(fit_pre$r_hat)
  r_post_used <- as.integer(fit_post$r_hat)
  
  Lambda_init_pool <- fit_full0$Lambda_init
  
  ## M2
  Lambda_M2_pre <- vector("list", 3L)
  Lambda_M2_post <- vector("list", 3L)
  Lambda_M2_pre[[1]] <- fit_pre$Lambda_proj[[1]]
  Lambda_M2_post[[1]] <- fit_post$Lambda_proj[[1]]
  
  Lambda_M2_pre[[2]] <- mode_informed_proj(
    X = X, dim_X = dim_X,
    st = 1L, ed = cp_hat, k = 2L,
    Lambda1_seg = fit_pre$Lambda_proj[[1]],
    Lambda_comp_pool = Lambda_init_pool[[3]],
    r_k = r_pre_used[2],
    centre = TRUE
  )
  
  Lambda_M2_post[[2]] <- mode_informed_proj(
    X = X, dim_X = dim_X,
    st = cp_hat + 1L, ed = Time, k = 2L,
    Lambda1_seg = fit_post$Lambda_proj[[1]],
    Lambda_comp_pool = Lambda_init_pool[[3]],
    r_k = r_post_used[2],
    centre = TRUE
  )
  
  Lambda_M2_pre[[3]] <- mode_informed_proj(
    X = X, dim_X = dim_X,
    st = 1L, ed = cp_hat, k = 3L,
    Lambda1_seg = fit_pre$Lambda_proj[[1]],
    Lambda_comp_pool = Lambda_init_pool[[2]],
    r_k = r_pre_used[3],
    centre = TRUE
  )
  
  Lambda_M2_post[[3]] <- mode_informed_proj(
    X = X, dim_X = dim_X,
    st = cp_hat + 1L, ed = Time, k = 3L,
    Lambda1_seg = fit_post$Lambda_proj[[1]],
    Lambda_comp_pool = Lambda_init_pool[[2]],
    r_k = r_post_used[3],
    centre = TRUE
  )
  
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
    true_pre_M2  <- true_match(true_pre[[k]],  Lambda_M2_pre[[k]])
    true_post_M2 <- true_match(true_post[[k]], Lambda_M2_post[[k]])
    
    M2_pre  <- tensorMiss::fle(Lambda_M2_pre[[k]],  true_pre_M2)
    M2_post <- tensorMiss::fle(Lambda_M2_post[[k]], true_post_M2)
    
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
      M3_pre = M3_pre,
      M3_post = M3_post
    )
  }))
  
  list(
    valid = TRUE,
    theta_true = theta_true,
    cp_hat = cp_hat,
    cp_hat_all = cp_hat_all,
    cp_err = cp_hat - theta_true,
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
    id_thd_coef = 3,
    data_setting = "s3_A2"
) {
  if (is.null(Time_list)) stop("Please provide Time_list.")
  if (is.null(dim_obs_list)) stop("Please provide dim_obs_list.")
  
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
      "[reest_s3A2] %d/%d: Time=%d, dim_obs=%s, setting=%s",
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
        r_hat_pre = r_hat_pre,
        r_hat_post = r_hat_post,
        r_hat = r_hat,
        detect_thd = detect_thd,
        V.diag = V.diag,
        lrv = lrv,
        id_thd_coef = id_thd_coef,
        seed = 900 + sim
      )
    }
    
    meta_df <- do.call(rbind, lapply(seq_len(nrep), function(sim) {
      x <- out_list[[sim]]
      
      valid_detect <- isTRUE(x$valid)
      valid_idt <- valid_detect &&
        isTRUE(x$changed_hat[1]) &&
        !isTRUE(x$changed_hat[2]) &&
        !isTRUE(x$changed_hat[3])
      
      data.frame(
        sim = sim,
        Time = Time_i,
        dim_obs = dim_chr,
        valid_detect = valid_detect,
        valid_idt = valid_idt,
        theta_true = x$theta_true,
        cp_hat = x$cp_hat,
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
      mode2_pre_M2  = get_mean(2, "M2_pre"),
      mode2_pre_M3  = get_mean(2, "M3_pre"),
      mode2_post_M1 = get_mean(2, "M1_post"),
      mode2_post_M2 = get_mean(2, "M2_post"),
      mode2_post_M3 = get_mean(2, "M3_post"),
      
      mode3_pre_M1  = get_mean(3, "M1_pre"),
      mode3_pre_M2  = get_mean(3, "M2_pre"),
      mode3_pre_M3  = get_mean(3, "M3_pre"),
      mode3_post_M1 = get_mean(3, "M1_post"),
      mode3_post_M2 = get_mean(3, "M2_post"),
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

library(dplyr)
library(tidyr)
library(ggplot2)

plot_reest_box <- function(reest_obj,
                           Time,
                           dim_obs,
                           mode = c(2, 3),
                           methods = c("M1", "M2", "M3"),
                           use_valid_idt_only = TRUE,
                           pool_type = c("stack", "average"),
                           y_lim = NULL,
                           outlier_size = 1.8,
                           box_width = 0.65) {
  pool_type <- match.arg(pool_type)
  
  methods <- unique(methods)
  methods <- match.arg(methods, choices = c("M1", "M2", "M3"), several.ok = TRUE)
  
  mode <- as.integer(mode)
  if (!all(mode %in% c(1L, 2L, 3L))) {
    stop("`mode` must be one or more of 1, 2, 3.")
  }
  
  loss_df <- if (use_valid_idt_only) {
    reest_obj$loss_mode_all_idt
  } else {
    reest_obj$loss_mode_all
  }
  
  if (is.null(loss_df) || nrow(loss_df) == 0L) {
    stop("No loss data available in `reest_obj`.")
  }
  
  dim_obs_chr <- paste(dim_obs, collapse = "x")
  
  df_sub <- loss_df %>%
    dplyr::filter(Time == !!Time, dim_obs == !!dim_obs_chr, mode %in% !!mode)
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found for the given `Time`, `dim_obs`, and `mode`.")
  }
  
  main_cols <- as.vector(outer(methods, c("pre", "post"), paste, sep = "_"))
  
  df_long_main <- df_sub %>%
    dplyr::select(sim, mode, dplyr::all_of(main_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(main_cols),
      names_to = c("method", "period"),
      names_sep = "_",
      values_to = "error"
    )
  
  if (pool_type == "stack") {
    df_pool <- df_long_main %>%
      dplyr::mutate(period = "pool")
  } else {
    pool_cols <- as.vector(outer(methods, c("pre", "post"), paste, sep = "_"))
    
    df_pool <- df_sub %>%
      dplyr::select(sim, mode, dplyr::all_of(pool_cols))
    
    for (m in methods) {
      df_pool[[m]] <- 0.5 * (df_pool[[paste0(m, "_pre")]] + df_pool[[paste0(m, "_post")]])
    }
    
    df_pool <- df_pool %>%
      dplyr::select(sim, mode, dplyr::all_of(methods)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(methods),
        names_to = "method",
        values_to = "error"
      ) %>%
      dplyr::mutate(period = "pool")
  }
  
  df_plot <- dplyr::bind_rows(df_long_main, df_pool) %>%
    dplyr::mutate(
      period = factor(period, levels = c("pre", "post", "pool"),
                      labels = c("Pre", "Post", "Pool")),
      method = factor(method, levels = c("M1", "M2", "M3")),
      mode = factor(mode, levels = c(1, 2, 3),
                    labels = c("Mode 1", "Mode 2", "Mode 3"))
    )
  
  method_colours <- c(
    "M1" = "#F4B6C2",  # soft pink
    "M2" = "#9ECAE1",  # soft blue
    "M3" = "#A1D99B"   # soft green
  )
  
  p <- ggplot(df_plot, aes(x = period, y = error, fill = method)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = box_width,
      outlier.size = outlier_size,
      outlier.alpha = 0.75
    ) +
    scale_fill_manual(values = method_colours[methods], drop = FALSE) +
    labs(
      x = NULL,
      y = "Loading estimation error",
      title = paste0(
        paste(levels(droplevels(df_plot$mode)), collapse = ", "),
        ": Time = ", Time,
        ", (p1,p2,p3) = (", dim_obs_chr, ")"
      ),
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  if (length(unique(df_plot$mode)) > 1L) {
    p <- p + facet_wrap(~ mode, nrow = 1)
  }
  
  if (!is.null(y_lim)) {
    p <- p + coord_cartesian(ylim = y_lim)
  }
  

  p
}


plot_reest_box_all <- function(reest_obj,
                               Time_list,
                               dim_obs_list,
                               mode = c(2, 3),
                               methods = c("M1", "M2", "M3"),
                               use_valid_idt_only = TRUE,
                               pool_type = c("stack", "average"),
                               y_lim = FALSE,
                               outlier_size = 1.2,
                               box_width = 0.65,
                               title_suffix = NULL) {
  pool_type <- match.arg(pool_type)
  
  methods <- unique(methods)
  methods <- match.arg(methods, choices = c("M1", "M2", "M3"), several.ok = TRUE)
  
  mode <- as.integer(mode)
  if (!all(mode %in% c(1L, 2L, 3L))) {
    stop("`mode` must be one or more of 1, 2, 3.")
  }
  
  loss_df <- if (use_valid_idt_only) {
    reest_obj$loss_mode_all_idt
  } else {
    reest_obj$loss_mode_all
  }
  
  if (is.null(loss_df) || nrow(loss_df) == 0L) {
    stop("No loss data available in `reest_obj`.")
  }
  
  if (!isTRUE(y_lim)) {
    y_lim <- "free_y"
  } else {
    y_lim <- "fixed"
  }
  
  dim_obs_chr_list <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
  
  df_sub <- loss_df %>%
    dplyr::filter(Time %in% !!Time_list,
                  dim_obs %in% !!dim_obs_chr_list,
                  mode %in% !!mode)
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found for the given `Time_list`, `dim_obs_list`, and `mode`.")
  }
  
  main_cols <- as.vector(outer(methods, c("pre", "post"), paste, sep = "_"))
  
  df_long_main <- df_sub %>%
    dplyr::select(sim, Time, dim_obs, mode, dplyr::all_of(main_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(main_cols),
      names_to = c("method", "period"),
      names_sep = "_",
      values_to = "error"
    )
  
  if (pool_type == "stack") {
    df_pool <- df_long_main %>%
      dplyr::mutate(period = "pool")
  } else {
    pool_cols <- as.vector(outer(methods, c("pre", "post"), paste, sep = "_"))
    
    df_pool <- df_sub %>%
      dplyr::select(sim, Time, dim_obs, mode, dplyr::all_of(pool_cols))
    
    for (m in methods) {
      df_pool[[m]] <- 0.5 * (df_pool[[paste0(m, "_pre")]] + df_pool[[paste0(m, "_post")]])
    }
    
    df_pool <- df_pool %>%
      dplyr::select(sim, Time, dim_obs, mode, dplyr::all_of(methods)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(methods),
        names_to = "method",
        values_to = "error"
      ) %>%
      dplyr::mutate(period = "pool")
  }
  
  df_plot <- dplyr::bind_rows(df_long_main, df_pool) %>%
    dplyr::mutate(
      period = factor(period, levels = c("pre", "post", "pool"),
                      labels = c("Pre", "Post", "Pool")),
      method = factor(method, levels = c("M1", "M2", "M3")),
      mode = factor(mode, levels = c(1, 2, 3),
                    labels = c("Mode 1", "Mode 2", "Mode 3")),
      Time = factor(Time, levels = Time_list,
                    labels = paste0("T = ", Time_list)),
      dim_obs = factor(dim_obs, levels = dim_obs_chr_list,
                       labels = paste0("(", gsub("x", ",", dim_obs_chr_list), ")"))
    )
  
  method_colours <- c(
    "M1" = "#F4B6C2",
    "M2" = "#9ECAE1",
    "M3" = "#A1D99B"
  )
  
  title_main <- paste(levels(droplevels(df_plot$mode)), collapse = ", ")
  if (!is.null(title_suffix) && nzchar(title_suffix)) {
    title_main <- paste0(title_main, ": ", title_suffix)
  }
  
  p <- ggplot(df_plot, aes(x = period, y = error, fill = method)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = box_width,
      outlier.size = outlier_size,
      outlier.alpha = 0.75
    ) +
    scale_fill_manual(values = method_colours[methods], drop = FALSE) +
    labs(
      x = NULL,
      y = "Loading estimation error",
      fill = NULL,
      title = title_main
    ) +
    facet_grid(dim_obs ~ Time, scales = y_lim) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold", angle = 270),
      strip.text.y.right = element_text(size = 10, face = "bold", angle = 270),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  p
}



check_outliers <- function(reest_obj,
                           Time,
                           dim_obs,
                           mode = 2,
                           methods = c("M1", "M2", "M3"),
                           error_lim = c(0.1, Inf),
                           use_valid_idt_only = TRUE,
                           include_pool = TRUE,
                           sort_by_max_error = TRUE) {
  methods <- match.arg(methods, choices = c("M1", "M2", "M3"), several.ok = TRUE)
  
  if (!is.numeric(error_lim) || length(error_lim) != 2L) {
    stop("error_lim must be a numeric vector of length 2.")
  }
  if (error_lim[1] > error_lim[2]) {
    stop("error_lim must satisfy error_lim[1] <= error_lim[2].")
  }
  
  loss_df <- if (use_valid_idt_only) {
    reest_obj$loss_mode_all_idt
  } else {
    reest_obj$loss_mode_all
  }
  
  meta_df <- reest_obj$meta_all
  dim_obs_chr <- paste(dim_obs, collapse = "x")
  
  df_sub <- loss_df %>%
    filter(Time == !!Time, dim_obs == !!dim_obs_chr, mode == !!mode)
  
  if (nrow(df_sub) == 0L) {
    stop("No matching rows found.")
  }
  
  cols_use <- as.vector(outer(methods, c("pre", "post"), paste, sep = "_"))
  
  df_check <- df_sub %>%
    select(sim, all_of(cols_use))
  
  if (include_pool) {
    for (m in methods) {
      df_check[[paste0(m, "_pool")]] <- 0.5 * (df_sub[[paste0(m, "_pre")]] + df_sub[[paste0(m, "_post")]])
    }
  }
  
  sim_keep <- df_check %>%
  pivot_longer(
    cols = -sim,
    names_to = c("method", "period"),
    names_sep = "_",
    values_to = "error"
  ) %>%
  filter(!is.na(error),
         error >= error_lim[1],
         error <= error_lim[2]) %>%
  group_by(sim) %>%
  summarise(max_error = max(error), .groups = "drop")

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
    filter(Time == !!Time, dim_obs == !!dim_obs_chr) %>%
    semi_join(sim_keep, by = "sim") %>%
    select(
      sim,
      r1_hat, r2_hat, r3_hat,
      r1_hat_pre, r2_hat_pre, r3_hat_pre,
      r1_hat_post, r2_hat_post, r3_hat_post
    )
  
  if (sort_by_max_error) {
    out <- sim_keep %>%
      arrange(desc(max_error)) %>%
      select(-max_error) %>%
      left_join(out, by = "sim")
  } else {
    out <- out %>% arrange(sim)
  }
  
  out
}