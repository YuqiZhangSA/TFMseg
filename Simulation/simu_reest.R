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
    trim_coef = 1/4,
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
      loss_mode = NULL,
      loss_overall = NULL
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
  
  true_pre  <- data_sim$loadings_hist[[1]]
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
  
  fit_full <- global_pca(
    X = X, dim_X = dim_X, r_hat = r_hat,
    st = 1L, ed = Time,
    centre = TRUE, proj = TRUE
  )
  
  true_match <- function(Ltrue, Lhat) {
    Ltrue[, seq_len(ncol(Lhat)), drop = FALSE]
  }
  
  loss_mode <- do.call(rbind, lapply(1:3, function(k) {
    true_pre_k  <- true_match(true_pre[[k]],  fit_pre$Lambda_proj[[k]])
    true_post_k <- true_match(true_post[[k]], fit_post$Lambda_proj[[k]])
    
    e1_pre  <- tensorMiss::fle(fit_pre$Lambda_proj[[k]],  true_pre_k)
    e1_post <- tensorMiss::fle(fit_post$Lambda_proj[[k]], true_post_k)
    
    if (changed_hat[k]) {
      e2_pre  <- tensorMiss::fle(fit_pre$Lambda_proj[[k]],  true_pre_k)
      e2_post <- tensorMiss::fle(fit_post$Lambda_proj[[k]], true_post_k)
    } else {
      true_full_pre_k  <- true_match(true_pre[[k]],  fit_full$Lambda_proj[[k]])
      true_full_post_k <- true_match(true_post[[k]], fit_full$Lambda_proj[[k]])
      
      e2_pre  <- tensorMiss::fle(fit_full$Lambda_proj[[k]], true_full_pre_k)
      e2_post <- tensorMiss::fle(fit_full$Lambda_proj[[k]], true_full_post_k)
    }
    
    data.frame(
      mode = k,
      true_changed = true_changed[k],
      changed_hat = changed_hat[k],
      id_correct = mode_id_correct[k],
      e1_pre = e1_pre,
      e1_post = e1_post,
      e2_pre = e2_pre,
      e2_post = e2_post
    )
  }))
  
  loss_overall <- data.frame(
    e1_pre = mean(loss_mode$e1_pre),
    e1_post = mean(loss_mode$e1_post),
    e2_pre = mean(loss_mode$e2_pre),
    e2_post = mean(loss_mode$e2_post)
  )
  
  list(
    valid = TRUE,
    theta_true = theta_true,
    cp_hat = cp_hat,
    cp_hat_all = cp_hat_all,
    cp_err = cp_hat - theta_true,
    changed_hat = changed_hat,
    mode_id_correct = mode_id_correct,
    ratio_Tp = id_out$ratio_Tp[, 1],
    loss_mode = loss_mode,
    loss_overall = loss_overall
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
    trim_coef = 1/4,
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
        trim_coef = trim_coef,
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
        mode3_id_correct = x$mode_id_correct[3]
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
    
    loss_overall_df <- do.call(rbind, lapply(seq_len(nrep), function(sim) {
      x <- out_list[[sim]]
      if (!isTRUE(x$valid) || is.null(x$loss_overall)) return(NULL)
      cbind(
        sim = sim,
        Time = Time_i,
        dim_obs = dim_chr,
        x$loss_overall
      )
    }))
    
    ## keep only valid_idt reps for error summaries
    valid_idt_keys <- meta_df %>%
      dplyr::filter(valid_idt) %>%
      dplyr::select(sim, Time, dim_obs)
    
    loss_mode_df_idt <- if (!is.null(loss_mode_df) && nrow(loss_mode_df)) {
      dplyr::semi_join(loss_mode_df, valid_idt_keys, by = c("sim", "Time", "dim_obs"))
    } else {
      loss_mode_df
    }
    
    loss_overall_df_idt <- if (!is.null(loss_overall_df) && nrow(loss_overall_df)) {
      dplyr::semi_join(loss_overall_df, valid_idt_keys, by = c("sim", "Time", "dim_obs"))
    } else {
      loss_overall_df
    }
    
    res_store[[i]] <- list(
      meta = meta_df,
      loss_mode = loss_mode_df,
      loss_overall = loss_overall_df,
      loss_mode_idt = loss_mode_df_idt,
      loss_overall_idt = loss_overall_df_idt,
      raw = out_list,
      setting = list(Time = Time_i, dim_obs = dim_obs_i)
    )
  }
  
  meta_all <- do.call(rbind, lapply(res_store, `[[`, "meta"))
  loss_mode_all <- do.call(rbind, lapply(res_store, `[[`, "loss_mode"))
  loss_overall_all <- do.call(rbind, lapply(res_store, `[[`, "loss_overall"))
  loss_mode_all_idt <- do.call(rbind, lapply(res_store, `[[`, "loss_mode_idt"))
  loss_overall_all_idt <- do.call(rbind, lapply(res_store, `[[`, "loss_overall_idt"))
  
  ## ---------------------------
  ## Table 1: loading errors
  ## computed on valid_idt only
  ## ---------------------------
  error_table <- do.call(rbind, lapply(res_store, function(obj) {
    meta_df <- obj$meta
    loss_mode_df <- obj$loss_mode_idt
    loss_overall_df <- obj$loss_overall_idt
    
    valid_detect_n <- sum(meta_df$valid_detect, na.rm = TRUE)
    valid_idt_n <- sum(meta_df$valid_idt, na.rm = TRUE)
    
    get_mean <- function(mode_k, var_nm) {
      if (!is.null(loss_mode_df) && nrow(loss_mode_df)) {
        mean(loss_mode_df[[var_nm]][loss_mode_df$mode == mode_k], na.rm = TRUE)
      } else {
        NA_real_
      }
    }
    
    overall_pre_e1 <- if (!is.null(loss_overall_df) && nrow(loss_overall_df)) {
      mean(loss_overall_df$e1_pre, na.rm = TRUE)
    } else NA_real_
    
    overall_post_e1 <- if (!is.null(loss_overall_df) && nrow(loss_overall_df)) {
      mean(loss_overall_df$e1_post, na.rm = TRUE)
    } else NA_real_
    
    overall_pre_e2 <- if (!is.null(loss_overall_df) && nrow(loss_overall_df)) {
      mean(loss_overall_df$e2_pre, na.rm = TRUE)
    } else NA_real_
    
    overall_post_e2 <- if (!is.null(loss_overall_df) && nrow(loss_overall_df)) {
      mean(loss_overall_df$e2_post, na.rm = TRUE)
    } else NA_real_
    
    data.frame(
      Time = unique(meta_df$Time),
      dim_obs = unique(meta_df$dim_obs),
      valid_detect = valid_detect_n,
      valid_idt = valid_idt_n,
      cp_true = unique(meta_df$theta_true),
      
      mode1_pre_e1 = get_mean(1, "e1_pre"),
      mode1_pre_e2 = get_mean(1, "e2_pre"),
      mode1_post_e1 = get_mean(1, "e1_post"),
      mode1_post_e2 = get_mean(1, "e2_post"),
      
      mode2_pre_e1 = get_mean(2, "e1_pre"),
      mode2_pre_e2 = get_mean(2, "e2_pre"),
      mode2_post_e1 = get_mean(2, "e1_post"),
      mode2_post_e2 = get_mean(2, "e2_post"),
      
      mode3_pre_e1 = get_mean(3, "e1_pre"),
      mode3_pre_e2 = get_mean(3, "e2_pre"),
      mode3_post_e1 = get_mean(3, "e1_post"),
      mode3_post_e2 = get_mean(3, "e2_post"),
      
      overall_pre_e1 = overall_pre_e1,
      overall_pre_e2 = overall_pre_e2,
      overall_post_e1 = overall_post_e1,
      overall_post_e2 = overall_post_e2
    )
  }))
  
  error_table <- error_table[order(error_table$Time, error_table$dim_obs), ]
  
  ## ---------------------------
  ## Table 2: identification performance
  ## computed on valid_detect only
  ## ---------------------------
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
    loss_overall_all = loss_overall_all,
    loss_mode_all_idt = loss_mode_all_idt,
    loss_overall_all_idt = loss_overall_all_idt,
    by_setting = res_store
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)

plot_reest <- function(error_table, period = c("pre", "post")) {
  period <- match.arg(period)
  
  cols_use <- if (period == "pre") {
    c("mode1_pre_e1", "mode1_pre_e2",
      "mode2_pre_e1", "mode2_pre_e2",
      "mode3_pre_e1", "mode3_pre_e2")
  } else {
    c("mode1_post_e1", "mode1_post_e2",
      "mode2_post_e1", "mode2_post_e2",
      "mode3_post_e1", "mode3_post_e2")
  }
  
  plot_df <- error_table %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols_use),
      names_to = c("mode", "period", "strategy"),
      names_pattern = "mode([123])_(pre|post)_(e[12])",
      values_to = "error"
    ) %>%
    dplyr::mutate(
      mode = factor(mode, levels = c("1", "2", "3"),
                    labels = c("Mode 1", "Mode 2", "Mode 3")),
      strategy = factor(strategy, levels = c("e1", "e2"),
                        labels = c("Direct-est", "Re-est")),
      x_num = dplyr::if_else(strategy == "Direct-est", 1, 2),
      setting = paste0("T=", Time, "\n", dim_obs)
    )
  
  line_df <- error_table %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols_use),
      names_to = c("mode", "period", "strategy"),
      names_pattern = "mode([123])_(pre|post)_(e[12])",
      values_to = "error"
    ) %>%
    tidyr::pivot_wider(
      names_from = strategy,
      values_from = error
    ) %>%
    dplyr::mutate(
      mode = factor(mode, levels = c("1", "2", "3"),
                    labels = c("Mode 1", "Mode 2", "Mode 3")),
      setting = paste0("T=", Time, "\n", dim_obs),
      line_col = dplyr::case_when(
        e1 > e2 ~ "Improved",
        e1 < e2 ~ "Worse",
        TRUE ~ "Equal"
      )
    )
  
  ggplot() +
    geom_segment(
      data = line_df,
      aes(
        x = 1, xend = 2,
        y = e1, yend = e2,
        group = setting,
        colour = line_col
      ),
      linewidth = 0.8,
      alpha = 0.8
    ) +
    geom_point(
      data = plot_df,
      aes(x = x_num, y = error, colour = strategy),
      size = 2.4
    ) +
    facet_wrap(~ mode, scales = "fixed") +
    #facet_wrap(~ mode, scales = "free_y") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Direct-est", "Re-est")
    ) +
    scale_colour_manual(
      values = c(
        "Direct-est" = "#1f78b4",
        "Re-est" = "#33a02c",
        "Improved" = "red",
        "Worse" = "blue",
        "Equal" = "grey50"
      ),
      breaks = c("Direct-est", "Re-est", "Improved", "Worse", "Equal")
    ) +
    labs(
      x = NULL,
      y = paste0(
        if (period == "pre") "Pre-change" else "Post-change",
        " loading estimation error"
      ),
      colour = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}