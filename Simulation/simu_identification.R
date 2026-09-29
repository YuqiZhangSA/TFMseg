min_len_delta <- function(cp_vec, j, st = 0L, ed, p) {
  cpj <- cp_vec[j]
  left  <- if (j == 1L) st else cp_vec[j - 1L]
  right <- if (j == length(cp_vec)) ed else cp_vec[j + 1L]
  Lpre  <- cpj - left
  Lpost <- right - cpj
  minLen <- max(1L, as.integer(min(Lpre, Lpost)))
  delta_minLenp <- (1 / sqrt(minLen)) + (1 / p)
  list(
    minLen = minLen,
    delta_minLenp = delta_minLenp,
    Lpre = as.integer(Lpre),
    Lpost = as.integer(Lpost)
  )
}

id_modes_1_one_setting <- function(
    cp_est_df = NULL,
    nrep = 100,
    dim_obs = c(10, 10, 10),
    dim_latent = c(3, 3, 3),
    Time = 400L,
    theta_coef = c(0.25, 0.5, 0.75),
    r_hat = NULL,
    dep = TRUE,
    trim_len = 1,
    use_true_cp = FALSE,
    data_setting = c("s0", "s1", "s1_var1", "s3_3I", "s3_A2", "s3_A1", "s2"),
    idio_coeff = TRUE
) {
  data_setting <- match.arg(data_setting)
  
  if (!isTRUE(use_true_cp) && is.null(cp_est_df)) {
    stop("When use_true_cp = FALSE, cp_est_df must be provided.")
  }
  
  need_fns <- c(
    "make_tensor_setting",
    "dgp_general",
    "global_pca",
    "mode_cov_from_G",
    "zeta_pair",
    "min_len_delta"
  )
  missing_fns <- need_fns[!vapply(need_fns, exists, logical(1), mode = "function")]
  
  if (length(missing_fns)) {
    stop("Missing required function(s): ", paste(missing_fns, collapse = ", "))
  }
  
  seed_start <- 900L
  
  if (!is.numeric(trim_len) || length(trim_len) != 1L ||
      is.na(trim_len) || trim_len <= 0 || trim_len > 1) {
    stop("trim_len must be in (0, 1].")
  }
  
  Time <- as.integer(Time)
  dim_obs <- as.integer(dim_obs)
  dim_latent <- as.integer(dim_latent)
  
  K <- length(dim_obs)
  p <- prod(dim_obs)
  
  theta <- floor(Time * theta_coef)
  theta <- sort(unique(theta[theta > 0L & theta < Time]))
  
  if (data_setting == "s0") {
    theta_eval <- as.integer(floor(Time / 2))
    q_true <- 1L
  } else {
    theta_eval <- theta
    q_true <- length(theta_eval)
    
    if (q_true < 1L) {
      stop("Need at least one true cp: provide non-empty theta_coef.")
    }
  }
  
  delta_Tp <- (1 / sqrt(Time)) + (1 / p)
  
  if (isTRUE(use_true_cp)) {
    sims <- seq_len(nrep)
  } else {
    need_cols <- c("sim", "scaled_cp", "Time", "dim_obs")
    miss <- setdiff(need_cols, names(cp_est_df))
    
    if (length(miss)) {
      stop("cp_est_df missing columns: ", paste(miss, collapse = ", "))
    }
    
    dim_obs_chr <- paste(dim_obs, collapse = "x")
    
    cp_est_df <- cp_est_df[
      cp_est_df$Time == Time & cp_est_df$dim_obs == dim_obs_chr,
      ,
      drop = FALSE
    ]
    
    if (!nrow(cp_est_df)) {
      stop("After filtering by Time and dim_obs, cp_est_df has 0 rows.")
    }
    
    sims <- sort(unique(cp_est_df$sim))
    sims <- sims[sims >= 1L & sims <= nrep]
    
    if (!length(sims)) {
      stop("No valid sims found after filtering.")
    }
  }
  
  zeta_by_cp <- vector("list", q_true)
  for (j in seq_len(q_true)) {
    zeta_by_cp[[j]] <- list()
  }
  
  min_window <- 2L
  
  for (sim_id in sims) {
    if (isTRUE(use_true_cp)) {
      cp_vec <- theta_eval
    } else {
      df_s <- cp_est_df[cp_est_df$sim == sim_id, , drop = FALSE]
      
      cp_vec <- as.integer(round(Time * df_s$scaled_cp))
      cp_vec <- sort(unique(cp_vec))
      cp_vec <- cp_vec[cp_vec > 1L & cp_vec < Time]
      
      if (length(cp_vec) >= q_true) {
        cp_vec <- cp_vec[seq_len(q_true)]
      }
    }
    
    if (length(cp_vec) < q_true) next
    
    set.seed(seed_start + sim_id)
    
    setting_args <- make_tensor_setting(
      data_setting = data_setting,
      theta = theta,
      dim_latent = dim_latent
    )
    
    X <- dgp_general(
      model = "tensor",
      dim_obs = dim_obs,
      dim_latent = dim_latent,
      Time = Time,
      dep = dep,
      coeff = 0.7,
      true_cp = setting_args$true_cp,
      type_mode = setting_args$type_mode,
      type_change = setting_args$type_change,
      shift_ind = setting_args$shift_ind,
      shift_mean = setting_args$shift_mean,
      shift_var = setting_args$shift_var,
      transform_list = setting_args$transform_list,
      add_factors = setting_args$add_factors,
      add_factors_coeff = setting_args$add_factors_coeff,
      dist = "Gaussian",
      idio_coeff = idio_coeff
    )$X
    
    est <- if (is.null(r_hat)) {
      global_pca(X = X, dim_X = dim_obs, centre = TRUE, proj = TRUE)
    } else {
      global_pca(X = X, dim_X = dim_obs, r_hat = r_hat, centre = TRUE, proj = TRUE)
    }
    
    G <- est$G_proj
    
    for (j in seq_len(q_true)) {
      cpj <- cp_vec[j]
      
      left <- if (j == 1L) 0L else cp_vec[j - 1L]
      right <- if (j == length(cp_vec)) Time else cp_vec[j + 1L]
      
      Lpre_avail <- cpj - left
      Lpost_avail <- right - cpj
      
      Lpre <- max(min_window, as.integer(floor(trim_len * Lpre_avail)))
      Lpost <- max(min_window, as.integer(floor(trim_len * Lpost_avail)))
      
      a_pre <- max(left, cpj - Lpre)
      b_pre <- cpj
      a_post <- cpj
      b_post <- min(right, cpj + Lpost)
      
      if (b_pre <= a_pre || b_post <= a_post) next
      
      len_info <- min_len_delta(
        cp_vec,
        j,
        st = 0L,
        ed = Time,
        p = p
      )
      
      for (k in seq_len(K)) {
        G_pre <- mode_cov_from_G(
          G,
          mode_k = k,
          a = a_pre,
          b = b_pre
        )
        
        G_post <- mode_cov_from_G(
          G,
          mode_k = k,
          a = a_post,
          b = b_post
        )
        
        z <- zeta_pair(G_pre, G_post)
        
        zeta_by_cp[[j]][[length(zeta_by_cp[[j]]) + 1L]] <- data.frame(
          sim = sim_id,
          cp_index = j,
          cp_time_hat = cp_vec[j],
          mode = k,
          zeta = z,
          ratio_Tp = z / delta_Tp,
          ratio_minLenp = z / len_info$delta_minLenp,
          ratio_logT = z / sqrt(log(Time)),
          minLen = len_info$minLen,
          Lpre = len_info$Lpre,
          Lpost = len_info$Lpost,
          window_pre = sprintf("(%d,%d]", a_pre, b_pre),
          window_post = sprintf("(%d,%d]", a_post, b_post),
          dim = paste(dim_obs, collapse = "x"),
          data_setting = data_setting,
          idio_coeff = idio_coeff,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  zeta_by_cp <- lapply(
    zeta_by_cp,
    function(lst) {
      if (length(lst)) do.call(rbind, lst) else data.frame()
    }
  )
  
  list(
    zeta_by_cp = zeta_by_cp,
    meta = list(
      setting = data_setting,
      Time = Time,
      dim_obs = dim_obs,
      dim_latent = dim_latent,
      p = p,
      theta_coef = theta_coef,
      true_cp = theta,
      theta_eval = theta_eval,
      delta_Tp = delta_Tp,
      delta_def = "delta_Tp = 1/sqrt(T) + 1/p; delta_minLenp(j) = 1/sqrt(minLen_j) + 1/p",
      trim_len = trim_len,
      use_true_cp = use_true_cp,
      dep = dep,
      idio_coeff = idio_coeff
    )
  )
}


id_modes_1 <- function(
    cp_est_df = NULL,
    nrep = 100,
    dim_obs_list = NULL,
    Time_list = NULL,
    dim_latent = c(3, 3, 3),
    theta_coef = c(0.25, 0.5, 0.75),
    r_hat = NULL,
    dep = TRUE,
    trim_len = 1,
    use_true_cp = FALSE,
    data_setting = c("s0", "s1", "s1_var1", "s3_3I", "s3_A2", "s3_A1", "s2"),
    idio_coeff = TRUE
) {
  data_setting <- match.arg(data_setting)
  
  if (!isTRUE(use_true_cp) && is.null(cp_est_df)) {
    stop("When use_true_cp = FALSE, cp_est_df must be provided.")
  }
  
  if (!isTRUE(use_true_cp)) {
    need_cols <- c("sim", "scaled_cp", "Time", "dim_obs")
    miss <- setdiff(need_cols, names(cp_est_df))
    
    if (length(miss)) {
      stop("cp_est_df missing columns: ", paste(miss, collapse = ", "))
    }
  }
  
  if (!is.null(dim_obs_list) && !isTRUE(use_true_cp)) {
    dim_obs_chr_list <- vapply(
      dim_obs_list,
      function(d) paste(as.integer(d), collapse = "x"),
      character(1)
    )
    
    cp_est_df <- cp_est_df[
      cp_est_df$dim_obs %in% dim_obs_chr_list,
      ,
      drop = FALSE
    ]
  }
  
  if (!is.null(Time_list) && !isTRUE(use_true_cp)) {
    cp_est_df <- cp_est_df[
      cp_est_df$Time %in% Time_list,
      ,
      drop = FALSE
    ]
  }
  
  if (isTRUE(use_true_cp)) {
    if (is.null(dim_obs_list) || is.null(Time_list)) {
      stop("When use_true_cp = TRUE, provide dim_obs_list and Time_list.")
    }
    
    settings <- expand.grid(
      Time = as.integer(Time_list),
      dim_obs = vapply(
        dim_obs_list,
        function(d) paste(as.integer(d), collapse = "x"),
        character(1)
      ),
      stringsAsFactors = FALSE
    )
    
    settings <- settings[
      order(settings$Time, settings$dim_obs),
      ,
      drop = FALSE
    ]
  } else {
    settings <- unique(cp_est_df[, c("Time", "dim_obs")])
    
    settings <- settings[
      order(settings$Time, settings$dim_obs),
      ,
      drop = FALSE
    ]
  }
  
  out <- vector("list", nrow(settings))
  
  for (i in seq_len(nrow(settings))) {
    Time_i <- as.integer(settings$Time[i])
    dim_chr <- settings$dim_obs[i]
    dim_i <- as.integer(strsplit(dim_chr, "x", fixed = TRUE)[[1]])
    
    message(sprintf(
      "[id_modes_1] %d/%d: Time=%d, dim_obs=%s, setting=%s, idio_coeff=%s",
      i,
      nrow(settings),
      Time_i,
      dim_chr,
      data_setting,
      idio_coeff
    ))
    
    out[[i]] <- id_modes_1_one_setting(
      cp_est_df = if (isTRUE(use_true_cp)) NULL else cp_est_df,
      nrep = nrep,
      dim_obs = dim_i,
      dim_latent = dim_latent,
      Time = Time_i,
      theta_coef = theta_coef,
      r_hat = r_hat,
      dep = dep,
      trim_len = trim_len,
      use_true_cp = use_true_cp,
      data_setting = data_setting,
      idio_coeff = idio_coeff
    )
    
    out[[i]]$meta$dim_obs_chr <- dim_chr
  }
  
  names(out) <- paste0(
    "T=", settings$Time,
    "|dim=", settings$dim_obs,
    "|set=", data_setting,
    "|idio_coeff=", idio_coeff
  )
  
  out
}



#------------------------------------------------------------
# plots for each cp
#------------------------------------------------------------

library(dplyr)
library(ggplot2)

plot_zeta_mode <- function(df,
                      y = c("zeta", "ratio_Tp", "ratio_minLenp"),
                      x = "sim",
                      colour = "mode",
                      keep_finite = TRUE,
                      low_mode = NULL,
                      low_bound = c("min", "q01", "q05", "q10", "q25"),
                      up_mode = NULL,
                      up_bound = c("max", "q99", "q95", "q90", "q75")) {
  y <- match.arg(y)
  low_bound <- match.arg(low_bound)
  up_bound  <- match.arg(up_bound)
  
  stopifnot(is.data.frame(df))
  if (!("mode" %in% names(df))) stop("df must contain a column named 'mode'.")
  if (!y %in% names(df)) stop("Column not found: y = ", y)
  if (!x %in% names(df)) stop("Column not found: x = ", x)
  
  dat <- df
  if (keep_finite) dat <- dat %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  
  get_bound <- function(vec, which_bound) {
    switch(
      which_bound,
      min = min(vec, na.rm = TRUE),
      max = max(vec, na.rm = TRUE),
      q01 = as.numeric(stats::quantile(vec, probs = 0.01, na.rm = TRUE)),
      q05 = as.numeric(stats::quantile(vec, probs = 0.05, na.rm = TRUE)),
      q10 = as.numeric(stats::quantile(vec, probs = 0.10, na.rm = TRUE)),
      q25 = as.numeric(stats::quantile(vec, probs = 0.25, na.rm = TRUE)),
      q75 = as.numeric(stats::quantile(vec, probs = 0.75, na.rm = TRUE)),
      q90 = as.numeric(stats::quantile(vec, probs = 0.90, na.rm = TRUE)),
      q95 = as.numeric(stats::quantile(vec, probs = 0.95, na.rm = TRUE)),
      q99 = as.numeric(stats::quantile(vec, probs = 0.99, na.rm = TRUE))
    )
  }
  
  lb <- NA_real_
  if (!is.null(low_mode)) {
    dat_low <- dat %>% filter(mode == low_mode)
    if (!nrow(dat_low)) stop("No rows for low_mode == ", low_mode)
    lb <- get_bound(dat_low[[y]], low_bound)
  }
  
  ub <- NA_real_
  if (!is.null(up_mode)) {
    dat_up <- dat %>% filter(mode == up_mode)
    if (!nrow(dat_up)) stop("No rows for up_mode == ", up_mode)
    ub <- get_bound(dat_up[[y]], up_bound)
  }
  
  p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(colour) && colour %in% names(dat)) {
    p <- p + geom_point(aes(colour = factor(.data[[colour]])), alpha = 0.6)
  } else {
    p <- p + geom_point(alpha = 0.6)
  }
  
  if (is.finite(lb)) {
    p <- p +
      geom_hline(yintercept = lb, linewidth = 0.9, linetype = 2) +
      annotate("text",
               x = Inf, y = lb,
               label = paste0("lower: mode ", low_mode, " (", low_bound, ") = ", signif(lb, 4)),
               hjust = 1.05, vjust = -0.4)
  }
  
  if (is.finite(ub)) {
    p <- p +
      geom_hline(yintercept = ub, linewidth = 0.9, linetype = 3) +
      annotate("text",
               x = -Inf, y = ub,
               label = paste0("upper: mode ", up_mode, " (", up_bound, ") = ", signif(ub, 4)),
               hjust = -0.05, vjust = -0.4)
  }
  
  p + labs(x = x, y = y, colour = if (!is.null(colour)) colour else NULL)
}



pool_by_cp <- function(res_all_true) {
  q <- length(res_all_true[[1]]$zeta_by_cp)
  out <- vector("list", q)
  
  for (j in seq_len(q)) {
    dfs_j <- list()
    ii <- 0L
    
    for (nm in names(res_all_true)) {
      obj <- res_all_true[[nm]]
      dfj <- obj$zeta_by_cp[[j]]
      if (is.null(dfj) || !nrow(dfj)) next
      
      ii <- ii + 1L
      dfj$Time <- obj$meta$Time
      dfj$dim_obs <- paste(obj$meta$dim_obs, collapse = "x")
      dfj$setting_key <- nm
      dfs_j[[ii]] <- dfj
    }
    
    out[[j]] <- if (length(dfs_j)) do.call(rbind, dfs_j) else data.frame()
  }
  
  out
}

#------------------------------------------------------------
# plots for pooled all together
#------------------------------------------------------------
#------------------------------------------------------------
# 1. True changed modes by setting
#------------------------------------------------------------
get_changed_modes <- function(data_setting, cp_index) {
  data_setting <- match.arg(
    data_setting,
    choices = c("s2", "s1", "s0", "s3_3I", "s3_A2", "s3_A1")
  )
  
  if (data_setting == "s1") {
    # cp1 -> mode 1, cp2 -> mode 2, cp3 -> mode 3
    out <- list(
      `1` = 1L,
      `2` = 2L,
      `3` = 3L
    )
    return(out[[as.character(cp_index)]])
  }
  
  if (data_setting == "s2") {
    out <- list(
      `1` = 1L,
      `2` = 2L,
      `3` = c(2L, 3L)
    )
    return(out[[as.character(cp_index)]])
  }
  
  if (data_setting %in% c("s3_3I", "s3_A2", "s3_A1")) {
    return(1L)
  }
  
  if (data_setting == "s0") {
    return(integer(0))
  }
  
  stop("Unknown data_setting: ", data_setting)
}

#------------------------------------------------------------
# 2. Pool all cps together and label changed / unchanged
#------------------------------------------------------------
pool_all <- function(res_all_true) {
  dfs <- list()
  ii <- 0L
  
  for (nm in names(res_all_true)) {
    obj <- res_all_true[[nm]]
    setting <- obj$meta$setting
    
    for (j in seq_along(obj$zeta_by_cp)) {
      dfj <- obj$zeta_by_cp[[j]]
      if (is.null(dfj) || !nrow(dfj)) next
      
      changed_modes_j <- get_changed_modes(setting, j)
      
      dfj <- dfj %>%
        mutate(
          Time = obj$meta$Time,
          dim_obs = paste(obj$meta$dim_obs, collapse = "x"),
          setting_key = nm,
          changed_truth = if_else(mode %in% changed_modes_j, "changed", "unchanged"),
          cp_label = paste0("cp", j)
        )
      
      ii <- ii + 1L
      dfs[[ii]] <- dfj
    }
  }
  
  if (!length(dfs)) return(data.frame())
  
  out <- bind_rows(dfs) %>%
    mutate(
      changed_truth = factor(changed_truth, levels = c("unchanged", "changed")),
      row_id = seq_len(n())
    )
  
  out
}

#------------------------------------------------------------
# 3. Pooled scatter plot
#------------------------------------------------------------
# plot_zeta_pooled <- function(df,
#                              y = c("zeta", "ratio_Tp", "ratio_minLenp"),
#                              x = c("row_id", "sim"),
#                              keep_finite = TRUE,
#                              up_bound = c("q75", "q90", "q95", "q99", "max"),
#                              point_alpha = 0.55,
#                              point_size = 1.8,
#                              facet_by = c("none", "Time", "dim_obs", "cp_label")) {
#   y <- match.arg(y)
#   x <- match.arg(x)
#   up_bound <- match.arg(up_bound)
#   facet_by <- match.arg(facet_by)
#   
#   stopifnot(is.data.frame(df))
#   req_cols <- c("changed_truth", y, x)
#   miss <- setdiff(req_cols, names(df))
#   if (length(miss)) {
#     stop("df is missing columns: ", paste(miss, collapse = ", "))
#   }
#   
#   dat <- df
#   if (keep_finite) {
#     dat <- dat %>%
#       dplyr::filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
#   }
#   
#   dat_unch <- dat %>% dplyr::filter(changed_truth == "unchanged")
#   if (!nrow(dat_unch)) stop("No unchanged rows found.")
#   
#   get_bound <- function(vec, which_bound) {
#     switch(
#       which_bound,
#       max = max(vec, na.rm = TRUE),
#       q75 = as.numeric(stats::quantile(vec, probs = 0.75, na.rm = TRUE)),
#       q90 = as.numeric(stats::quantile(vec, probs = 0.90, na.rm = TRUE)),
#       q95 = as.numeric(stats::quantile(vec, probs = 0.95, na.rm = TRUE)),
#       q99 = as.numeric(stats::quantile(vec, probs = 0.99, na.rm = TRUE))
#     )
#   }
#   
#   ub <- get_bound(dat_unch[[y]], up_bound)
#   
#   p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]], colour = changed_truth)) +
#     geom_point(alpha = point_alpha, size = point_size) +
#     geom_hline(yintercept = ub, linewidth = 0.9, linetype = 2) +
#     annotate(
#       "text",
#       x = -Inf, y = ub,
#       label = paste0("upper bound (", up_bound, ") = ", signif(ub, 4)),
#       hjust = -0.05, vjust = -0.4
#     ) +
#     labs(
#       x = NULL,
#       y = expression(Xi/(1/sqrt(T) + 1/p)),
#       colour = NULL
#     ) +
#     theme_classic() +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       legend.position = "bottom",
#       legend.text = element_text(size = 12),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank()
#     )
#   
#   if (facet_by != "none") {
#     p <- p + facet_wrap(stats::as.formula(paste("~", facet_by)), scales = "free_x")
#   }
#   
#   p
# }


plot_zeta_pooled <- function(df,
                             y = c("zeta", "ratio_Tp", "ratio_minLenp"),
                             x = c("row_id", "sim"),
                             keep_finite = TRUE,
                             up_bound = c("q75", "q90", "q95", "q99", "max"),
                             point_alpha = 0.55,
                             point_size = 1.8,
                             show_setting_boundaries = TRUE,
                             show_dim_top = TRUE,
                             show_time_bracket = TRUE,
                             dim_text_size = 2.8,
                             time_text_size = 2.8,
                             ub_text_size = 3.2,
                             axis_text_size = 9,
                             axis_title_size = 9,
                             legend_text_size = 10,
                             legend_key_size = 1.0) {
  y <- match.arg(y)
  x <- match.arg(x)
  up_bound <- match.arg(up_bound)
  
  stopifnot(is.data.frame(df))
  req_cols <- c("changed_truth", y, x, "Time", "dim_obs")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) {
    stop("df is missing columns: ", paste(miss, collapse = ", "))
  }
  
  dat <- df
  
  if (keep_finite) {
    dat <- dat %>%
      dplyr::filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  }
  
  dat_unch <- dat %>% dplyr::filter(changed_truth == "unchanged")
  if (!nrow(dat_unch)) stop("No unchanged rows found.")
  
  get_bound <- function(vec, which_bound) {
    switch(
      which_bound,
      max = max(vec, na.rm = TRUE),
      q75 = as.numeric(stats::quantile(vec, probs = 0.75, na.rm = TRUE)),
      q90 = as.numeric(stats::quantile(vec, probs = 0.90, na.rm = TRUE)),
      q95 = as.numeric(stats::quantile(vec, probs = 0.95, na.rm = TRUE)),
      q99 = as.numeric(stats::quantile(vec, probs = 0.99, na.rm = TRUE))
    )
  }
  
  ub <- get_bound(dat_unch[[y]], up_bound)
  
  setting_df <- dat %>%
    dplyr::group_by(Time, dim_obs) %>%
    dplyr::summarise(
      xmin = min(.data[[x]], na.rm = TRUE),
      xmax = max(.data[[x]], na.rm = TRUE),
      xmid = (xmin + xmax) / 2,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dim_obs_lab = paste0("(", gsub("x", ",", dim_obs), ")")
    ) %>%
    dplyr::arrange(xmin)
  
  time_df <- setting_df %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(
      xmin = min(xmin),
      xmax = max(xmax),
      xmid = (xmin + xmax) / 2,
      .groups = "drop"
    ) %>%
    dplyr::arrange(xmin)
  
  y_min <- min(dat[[y]], na.rm = TRUE)
  y_max <- max(dat[[y]], na.rm = TRUE)
  y_rng <- y_max - y_min
  
  y_top_lab      <- y_max + 0.04 * y_rng
  y_bot_bracket  <- y_min - 0.010 * y_rng
  y_bot_text     <- y_min - 0.050 * y_rng
  bracket_tick_h <- 0.030 * y_rng
  
  p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]], colour = changed_truth)) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_hline(yintercept = ub, linewidth = 0.9, linetype = 2) +
    annotate(
      "text",
      x = min(dat[[x]], na.rm = TRUE),
      y = ub,
      label = paste0("upper bound (", up_bound, ") = ", signif(ub, 4)),
      hjust = 0,
      vjust = -0.35,
      size = ub_text_size
    ) +
    labs(
      x = NULL,
      y = expression(Xi / (1/sqrt(T) + 1/p)),
      colour = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = axis_text_size),
      axis.title.y = element_text(size = axis_title_size),
      legend.position = "bottom",
      legend.text = element_text(size = legend_text_size),
      legend.key.size = grid::unit(legend_key_size, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 28, r = 12, b = 12, l = 12)
    ) +
    coord_cartesian(
      ylim = c(y_min - 0.060 * y_rng, y_max + 0.060 * y_rng),
      clip = "off"
    )
  
  if (show_setting_boundaries && nrow(setting_df) > 1) {
    boundary_x <- head(setting_df$xmax, -1) + 0.5
    p <- p +
      geom_vline(
        xintercept = boundary_x,
        linetype = 3,
        linewidth = 0.4,
        colour = "grey60"
      )
  }
  
  if (show_dim_top) {
    p <- p +
      geom_text(
        data = setting_df,
        aes(x = xmid, y = y_top_lab, label = dim_obs_lab),
        inherit.aes = FALSE,
        angle = 90,
        hjust = 1,      # align by right bracket at the top
        vjust = 0.5,
        size = dim_text_size
      )
  }
  
  if (show_time_bracket) {
    p <- p +
      geom_segment(
        data = time_df,
        aes(x = xmin, xend = xmax, y = y_bot_bracket, yend = y_bot_bracket),
        inherit.aes = FALSE,
        linewidth = 0.45,
        colour = "grey35"
      ) +
      geom_segment(
        data = time_df,
        aes(x = xmin, xend = xmin,
            y = y_bot_bracket - bracket_tick_h / 2,
            yend = y_bot_bracket + bracket_tick_h / 2),
        inherit.aes = FALSE,
        linewidth = 0.45,
        colour = "grey35"
      ) +
      geom_segment(
        data = time_df,
        aes(x = xmax, xend = xmax,
            y = y_bot_bracket - bracket_tick_h / 2,
            yend = y_bot_bracket + bracket_tick_h / 2),
        inherit.aes = FALSE,
        linewidth = 0.45,
        colour = "grey35"
      ) +
      geom_text(
        data = time_df,
        aes(x = xmid, y = y_bot_text, label = paste0("T = ", Time)),
        inherit.aes = FALSE,
        size = time_text_size,
        colour = "grey20"
      )
  }
  
  p
}