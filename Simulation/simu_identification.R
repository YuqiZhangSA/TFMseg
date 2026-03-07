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
    data_setting = c("s2", "s1", "s0", "s3_3I", "s3_A2", "s3_A1")
) {
  data_setting <- match.arg(data_setting)
  
  if (!isTRUE(use_true_cp) && is.null(cp_est_df)) {
    stop("When use_true_cp = FALSE, cp_est_df must be provided.")
  }
  
  seed_start <- 900
  
  if (!is.numeric(trim_len) || trim_len <= 0 || trim_len > 1) {
    stop("trim_len must be in (0,1].")
  }
  
  Time <- as.integer(Time)
  dim_obs <- as.integer(dim_obs)
  dim_latent <- as.integer(dim_latent)
  
  K <- length(dim_obs)
  p <- prod(dim_obs)
  
  theta <- floor(Time * theta_coef)
  theta <- sort(unique(theta[theta > 0L & theta < Time]))
  
  if (data_setting == "s0") {
    theta <- as.integer(floor(Time / 2))
    q_true <- 1L
  } else {
    q_true <- length(theta)
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
  for (j in seq_len(q_true)) zeta_by_cp[[j]] <- list()
  
  min_window <- 2L
  
  for (sim_id in sims) {
    if (isTRUE(use_true_cp)) {
      cp_vec <- theta
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
    r0 <- 3L
    
    if (data_setting == "s1") {
      A3 <- matrix(rnorm(r0^2, sd = 1 / sqrt(r0)), nrow = r0)
      A3[lower.tri(A3)] <- t(A3)[lower.tri(A3)]
      
      A2 <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 0
      ), nrow = 3, byrow = TRUE)
      
      A1 <- matrix(0, 3, 3)
      A1[1, 1] <- 0.5
      A1[2, 1] <- rnorm(1, 0, 1); A1[2, 2] <- 1
      A1[3, 1] <- rnorm(1, 0, 1); A1[3, 2] <- rnorm(1, 0, 1); A1[3, 3] <- 1.5
      
      transform_list <- list(
        list(A1, diag(3), diag(3)),
        list(diag(3), A2, diag(3)),
        list(diag(3), diag(3), A3)
      )
      type_mode <- list(1, 2, 3)
      type_change <- list("l", "l", "l")
      true_cp <- theta
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else if (data_setting == "s3_3I") {
      stopifnot(length(theta) == 1L)
      
      r1 <- dim_latent[1]
      transform_list <- list(
        list(3 * diag(r1), diag(dim_latent[2]), diag(dim_latent[3]))
      )
      type_mode <- list(1)
      type_change <- list("l")
      true_cp <- theta
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else if (data_setting == "s3_A2") {
      stopifnot(length(theta) == 1L)
      stopifnot(all(dim_latent >= 3))
      
      A2 <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 0
      ), nrow = 3, byrow = TRUE)
      
      transform_list <- list(
        list(A2, diag(dim_latent[2]), diag(dim_latent[3]))
      )
      type_mode <- list(1)
      type_change <- list("l")
      true_cp <- theta
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else if (data_setting == "s3_A1") {
      stopifnot(length(theta) == 1L)
      stopifnot(all(dim_latent >= 3))
      
      A1 <- matrix(0, 3, 3)
      A1[1, 1] <- 0.5
      A1[2, 1] <- rnorm(1, 0, 1); A1[2, 2] <- 1
      A1[3, 1] <- rnorm(1, 0, 1); A1[3, 2] <- rnorm(1, 0, 1); A1[3, 3] <- 1.5
      
      transform_list <- list(
        list(A1, diag(dim_latent[2]), diag(dim_latent[3]))
      )
      type_mode <- list(1)
      type_change <- list("l")
      true_cp <- theta
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else if (data_setting == "s2") {
      if (length(theta) != 3L) {
        stop("For s2 set theta_coef to length 3.")
      }
      stopifnot(length(dim_latent) == 3L)
      stopifnot(all(dim_latent >= 3))
      
      A3_3 <- matrix(rnorm(r0^2, sd = 1 / sqrt(r0)), nrow = r0)
      A3_3[lower.tri(A3_3)] <- t(A3_3)[lower.tri(A3_3)]
      
      A2_3 <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 0
      ), nrow = 3, byrow = TRUE)
      
      A1_3 <- matrix(0, 3, 3)
      A1_3[1, 1] <- 0.5
      A1_3[2, 1] <- rnorm(1, 0, 1); A1_3[2, 2] <- 1
      A1_3[3, 1] <- rnorm(1, 0, 1); A1_3[3, 2] <- rnorm(1, 0, 1); A1_3[3, 3] <- 1.5
      
      A3_mode3 <- diag(dim_latent[3])
      A3_mode3[1:3, 1:3] <- A3_3
      
      A2_mode2 <- diag(dim_latent[2])
      A2_mode2[1:3, 1:3] <- A2_3
      
      A1_mode1 <- diag(dim_latent[1])
      A1_mode1[1:3, 1:3] <- A1_3
      
      I3_mode3 <- diag(dim_latent[3])
      I3_mode3[1:3, 1:3] <- 3 * diag(3)
      
      vals <- pmax(1 - 0.4 * (0:(dim_latent[2] - 1)), 0.1)
      M2_mode2 <- diag(vals)
      
      transform_list <- list(
        list(A1_mode1, diag(dim_latent[2]), I3_mode3),
        list(diag(dim_latent[1]), A2_mode2, diag(dim_latent[3])),
        list(diag(dim_latent[1]), M2_mode2, A3_mode3)
      )
      type_mode <- list(c(1, 3), 2, c(2, 3))
      type_change <- list("l", "l", "l")
      true_cp <- theta
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else if (data_setting == "s0") {
      transform_list <- NULL
      type_mode <- list()
      type_change <- list()
      true_cp <- integer(0)
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else {
      stop("Unknown data_setting: ", data_setting)
    }
    
    X <- dgp_general(
      model = "tensor",
      dim_obs = dim_obs,
      dim_latent = dim_latent,
      Time = Time,
      dep = dep,
      coeff = 0.7,
      true_cp = true_cp,
      type_mode = type_mode,
      type_change = type_change,
      shift_ind = shift_ind,
      shift_mean = shift_mean,
      shift_var = shift_var,
      transform_list = transform_list,
      add_factors = add_factors,
      add_factors_coeff = add_factors_coeff,
      dist = "Gaussian"
    )$X
    
    est <- if (is.null(r_hat)) {
      global_pca(X = X, dim_X = dim_obs, centre = TRUE, proj = TRUE)
    } else {
      global_pca(X = X, dim_X = dim_obs, r_hat = r_hat, centre = TRUE, proj = TRUE)
    }
    
    G <- est$G_proj
    
    for (j in seq_len(q_true)) {
      cpj <- cp_vec[j]
      left  <- if (j == 1L) 0L else cp_vec[j - 1L]
      right <- if (j == length(cp_vec)) Time else cp_vec[j + 1L]
      
      Lpre_avail  <- cpj - left
      Lpost_avail <- right - cpj
      
      Lpre  <- max(min_window, as.integer(floor(trim_len * Lpre_avail)))
      Lpost <- max(min_window, as.integer(floor(trim_len * Lpost_avail)))
      
      a_pre  <- max(left, cpj - Lpre)
      b_pre  <- cpj
      a_post <- cpj
      b_post <- min(right, cpj + Lpost)
      
      if (b_pre <= a_pre || b_post <= a_post) next
      
      len_info <- min_len_delta(cp_vec, j, st = 0L, ed = Time, p = p)
      
      for (k in seq_len(K)) {
        G_pre  <- mode_cov_from_G(G, mode_k = k, a = a_pre,  b = b_pre)
        G_post <- mode_cov_from_G(G, mode_k = k, a = a_post, b = b_post)
        
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
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  zeta_by_cp <- lapply(
    zeta_by_cp,
    function(lst) if (length(lst)) do.call(rbind, lst) else data.frame()
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
      delta_Tp = delta_Tp,
      delta_def = "delta_Tp = 1/sqrt(T) + 1/p; delta_minLenp(j) = 1/sqrt(minLen_j) + 1/p",
      trim_len = trim_len,
      use_true_cp = use_true_cp
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
    trim_len = 1,
    use_true_cp = FALSE,
    data_setting = c("s2", "s1", "s0", "s3_3I", "s3_A2", "s3_A1")
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
    cp_est_df <- cp_est_df[cp_est_df$dim_obs %in% dim_obs_chr_list, , drop = FALSE]
  }
  
  if (!is.null(Time_list) && !isTRUE(use_true_cp)) {
    cp_est_df <- cp_est_df[cp_est_df$Time %in% Time_list, , drop = FALSE]
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
    settings <- settings[order(settings$Time, settings$dim_obs), , drop = FALSE]
  } else {
    settings <- unique(cp_est_df[, c("Time", "dim_obs")])
    settings <- settings[order(settings$Time, settings$dim_obs), , drop = FALSE]
  }
  
  out <- vector("list", nrow(settings))
  
  for (i in seq_len(nrow(settings))) {
    Time_i <- as.integer(settings$Time[i])
    dim_chr <- settings$dim_obs[i]
    dim_i <- as.integer(strsplit(dim_chr, "x", fixed = TRUE)[[1]])
    
    message(sprintf(
      "[id_modes_1] %d/%d: Time=%d, dim_obs=%s, setting=%s",
      i, nrow(settings), Time_i, dim_chr, data_setting
    ))
    
    out[[i]] <- id_modes_1_one_setting(
      cp_est_df = if (isTRUE(use_true_cp)) NULL else cp_est_df,
      nrep = nrep,
      dim_obs = dim_i,
      dim_latent = dim_latent,
      Time = Time_i,
      theta_coef = theta_coef,
      r_hat = r_hat,
      trim_len = trim_len,
      use_true_cp = use_true_cp,
      data_setting = data_setting
    )
    
    out[[i]]$meta$dim_obs_chr <- dim_chr
  }
  
  nm <- paste0("T=", settings$Time, "|dim=", settings$dim_obs, "|set=", data_setting)
  names(out) <- nm
  out
}





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
