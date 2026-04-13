simu_comparison <- function(method = c("TFMseg", "TFMseg.vec", "FMseg", "LR"),
                            nrep = 100,
                            data_setting = c("s0", "s1", "s1_var1", "s3_3I", "s3_A2", "s3_A1", "s2"),
                            dim_obs,
                            dim_latent,
                            Time,
                            dep = TRUE,
                            m = NULL,
                            theta_coef = NULL,
                            r_hat,
                            trim_coef = 1/4,
                            detect_thd = NULL,
                            thd.type = c("fixed", "oracle"),
                            V.diag = TRUE,
                            lrv = TRUE,
                            rvs = FALSE,
                            miss = FALSE) {
  
  method <- match.arg(method)
  thd.type <- match.arg(thd.type)
  data_setting <- match.arg(data_setting)
  
  seed_start <- 900
  
  need_fns <- c("dgp_general", "global_pca", "TFMseg", "bs_LR_globaltrim")
  missing_fns <- need_fns[!vapply(need_fns, exists, logical(1), mode = "function")]
  if (length(missing_fns)) {
    stop("Missing required function(s): ", paste(missing_fns, collapse = ", "))
  }
  
  if (is.null(theta_coef)) theta_coef <- c(0.25, 0.5, 0.75)
  theta <- floor(Time * theta_coef)
  theta <- sort(unique(theta[theta > 0 & theta < Time]))
  
  trim <- floor(trim_coef * (Time / log(Time)))
  
  frequency_table <- matrix(0L, nrow = nrep, ncol = 5L)
  accuracy_count <- if (length(theta)) matrix(0L, nrep, length(theta)) else matrix(NA_real_, nrep, 0L)
  detected_cp_list <- vector("list", nrep)
  cp_est_scaled <- vector("list", nrep)
  time_sec <- numeric(nrep)
  
  for (sim in seq_len(nrep)) {
    set.seed(seed_start + sim)
    
    if (data_setting == "s1") {
      r0 <- 3
      
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
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
      
    } else if (data_setting == "s1_var1") {
      r0 <- 3
      
      A3 <- matrix(rnorm(r0^2, sd = 1), nrow = r0)
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
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
      
    } else if (data_setting == "s0") {
      true_cp <- integer(0)
      type_mode <- list()
      type_change <- list()
      transform_list <- NULL
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
    } else if (data_setting == "s3_3I") {
      stopifnot(length(theta) == 1L)
      r0 <- dim_latent[1]
      
      transform_list <- list(
        list(3 * diag(r0), diag(dim_latent[2]), diag(dim_latent[3]))
      )
      
      type_mode <- list(1)
      type_change <- list("l")
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
      
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
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
      
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
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
      
    } else if (data_setting == "s2") {
      if (length(theta) != 3L) {
        stop("For s2 set theta_coef to length 3.")
      }
      stopifnot(length(dim_latent) == 3L)
      stopifnot(all(dim_latent >= 3))
      
      embed3 <- function(B3, r) {
        M <- diag(r)
        M[1:3, 1:3] <- B3
        M
      }
      
      make_M2 <- function(r) {
        vals <- pmax(1 - 0.4 * (0:(r - 1)), 0.1)
        diag(vals)
      }
      
      r0 <- 3
      
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
      
      A3_mode3 <- embed3(A3_3, dim_latent[3])
      A2_mode2 <- embed3(A2_3, dim_latent[2])
      A1_mode1 <- embed3(A1_3, dim_latent[1])
      I3_mode3 <- embed3(3 * diag(3), dim_latent[3])
      M2_mode2 <- make_M2(dim_latent[2])
      
      transform_list <- list(
        list(A1_mode1, diag(dim_latent[2]), I3_mode3),
        list(diag(dim_latent[1]), A2_mode2, diag(dim_latent[3])),
        list(diag(dim_latent[1]), M2_mode2, A3_mode3)
      )
      
      type_mode <- list(c(1, 3), 2, c(2, 3))
      type_change <- list("l", "l", "l")
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
    } else {
      stop("Unknown data_setting: ", data_setting)
    }
    
    m_true <- length(true_cp)
    
    data_sim <- dgp_general(
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
    )
    
    X <- data_sim$X
    K <- length(dim(X)) - 1L
    dim_X <- dim(X)[seq_len(K)] #(p1,p2,p3)
    
    if (isTRUE(rvs)) X <- X[, , , Time:1, drop = FALSE]
    
    if (method == "TFMseg") {
      start_time <- Sys.time()
      
      if (!isTRUE(miss)) {
        if (is.null(r_hat)) {
          est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
          G <- est_load$G_proj
          G_dim <- as.vector(est_load$r_hat)
        } else {
          est_load <- global_pca(X = X, dim_X = dim_X, r_hat = r_hat, proj = TRUE)
          G <- est_load$G_proj
          G_dim <- as.vector(r_hat)
        }
      } else {
        set.seed(seed_start + sim)
        X_miss <- tensorMiss::miss_gen(aperm(X, c(4, 1, 2, 3)), type = "simul")
        
        if (is.null(r_hat)) {
          est_load <- tensorMiss::miss_factor_est(X_miss, r = 0)
          G <- aperm(est_load$Ft, c(2, 3, 4, 1))
          G_dim <- as.vector(est_load$r)
        } else {
          est_load <- tensorMiss::miss_factor_est(X_miss, r = r_hat)
          G <- aperm(est_load$Ft, c(2, 3, 4, 1))
          G_dim <- as.vector(r_hat)
        }
      }
      
      if (thd.type == "oracle") {
        out <- TFMseg(
          G, G_dim,
          trim = trim,
          method = "oracle",
          m = m_true,
          V.diag = V.diag,
          lrv = lrv
        )
      } else {
        dr <- sum(G_dim * (G_dim + 1L) / 2L)
        
        if (is.null(detect_thd)) {
          thd <- 209.8954613 * sqrt(log(Time)) +
            0.7127491 * sqrt(dr) +
            1566.8875353 * sqrt(1 / log(Time)) +
            1572.7173337 * log(log(Time)) / sqrt(log(Time)) +
            (-2298.3882769)
        } else {
          thd <- detect_thd
        }
        
        out <- TFMseg(
          G, G_dim,
          method = "fixed",
          threshold = thd,
          V.diag = V.diag,
          lrv = lrv
        )
      }
      
      detected_cp <- out$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(lubridate::as.period(end_time - start_time, unit = "sec"))
      
    } else if (method == "TFMseg.vec") {
      start_time <- Sys.time()
      
      if (is.null(r_hat)) {
        est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
        G_core <- est_load$G_proj
        G_dim_vec <- as.vector(est_load$r_hat)
      } else {
        est_load <- global_pca(X = X, dim_X = dim_X, r_hat = r_hat, proj = TRUE)
        G_core <- est_load$G_proj
        G_dim_vec <- as.vector(r_hat)
      }
      
      G <- matrix(G_core, nrow = prod(G_dim_vec), ncol = Time)
      G_dim <- prod(G_dim_vec)
      
      out <- TFMseg(
        G, G_dim,
        trim = trim,
        method = "oracle",
        m = m_true,
        V.diag = V.diag,
        lrv = lrv
      )
      
      detected_cp <- out$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(lubridate::as.period(end_time - start_time, unit = "sec"))
      
    } else if (method == "FMseg") {
      start_time <- Sys.time()
      
      if (is.null(r_hat)) {
        est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
        G_dim <- as.vector(est_load$r_hat)
      } else {
        est_load <- global_pca(X = X, dim_X = dim_X, r_hat = r_hat, proj = TRUE)
        G_dim <- as.vector(r_hat)
      }
      
      Xd <- matrix(aperm(X, c(K + 1L, seq_len(K))), nrow = Time)
      X_vec <- t(Xd)
      
      r_vec <- prod(G_dim)
      lbd <- floor(Time / log(Time))
      
      detected_cp <- FMSeg(
        x = X_vec,
        r = r_vec,
        m = m_true,
        trim = trim,
        method = "oracle",
        lbd = lbd
      )$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(lubridate::as.period(end_time - start_time, unit = "sec"))
      
    } else {
      #source("/Users/yuqi/Documents/GitHub/TFMseg/Simulation/compare_detect_LR.R")
      start_time <- Sys.time()
      
      Xd <- matrix(aperm(X, c(K + 1L, seq_len(K))), nrow = Time)
      min_size <- round(Time / log(Time))
      
      lr_out <- bs_LR_globaltrim(
        Xd,
        tau_global = 0.1, #min(diff(theta_coef))
        min_size = min_size,
        r_est = NULL,
        seed = seed_start + sim
      )
      
      detected_cp <- lr_out$cps %||% integer(0)
      detected_cp <- if (length(detected_cp)) sort(unique(detected_cp)) else integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(lubridate::as.period(end_time - start_time, unit = "sec"))
    }
    
    detected_cp_sorted <- if (length(detected_cp)) sort(detected_cp) else integer(0)
    detected_cp_list[[sim]] <- detected_cp_sorted
    cp_est_scaled[[sim]] <- detected_cp_sorted / Time
    
    m_est_diff <- length(detected_cp_sorted) - m_true
    frequency_table[sim, 1] <- as.integer(m_est_diff <= -2)
    frequency_table[sim, 2] <- as.integer(m_est_diff == -1)
    frequency_table[sim, 3] <- as.integer(m_est_diff == 0)
    frequency_table[sim, 4] <- as.integer(m_est_diff == 1)
    frequency_table[sim, 5] <- as.integer(m_est_diff >= 2)
    
    if (length(theta)) {
      acc_win <- max(1L, as.integer(round(2 * log(Time))))
      for (j in seq_along(theta)) {
        accuracy_count[sim, j] <- as.integer(
          any(abs(detected_cp_sorted - theta[j]) <= acc_win)
        )
      }
    }
  }
  
  freq_summary <- colSums(frequency_table) / nrep
  accuracy_summary <- if (ncol(accuracy_count)) colMeans(accuracy_count) else numeric(0)
  
  rep_raw <- data.frame(
    rep_id = seq_len(nrep),
    le_m2 = frequency_table[, 1],
    m1 = frequency_table[, 2],
    m0 = frequency_table[, 3],
    p1 = frequency_table[, 4],
    ge_p2 = frequency_table[, 5],
    stringsAsFactors = FALSE
  )
  
  if (ncol(accuracy_count) == 0) {
    rep_raw$acc1 <- NA_real_
    rep_raw$acc2 <- NA_real_
    rep_raw$acc3 <- NA_real_
  } else {
    rep_raw$acc1 <- if (ncol(accuracy_count) >= 1) accuracy_count[, 1] else NA_real_
    rep_raw$acc2 <- if (ncol(accuracy_count) >= 2) accuracy_count[, 2] else NA_real_
    rep_raw$acc3 <- if (ncol(accuracy_count) >= 3) accuracy_count[, 3] else NA_real_
  }
  
  rep_raw$time <- time_sec
  
  list(
    cp_est_scaled = cp_est_scaled,
    freq_summary = freq_summary,
    accuracy_summary = accuracy_summary,
    rep_raw = rep_raw,
    V.diag = V.diag,
    Threshold_type = thd.type,
    time_sec = time_sec
  )
}





simu_ret <- function(methods = c("TFMseg", "TFMseg.vec", "FMseg", "LR"),
                     detect_thd = NULL,
                     lrv = TRUE,
                     simu_set_detect = NULL,
                     Time_list = NULL,
                     dim_obs_list = NULL,
                     dim_latent = c(3, 3, 3),
                     data_setting = c("s0", "s1", "s1_var1", "s3_3I", "s3_A2", "s3_A1", "s2"),
                     dep = TRUE,
                     theta_coef = NULL,
                     r_hat = NULL,
                     trim_coef = 1/4,
                     nrep = 100,
                     rvs = FALSE,
                     miss = FALSE) {
  
  data_setting <- match.arg(data_setting)
  
  if (is.null(simu_set_detect)) {
    if (is.null(Time_list) || is.null(dim_obs_list)) {
      stop("Provide either simu_set_detect, or both Time_list and dim_obs_list.")
    }
    
    simu_set_detect <- unlist(
      lapply(dim_obs_list, function(dobs) {
        lapply(Time_list, function(Time) {
          list(
            Time = as.integer(Time),
            dim_obs = as.integer(dobs),
            dim_latent = as.integer(dim_latent)
          )
        })
      }),
      recursive = FALSE
    )
  }
  
  results <- list()
  total_jobs <- length(methods) * length(simu_set_detect)
  job_id <- 0L
  
  for (method in methods) {
    thd.type_iter <- if (method == "TFMseg") {
      "fixed"
    } else if (method == "TFMseg.vec") {
      "oracle"
    } else {
      "--"
    }
    
    V_shap_iter <- switch(
      method,
      "FMseg" = "diag",
      "LR" = "--",
      "diag"
    )
    
    for (thd.type in thd.type_iter) {
      for (V_shap in V_shap_iter) {
        for (setting in simu_set_detect) {
          
          job_id <- job_id + 1L
          
          Time <- setting$Time
          dim_obs <- setting$dim_obs
          dim_latent_i <- setting$dim_latent
          
          theta_local <- if (is.null(theta_coef)) c(0.25, 0.5, 0.75) else theta_coef
          theta_int <- floor(Time * theta_local)
          theta_int <- sort(unique(theta_int[theta_int > 0 & theta_int < Time]))
          m_local <- length(theta_int)
          
          message(sprintf(
            "[simu_detect] %d/%d: Time=%d, dim_obs=%s, setting=%s, method=%s",
            job_id, total_jobs, Time, paste(dim_obs, collapse = "x"),
            data_setting, method
          ))
          
          res <- simu_comparison(
            method = method,
            nrep = nrep,
            data_setting = data_setting,
            dim_obs = dim_obs,
            dim_latent = dim_latent_i,
            Time = Time,
            dep = dep,
            m = m_local,
            theta_coef = theta_local,
            r_hat = r_hat,
            trim_coef = trim_coef,
            detect_thd = detect_thd,
            thd.type = if (method %in% c("TFMseg", "TFMseg.vec")) thd.type else "fixed",
            V.diag = (V_shap == "diag"),
            lrv = lrv,
            rvs = rvs,
            miss = miss
          )
          
          results <- append(results, list(list(
            method = method,
            Time = Time,
            dim_obs = dim_obs,
            freq_summary = res$freq_summary,
            accuracy_summary = res$accuracy_summary,
            cp_est_scaled = res$cp_est_scaled,
            rep_raw = res$rep_raw,
            V_shap = V_shap,
            Threshold_type = if (!is.null(res$Threshold_type)) res$Threshold_type else thd.type,
            time_sec = res$time_sec
          )))
        }
      }
    }
  }
  
  summary_table <- do.call(rbind, lapply(results, function(res) {
    data.frame(
      Time = res$Time,
      dim_obs = paste(res$dim_obs, collapse = "x"),
      Method = res$method,
      Threshold = ifelse(res$method == "LR", "--", res$Threshold_type),
      LRV = ifelse(res$V_shap == "--", "--",
                   ifelse(res$V_shap == "diag", "Diagonal", "Full")),
      `<= -2` = res$freq_summary[1],
      `-1` = res$freq_summary[2],
      `0` = res$freq_summary[3],
      `1` = res$freq_summary[4],
      `>= 2` = res$freq_summary[5],
      Accuracy_j1 = res$accuracy_summary[1] %||% NA_real_,
      Accuracy_j2 = res$accuracy_summary[2] %||% NA_real_,
      Accuracy_j3 = res$accuracy_summary[3] %||% NA_real_,
      Mean_time_s = mean(res$time_sec, na.rm = TRUE),
      SD_time_s = stats::sd(res$time_sec, na.rm = TRUE)
    )
  }))
  
  names(summary_table) <- c(
    "Time", "dim_obs", "Method", "Threshold", "LRV",
    "$\\leq -2$", "$-1$", "$0$", "$1$", "$\\geq 2$",
    "$j = 1$", "$j = 2$", "$j = 3$",
    "Mean time (s)", "SD time (s)"
  )
  
  summary_table$Method <- factor(
    summary_table$Method,
    levels = c("TFMseg", "TFMseg.vec", "FMseg", "LR")
  )
  
  summary_table <- summary_table[order(
    summary_table$Time,
    summary_table$dim_obs,
    summary_table$Method,
    summary_table$Threshold,
    summary_table$LRV
  ), ]
  
  table_output <- knitr::kable(summary_table, "html", escape = FALSE, row.names = FALSE) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    ) %>%
    kableExtra::add_header_above(
      c(" " = 5, "$\\widehat{m} - m$" = 5, "Accuracy" = 3, "Runtime (s)" = 2),
      escape = FALSE
    ) %>%
    kableExtra::add_footnote(
      paste0("Summary of change point estimation over ", nrep, " realisations"),
      notation = "none"
    )
  
  all_cp_est_scaled <- do.call(c, lapply(results, function(res) res$cp_est_scaled))
  
  cp_est_df <- do.call(rbind, lapply(results, function(r) {
    do.call(rbind, lapply(seq_along(r$cp_est_scaled), function(i) {
      cp_vals <- r$cp_est_scaled[[i]]
      if (length(cp_vals) == 0) {
        data.frame(
          scaled_cp = NA_real_,
          sim = i,
          Time = r$Time,
          dim_obs = paste(r$dim_obs, collapse = "x"),
          V_shap = ifelse(r$V_shap == "diag", "Diagonal",
                          ifelse(r$V_shap == "full", "Full", "--")),
          method = r$method,
          threshold = r$Threshold_type
        )
      } else {
        data.frame(
          scaled_cp = cp_vals,
          sim = i,
          Time = r$Time,
          dim_obs = paste(r$dim_obs, collapse = "x"),
          V_shap = ifelse(r$V_shap == "diag", "Diagonal",
                          ifelse(r$V_shap == "full", "Full", "--")),
          method = r$method,
          threshold = r$Threshold_type
        )
      }
    }))
  }))
  
  cp_est_df <- dplyr::mutate(
    cp_est_df,
    Group = dplyr::case_when(
      method == "TFMseg" & threshold == "fixed"  & V_shap == "Diagonal" ~ "TFMseg: fixed, diag",
      method == "TFMseg" & threshold == "fixed"  & V_shap == "Full" ~ "TFMseg: fixed, full",
      method == "TFMseg" & threshold == "oracle" & V_shap == "Diagonal" ~ "TFMseg: oracle, diag",
      method == "TFMseg" & threshold == "oracle" & V_shap == "Full" ~ "TFMseg: oracle, full",
      method == "TFMseg.vec" & threshold == "oracle" & V_shap == "Diagonal" ~ "TFMseg.vec: oracle, diag",
      method == "TFMseg.vec" & threshold == "oracle" & V_shap == "Full" ~ "TFMseg.vec: oracle, full",
      method == "FMseg" ~ "FMseg: oracle m, diag",
      method == "LR" ~ "LR",
      TRUE ~ paste(method, threshold, V_shap, sep = ", ")
    )
  )
  
  df_raw <- do.call(rbind, lapply(results, function(r) {
    out <- r$rep_raw
    out$Time <- r$Time
    out$dim_obs <- paste(r$dim_obs, collapse = "x")
    out$method <- r$method
    out
  }))
  
  df_raw <- df_raw[, c(
    "Time", "dim_obs", "method", "rep_id",
    "le_m2", "m1", "m0", "p1", "ge_p2",
    "acc1", "acc2", "acc3", "time"
  )]
  
  check_counts <- df_raw %>%
    dplyr::count(Time, dim_obs, method, name = "n") %>%
    dplyr::arrange(Time, dim_obs, method)
  
  list(
    summary = summary_table,
    results = results,
    table = table_output,
    cp_est_scaled = all_cp_est_scaled,
    cp_est_df = cp_est_df,
    df_raw = df_raw,
    check_counts = check_counts
  )
}







extract_accurate <- function(res,
                             Time_list = NULL,
                             dim_obs_list = NULL,
                             source_names = NULL,
                             q = 3L) {
  stopifnot(is.list(res), !is.null(res$df_raw), is.data.frame(res$df_raw))
  stopifnot(!is.null(res$cp_est_df), is.data.frame(res$cp_est_df))
  
  if (!is.null(q)) {
    if (!is.numeric(q) || length(q) != 1L || is.na(q) ||
        q < 1 || q != as.integer(q)) {
      stop("`q` must be NULL or a positive integer.")
    }
    q <- as.integer(q)
  }
  
  df_raw <- res$df_raw
  cp_df  <- res$cp_est_df
  
  if (!is.null(Time_list)) {
    df_raw <- df_raw %>% dplyr::filter(Time %in% Time_list)
    cp_df  <- cp_df  %>% dplyr::filter(Time %in% Time_list)
  }
  
  if (!is.null(source_names)) {
    df_raw <- df_raw %>% dplyr::filter(method %in% source_names)
    cp_df  <- cp_df  %>% dplyr::filter(method %in% source_names)
  }
  
  if (!is.null(dim_obs_list)) {
    dim_keep <- vapply(dim_obs_list, paste, collapse = "x", FUN.VALUE = character(1))
    df_raw <- df_raw %>% dplyr::filter(dim_obs %in% dim_keep)
    cp_df  <- cp_df  %>% dplyr::filter(dim_obs %in% dim_keep)
  }
  
  df_good <- df_raw
  
  if (!is.null(q)) {
    acc_cols <- paste0("acc", seq_len(q))
    miss_acc <- setdiff(acc_cols, names(df_good))
    
    if (length(miss_acc) > 0L) {
      stop(
        sprintf(
          "For q = %d, the following columns are missing in res$df_raw: %s",
          q, paste(miss_acc, collapse = ", ")
        )
      )
    }
    
    df_good <- df_good %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(acc_cols), ~ . == 1))
  }
  
  keys <- df_good %>%
    dplyr::transmute(Time, dim_obs, method, sim = rep_id) %>%
    dplyr::distinct()
  
  cp_good <- cp_df %>%
    dplyr::semi_join(keys, by = c("Time", "dim_obs", "method", "sim")) %>%
    dplyr::arrange(Time, dim_obs, method, sim, scaled_cp)
  
  out <- res
  out$df_raw <- df_good
  out$cp_est_df <- cp_good
  
  out
}
