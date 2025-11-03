library(lubridate)
library(dplyr)
library(kableExtra)

`%||%` <- function(x, y) if (is.null(x)) y else x

simu_comparison <- function(method = c("TNotSBS", "FMseg", "LR"), 
                            nrep = 100, 
                            seed_start = 888, 
                            data_setting = c("s1","s2"),
                            dim_obs,
                            dim_latent,
                            Time,
                            dist = c("Gaussian", "heavy"),
                            coeff = 0.7,
                            dep = TRUE,
                            m = NULL, 
                            theta_coef = NULL, 
                            r_hat,
                            trim = round(2*log(Time)), 
                            threshold_coef = NULL, 
                            thd.type = c("fixed", "oracle"), 
                            V.diag = TRUE, 
                            lrv = TRUE,
                            acc_coef = 2) {
  
  method   <- match.arg(method)
  thd.type <- match.arg(thd.type)
  data_setting <- match.arg(data_setting)
  
  if (is.null(theta_coef)) theta_coef <- c(0.25, 0.5, 0.75)
  
  if (!all(theta_coef == c(0, 0, 0))) {
    theta <- floor(Time * theta_coef)
  } else {
    theta <- integer(0)
  }
  if (is.null(m)) m <- length(theta)
  
  frequency_table <- matrix(0L, nrow = nrep, ncol = 5L)
  accuracy_count  <- if (length(theta)) matrix(0L, nrep, length(theta)) else matrix(NA_real_, nrep, 0L)
  detected_cp_list <- vector("list", nrep)
  cp_est_scaled <- vector("list", nrep)
  time_sec <- numeric(nrep) 
  
  for (sim in seq_len(nrep)) {
    set.seed(seed_start + sim)
    
    if (data_setting == "s1") {
      type_mode <- list(c(1), c(2), c(1,3))
      type_change <- list("l", "f", c("f", "l"))
      shift_ind <- list(list(c(dim_obs[1]/2,dim_latent[1]), NULL, NULL),
                        list(NULL, NULL, NULL), 
                        list(NULL, NULL, c(dim_obs[3]/2,dim_latent[3]/2)))
      shift_mean <- c(1, 0, 0)
      shift_var <- c(2^2, 0, 1^2)
      transform_list <- NULL
      add_factors <- list(c(0, 0, 0),  c(0, 3, 0), c(1, 0, 0))
      add_factors_coeff <- list(c(0, 0, 0),  c(0, 0.6, 0), c(0.3, 0, 0))
      true_cp <- theta
    } else if (data_setting == "s2") {
      r0 <- 3
      C1 <- matrix(rnorm(r0^2, mean = 0, sd = 1/sqrt(r0)), nrow = r0)
      C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
      C2 <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 0
      ), nrow = 3, byrow = TRUE)
      C3 <- matrix(0, 3, 3)
      C3[1,1] <- 0.5; C3[2,1] <- rnorm(1, 0, 1); C3[2,2] <- 1; 
      C3[3,1] <- rnorm(1, 0, 1); C3[3,2] <- rnorm(1, 0, 1); C3[3,3] <- 1.5
      transform_list <- list(
        list(C3, diag(3), diag(3)),
        list(diag(3), C2, diag(3)), 
        list(diag(3), diag(3), C1)
      )
      type_mode <- list(1, 2, 3)
      type_change <- list("l", "l", "l")
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      true_cp <- theta
    }
    
    data_sim <- dgp_general(
      model = "tensor",
      dim_obs = dim_obs,
      dim_latent  = dim_latent,
      Time = Time,
      dep = dep,
      coeff = coeff,
      true_cp = true_cp,
      type_mode = type_mode,
      type_change = type_change,
      shift_ind = shift_ind,
      shift_mean = shift_mean,
      shift_var = shift_var,
      transform_list = transform_list,
      add_factors = add_factors,
      add_factors_coeff = add_factors_coeff,
      dist = dist
    )
    
    X <- data_sim$X
    dim_X <- dim(X)[1:3]
    
    if (method == "TNotSBS") {
      start_time <- Sys.time()
      
      if (is.null(r_hat)) {
        est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
        G <- est_load$G_proj
        G_dim <- as.vector(est_load$r_hat)
      } else {
        est_load <- global_pca(X = X, dim_X = dim_X, r_hat = r_hat, proj = TRUE)
        G <- est_load$G_proj
        G_dim <- as.vector(r_hat)
      }
      
      if (thd.type == "oracle") {
        out <- TNotSBS(G, G_dim, method = "oracle", m = 3, V.diag = V.diag, lrv = lrv)
      } else {  # fixed
        dr <- sum(G_dim * (G_dim + 1L) / 2L)
        int_len <- max(2L, round(6 * log(Time)))
        if (is.null(threshold_coef) || length(threshold_coef) < 3)
          stop("For thd.type = 'fixed', provide threshold_coef of length 3.")
        thd <- pmax(
          exp(threshold_coef[1] * log(log(Time/int_len)) + threshold_coef[2] * log(dr)),
          threshold_coef[3] * log(Time)
        )
        out <- TNotSBS(G, G_dim, method = "fixed", threshold = thd, V.diag = V.diag, lrv = lrv)
      }
      detected_cp <- out$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
      
    } else if (method == "FMseg") {
      start_time <- Sys.time()
      
      if (is.null(r_hat)) {
        est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
        G <- est_load$G_proj
        G_dim <- as.vector(est_load$r_hat)
      } else {
        est_load <- global_pca(X = X, dim_X = dim_X, r_hat = r_hat, proj = TRUE)
        G <- est_load$G_proj
        G_dim <- as.vector(r_hat)
      }
      
      X_vec <- t(matrix(aperm(X, c(4, 1, 2, 3)), nrow = Time))   # p-by-T
      r_vec <- prod(G_dim)
      p_vec <- nrow(X_vec)
      lbd <- round(Time^(max(2/5, 1 - min(1, log(p_vec)/log(Time)))) * log(Time)^1.1)
      detected_cp <- FMSeg(x = X_vec, r = r_vec, m = 3, trim = trim, method = "oracle", lbd = lbd)$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
      
    } else {  # LR
      start_time <- Sys.time()
      
      Xd <- matrix(aperm(X, c(4, 1, 2, 3)), nrow = Time)
      #detected_cp <- recurse_LR(Xd, tau1 = 0.2, offset = 0)$cps %||% integer(0)
      #detected_cp <- recurse_LR(Xd, tau1 = min(diff(theta_coef)), offset = 0)$cps %||% integer(0)
      #detected_cp <- sort(detected_cp)
      lr_out <- bs_LR_globaltrim(Xd, tau_global = min(diff(theta_coef)), min_size= round(min(diff(theta_coef))*Time), #min_size= trim,
                                 r_est = NULL, tau_for_cv = 0.1)
      detected_cp <- (lr_out$cps %||% integer(0))
      detected_cp <- if (length(detected_cp)) sort(unique(detected_cp)) else integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
    }
    
    detected_cp_sorted <- if (length(detected_cp)) sort(detected_cp) else integer(0)
    detected_cp_list[[sim]] <- detected_cp_sorted
    cp_est_scaled[[sim]] <- detected_cp_sorted / Time
    
    m_est_diff <- length(detected_cp_sorted) - m
    frequency_table[sim, 1] <- as.integer(m_est_diff <= -2)
    frequency_table[sim, 2] <- as.integer(m_est_diff == -1)
    frequency_table[sim, 3] <- as.integer(m_est_diff == 0)
    frequency_table[sim, 4] <- as.integer(m_est_diff == 1)
    frequency_table[sim, 5] <- as.integer(m_est_diff >= 2)
    
    if (length(theta)) {
      for (j in seq_along(theta)) {
        accuracy_count[sim, j] <- as.integer(any(abs(detected_cp_sorted - theta[j]) <= acc_coef * round(log(Time))))
      }
    }
  }
  
  freq_summary <- colSums(frequency_table) / nrep
  accuracy_summary <- if (ncol(accuracy_count)) colMeans(accuracy_count) else numeric(0)
  
  list(
    cp_est_scaled = cp_est_scaled,
    freq_summary  = freq_summary,
    accuracy_summary = accuracy_summary,
    V.diag = V.diag,
    Threshold_type = thd.type,
    time_sec = time_sec 
  )
}


simu_ret <- function(methods = c("TNotSBS", "FMseg", "LR"), 
                     thd.type_values = c("fixed", "oracle"), 
                     V_shap_values = c("diag", "full"), lrv = TRUE,
                     simu_set_detect, 
                     dist = "Gaussian",
                     data_setting = c("s1","s2"),
                     dep = TRUE,
                     theta_coef = NULL, 
                     r_hat = NULL,
                     trim = round(2*log(Time)), 
                     m = 3, 
                     nrep, 
                     threshold_coef = NULL,
                     acc_coef = 2) {
  
  results <- list()
  data_setting <- match.arg(data_setting)
  
  for (method in methods) {
    thd.type_iter <- if (method == "TNotSBS") thd.type_values else "--"
    V_shap_iter <- switch(method, "FMseg" = "diag", "LR" = "--", V_shap_values)
    
    for (thd.type in thd.type_iter) {
      for (V_shap in V_shap_iter) {
        for (setting in simu_set_detect) {
          Time <- setting$Time
          dim_obs <- setting$dim_obs
          dim_latent <- setting$dim_latent
          
          cat(sprintf("\nRunning simulation: Method = %s, T = %d, dim_obs = %s, V_shap = %s, thd.type = %s\n", 
                      method, Time, paste(dim_obs, collapse = "x"), V_shap, thd.type))
          
          res <- simu_comparison(method = method, 
                                 nrep = nrep, 
                                 data_setting = data_setting,
                                 dim_obs = dim_obs,
                                 dim_latent = dim_latent,
                                 Time = Time,
                                 dist = dist,
                                 dep = dep,
                                 m = m, 
                                 theta_coef = theta_coef, 
                                 r_hat = r_hat,
                                 trim = round(2*log(Time)), 
                                 threshold_coef = threshold_coef, 
                                 thd.type = if (method == "TNotSBS") thd.type else "fixed", 
                                 V.diag = (V_shap == "diag"), 
                                 lrv = lrv,
                                 acc_coef = acc_coef)
          
          results <- append(results, list(list(
            method           = method,
            Time             = Time,
            dim_obs          = dim_obs,
            freq_summary     = res$freq_summary,
            accuracy_summary = res$accuracy_summary,
            cp_est_scaled    = res$cp_est_scaled,
            V_shap           = V_shap,
            Threshold_type   = if (!is.null(res$Threshold_type)) res$Threshold_type else thd.type,
            time_sec         = res$time_sec 
          )))
        }
      }
    }
  }
  
  summary_table <- do.call(rbind, lapply(results, function(res) {
    data.frame(
      Time        = res$Time,
      dim_obs     = paste(res$dim_obs, collapse = "x"),
      Method      = res$method,
      Threshold   = ifelse(res$method == "LR", "--", res$Threshold_type),
      LRV         = ifelse(res$V_shap == "--", "--",
                           ifelse(res$V_shap == "diag", "Diagonal", "Full")),
      `<= -2`     = res$freq_summary[1],
      `-1`        = res$freq_summary[2],
      `0`         = res$freq_summary[3],
      `1`         = res$freq_summary[4],
      `>= 2`      = res$freq_summary[5],
      Accuracy_j1 = res$accuracy_summary[1] %||% NA_real_,
      Accuracy_j2 = res$accuracy_summary[2] %||% NA_real_,
      Accuracy_j3 = res$accuracy_summary[3] %||% NA_real_,
      Mean_time_s = mean(res$time_sec, na.rm = TRUE),
      SD_time_s   = stats::sd(res$time_sec, na.rm = TRUE)
    )
  }))
  
  names(summary_table) <- c(
    "Time", "dim_obs", "Method", "Threshold", "LRV",
    "$\\leq -2$", "$-1$", "$0$", "$1$", "$\\geq 2$",
    "$j = 1$", "$j = 2$", "$j = 3$",
    "Mean time (s)", "SD time (s)"
  )
  
  summary_table$Method <- factor(summary_table$Method, levels = c("TNotSBS", "FMseg", "LR"))
  
  summary_table <- summary_table[order(summary_table$Time,
                                       summary_table$dim_obs,
                                       summary_table$Method,
                                       summary_table$Threshold,
                                       summary_table$LRV), ]
  
  table_output <- kable(summary_table, "html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
    add_header_above(c(" " = 5, "$\\widehat{m} - m$" = 5, "Accuracy" = 3, "Runtime (s)" = 2), escape = FALSE) %>%
    add_footnote(paste0("Summary of change point estimation by TNotSBS (and/or FMseg, LR) over ", 
                        nrep, " realisations"), notation = "none")
  
  all_cp_est_scaled <- do.call(c, lapply(results, function(res) res$cp_est_scaled))
  cp_est_df <- do.call(rbind, lapply(results, function(r) {
    do.call(rbind, lapply(seq_along(r$cp_est_scaled), function(i) {
      cp_vals <- r$cp_est_scaled[[i]]
      if (length(cp_vals) == 0) {
        data.frame(scaled_cp = NA_real_, sim = i, Time = r$Time, dim_obs = paste(r$dim_obs, collapse = "x"),
                   V_shap = ifelse(r$V_shap == "diag", "Diagonal", ifelse(r$V_shap == "full", "Full", "--")),
                   method = r$method, threshold = r$Threshold_type)
      } else {
        data.frame(scaled_cp = cp_vals, sim = i, Time = r$Time, dim_obs = paste(r$dim_obs, collapse = "x"),
                   V_shap = ifelse(r$V_shap == "diag", "Diagonal", ifelse(r$V_shap == "full", "Full", "--")),
                   method = r$method, threshold = r$Threshold_type)
      }
    }))
  }))
  
  cp_est_df <- cp_est_df %>%
    mutate(Group = case_when(
      method == "TNotSBS" & threshold == "fixed"  & V_shap == "Diagonal" ~ "TNotSBS: fixed, diag",
      method == "TNotSBS" & threshold == "fixed"  & V_shap == "Full"     ~ "TNotSBS: fixed, full",
      method == "TNotSBS" & threshold == "oracle" & V_shap == "Diagonal" ~ "TNotSBS: oracle, diag",
      method == "TNotSBS" & threshold == "oracle" & V_shap == "Full"     ~ "TNotSBS: oracle, full",
      method == "FMseg"                                                          ~ "FMseg: fixed, diag",
      method == "LR"                                                             ~ "LR",
      TRUE                                                                       ~ paste(method, threshold, V_shap, sep = ", ")
    ))
  
  list(summary = summary_table, 
       results = results, 
       table = table_output, 
       cp_est_scaled = all_cp_est_scaled,
       cp_est_df = cp_est_df)
}
