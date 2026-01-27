library(lubridate)
library(dplyr)
library(kableExtra)

`%||%` <- function(x, y) if (is.null(x)) y else x

# -------------------------------------------------------------------
# 0
# -------------------------------------------------------------------
# strong in the sense of \Vert HH^\trans - I \Vert_F, check with the following
# frob_dev <- function(H) norm(H %*% t(H) - diag(nrow(H)), type = "F")
# frob_dev(make_rot_mat_mode1(3, "weak")); frob_dev(make_rot_mat_mode1(3, "strong"))
make_rot_mat_mode1 <- function(r, strength = c("weak", "strong")) {
  strength <- match.arg(strength)
  if (r < 2) stop("Need dim_latent[1] >= 2 for a rotation/transform.")
  
  if (r == 2) {
    return(if (strength == "weak")
      matrix(c(1, 0.15,
               0, 1), 2, 2, byrow = TRUE)
      else
        matrix(c(1, 0.80,
                 0, 1), 2, 2, byrow = TRUE)
    )
  }
  
  H3_weak <- matrix(c(
    1, 0.15, 0,
    0, 1, 0.15,
    0, 0, 1
  ), 3, 3, byrow = TRUE)
  
  H3_strong <- matrix(c(
    1, 0.80, 0.75,
    0, 1, 0.80,
    0, 0, 1
  ), 3, 3, byrow = TRUE)
  
  H <- diag(r)
  H[1:3, 1:3] <- if (strength == "weak") H3_weak else H3_strong
  H
}

# -------------------------------------------------------------------
# 1
# -------------------------------------------------------------------
simu_comparison <- function(method = c("TFMseg", "TFMseg.vec", "FMseg", "LR"),
                            nrep = 100,
                            seed_start = 900,
                            data_setting = c("s1","s2","s2-2","s3","s4","s5","s6"),
                            dim_obs,
                            dim_latent,
                            Time,
                            dist = c("Gaussian", "heavy"),
                            coeff = 0.7,
                            dep = TRUE,
                            m = NULL,                
                            theta_coef = NULL,
                            r_hat,
                            trim_coef = 3,
                            threshold_coef = NULL,
                            thd.type = c("fixed", "oracle"),
                            V.diag = TRUE,
                            lrv = TRUE,
                            acc_coef = 2) {
  
  method <- match.arg(method)
  thd.type <- match.arg(thd.type)
  data_setting <- match.arg(data_setting)
  dist <- match.arg(dist)
  
  need_fns <- c("dgp_general", "global_pca", "TFMseg", "bs_LR_globaltrim")
  missing_fns <- need_fns[!vapply(need_fns, exists, logical(1), mode = "function")]
  if (length(missing_fns)) {
    stop("Missing required function(s): ", paste(missing_fns, collapse = ", "))
  }
  
  if (is.null(theta_coef)) theta_coef <- c(0.25, 0.5, 0.75)
  theta <- floor(Time * theta_coef)
  theta <- sort(unique(theta[theta > 0 & theta < Time]))
  
  # s3–s6 for single CP
  if (data_setting %in% c("s3","s4","s5","s6") && length(theta) != 1L) {
    stop("For s3–s6 please set theta_coef to a single value, e.g. theta_coef = 0.5.")
  }
  
  trim <- round(trim_coef * (Time / log(Time)))
  
  frequency_table <- matrix(0L, nrow = nrep, ncol = 5L)
  accuracy_count <- if (length(theta)) matrix(0L, nrep, length(theta)) else matrix(NA_real_, nrep, 0L)
  detected_cp_list <- vector("list", nrep)
  cp_est_scaled <- vector("list", nrep)
  time_sec <- numeric(nrep)
  
  for (sim in seq_len(nrep)) {
    set.seed(seed_start + sim)
    
    # -----------------------------------------------------------------
    # Settings
    # -----------------------------------------------------------------
    if (data_setting == "s1") {
      type_mode   <- list(c(1), c(2), c(1,3))
      type_change <- list("l", "f", c("f", "l"))
      shift_ind <- list(list(c(dim_obs[1]/2, dim_latent[1]), NULL, NULL),
                        list(NULL, NULL, NULL),
                        list(NULL, NULL, c(dim_obs[3]/2, dim_latent[3]/2)))
      shift_mean <- c(1, 0, 0)
      shift_var  <- c(2^2, 0, 1^2)
      transform_list <- NULL
      add_factors <- list(c(0, 0, 0),  c(0, 3, 0), c(1, 0, 0))
      add_factors_coeff <- list(c(0, 0, 0),  c(0, 0.6, 0), c(0.3, 0, 0))
      true_cp <- theta
      
    } else if (data_setting == "s2-2") {
      r0 <- 3
      C1 <- matrix(rnorm(r0^2, sd = 1), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
      # increase variance
      C2 <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 0
      ), nrow = 3, byrow = TRUE)
      
      C3 <- matrix(0, 3, 3)
      C3[1,1] <- 0.5; C3[2,1] <- rnorm(1, 0, 1); C3[2,2] <- 1
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
      
    }else if (data_setting == "s2") {
      r0 <- 3
      C1 <- matrix(rnorm(r0^2, sd = 1/sqrt(r0)), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
      
      C2 <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 0
      ), nrow = 3, byrow = TRUE)
      
      C3 <- matrix(0, 3, 3)
      C3[1,1] <- 0.5; C3[2,1] <- rnorm(1, 0, 1); C3[2,2] <- 1
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
      
    } else if (data_setting %in% c("s3","s4")) {
      
      type_mode <- list(1)     # only mode 1 changes
      type_change <- list("l") 
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
      R1 <- make_rot_mat_mode1(dim_latent[1],
                               strength = if (data_setting == "s3") "weak" else "strong")
      
      transform_list <- list(list(R1, diag(dim_latent[2]), diag(dim_latent[3])))
      true_cp <- theta
      
    } else if (data_setting %in% c("s5","s6")) {
      stopifnot(length(theta) == 1L)
      
      stopifnot(length(dim_latent) == 2L)
      #if (any(dim_latent != c(5, 5))) warning("s5/s6 designed for dim_latent = c(5,5).")
      
      type_mode <- list(1)
      type_change <- list("l")
      shift_ind <- NULL
      shift_mean <- NULL
      shift_var <- NULL
      add_factors <- NULL
      add_factors_coeff <- NULL
      
      diag_mat <- function(r, a) {
        vals <- pmax(1 - a * (0:(r - 1)), 0.1)
        diag(vals)
      }
      
      R_weak <- diag_mat(dim_latent[1], 0.1)
      R_strong <- diag_mat(dim_latent[1], 0.2)
      # check magnitude change
      #frob_dev(diag_mat(2, 0.1))
      #frob_dev(diag_mat(2, 0.2))
      R1 <- if (data_setting == "s5") R_weak else R_strong
      
      transform_list <- list(list(R1, diag(dim_latent[2])))
      true_cp <- theta
    }
    
    m_true <- length(true_cp)
    
    model_used <- if (data_setting %in% c("s5","s6")) "matrix" else "tensor"
    
    data_sim <- dgp_general(
      model = model_used,
      dim_obs = dim_obs,
      dim_latent = dim_latent,
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
    K <- length(dim(X)) - 1L
    dim_X <- dim(X)[seq_len(K)]
    
    # -----------------------------------------------------------------
    # Different methods
    # -----------------------------------------------------------------
    if (method == "TFMseg") {
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
        out <- TFMseg(G, G_dim,
                       trim = trim,
                       method = "oracle",
                       m = m_true,
                       V.diag = V.diag,
                       lrv = lrv)
      } else {
        dr <- sum(G_dim * (G_dim + 1L) / 2L)
        if (is.null(threshold_coef) || length(threshold_coef) < 3)
          stop("For thd.type = 'fixed', provide threshold_coef of length 3.")
        thd <- pmax(
          exp(threshold_coef[1] * log(log(Time)) + threshold_coef[2] * log(dr)),
          threshold_coef[3] * log(Time)
        )
        out <- TFMseg(G, G_dim,
                       method = "fixed",
                       threshold = thd,
                       V.diag = V.diag,
                       lrv = lrv)
      }
      
      detected_cp <- out$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
      
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
      
      out <- TFMseg(G, G_dim,
                     trim = trim,
                     method = "oracle",
                     m = m_true,
                     V.diag = V.diag,
                     lrv = lrv) #check if works for vector
      
      detected_cp <- out$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
      
    } else if (method == "FMseg") {
      start_time <- Sys.time()
      
      if (is.null(r_hat)) {
        est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
        G_dim <- as.vector(est_load$r_hat)
      } else {
        est_load <- global_pca(X = X, dim_X = dim_X, r_hat = r_hat, proj = TRUE)
        G_dim <- as.vector(r_hat)
      }
      
      Xd <- matrix(aperm(X, c(K + 1L, seq_len(K))), nrow = Time)  # Time by p_vec
      X_vec <- t(Xd)                                               # p_vec by Time
      
      r_vec <- prod(G_dim)
      lbd <- floor(Time/log(Time))
      
      # IMPORTANT: oracle m must be m_true
      detected_cp <- FMSeg(x = X_vec, r = r_vec, m = m_true, trim = trim,
                           method = "oracle", lbd = lbd)$est.cp %||% integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
      
    } else {  # LR
      start_time <- Sys.time()
      
      Xd <- matrix(aperm(X, c(K + 1L, seq_len(K))), nrow = Time)
      
      tau_global <- if (length(theta_coef) >= 2) min(diff(theta_coef)) else 0.2
      min_size <- round(tau_global * Time)
      
      lr_out <- bs_LR_globaltrim(
        Xd,
        tau_global = tau_global,
        min_size = min_size,
        r_est = NULL,
        tau_for_cv = 0.1
      )
      
      detected_cp <- (lr_out$cps %||% integer(0))
      detected_cp <- if (length(detected_cp)) sort(unique(detected_cp)) else integer(0)
      
      end_time <- Sys.time()
      time_sec[sim] <- as.numeric(as.period(end_time - start_time, unit = "sec"))
    }
    
    # -----------------------------------------------------------------
    # Summaries: cp num dist + accuracy
    # -----------------------------------------------------------------
    detected_cp_sorted <- if (length(detected_cp)) sort(detected_cp) else integer(0)
    detected_cp_list[[sim]] <- detected_cp_sorted
    cp_est_scaled[[sim]] <- detected_cp_sorted / Time
    
    # m-hat - m (m must be the truth here)
    m_est_diff <- length(detected_cp_sorted) - m_true
    frequency_table[sim, 1] <- as.integer(m_est_diff <= -2)
    frequency_table[sim, 2] <- as.integer(m_est_diff == -1)
    frequency_table[sim, 3] <- as.integer(m_est_diff == 0)
    frequency_table[sim, 4] <- as.integer(m_est_diff == 1)
    frequency_table[sim, 5] <- as.integer(m_est_diff >= 2)
    
    # Accuracy per true CP: within acc_win = acc_coef * log(Time)
    if (length(theta)) {
      acc_win <- max(1L, as.integer(round(acc_coef * log(Time))))
      for (j in seq_along(theta)) {
        accuracy_count[sim, j] <- as.integer(
          any(abs(detected_cp_sorted - theta[j]) <= acc_win)
        )
      }
    }
  }
  
  freq_summary <- colSums(frequency_table) / nrep
  accuracy_summary <- if (ncol(accuracy_count)) colMeans(accuracy_count) else numeric(0)
  
  list(
    cp_est_scaled = cp_est_scaled,
    freq_summary = freq_summary,
    accuracy_summary = accuracy_summary,
    V.diag = V.diag,
    Threshold_type = thd.type,
    time_sec = time_sec
  )
}

# -------------------------------------------------------------------
# 2
# -------------------------------------------------------------------
simu_ret <- function(methods = c("TFMseg", "TFMseg.vec", "FMseg", "LR"),
                     thd.type_values = c("fixed", "oracle"),
                     V_shap_values = c("diag", "full"),
                     lrv = TRUE,
                     simu_set_detect,
                     dist = "Gaussian",
                     data_setting = c("s1","s2","s2-2","s3","s4","s5","s6"),
                     dep = TRUE,
                     theta_coef = NULL,
                     r_hat = NULL,
                     trim_coef = 3,
                     nrep,
                     threshold_coef = NULL,
                     acc_coef = 2) {
  
  data_setting <- match.arg(data_setting)
  
  # Guard: s3–s6 require single CP
  if (data_setting %in% c("s3","s4","s5","s6")) {
    if (is.null(theta_coef) || length(theta_coef) != 1L) {
      stop("For s3–s6 you must set theta_coef to a single value (e.g. 0.5).")
    }
  }
  
  results <- list()
  
  for (method in methods) {
    thd.type_iter <- if (method == "TFMseg") {
      thd.type_values
    } else if (method == "TFMseg.vec") {
      "oracle"
    } else {
      "--"
    }
    
    V_shap_iter <- switch(method,
                          "FMseg" = "diag",
                          "LR"    = "--",
                          V_shap_values
    )
    
    for (thd.type in thd.type_iter) {
      for (V_shap in V_shap_iter) {
        for (setting in simu_set_detect) {
          
          Time <- setting$Time
          dim_obs <- setting$dim_obs
          dim_latent <- setting$dim_latent
          
          theta_local <- if (is.null(theta_coef)) c(0.25, 0.5, 0.75) else theta_coef
          theta_int <- floor(Time * theta_local)
          theta_int <- sort(unique(theta_int[theta_int > 0 & theta_int < Time]))
          m_local <- length(theta_int)
          
          cat(sprintf(
            "\nRunning simulation: Method = %s, T = %d, dim_obs = %s, V_shap = %s, thd.type = %s, m_true = %d\n",
            method, Time, paste(dim_obs, collapse = "x"), V_shap, thd.type, m_local
          ))
          
          res <- simu_comparison(
            method = method,
            nrep = nrep,
            seed_start = 900,
            data_setting = data_setting,
            dim_obs = dim_obs,
            dim_latent = dim_latent,
            Time = Time,
            dist = dist,
            coeff = 0.7,
            dep = dep,
            m = m_local,      
            theta_coef = theta_local,
            r_hat = r_hat,
            trim_coef = trim_coef,
            threshold_coef = threshold_coef,
            thd.type = if (method %in% c("TFMseg","TFMseg.vec")) thd.type else "fixed",
            V.diag = (V_shap == "diag"),
            lrv = lrv,
            acc_coef = acc_coef
          )
          
          results <- append(results, list(list(
            method = method,
            Time = Time,
            dim_obs = dim_obs,
            freq_summary = res$freq_summary,
            accuracy_summary = res$accuracy_summary,
            cp_est_scaled = res$cp_est_scaled,
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
  
  summary_table$Method <- factor(summary_table$Method,
                                 levels = c("TFMseg", "TFMseg.vec", "FMseg", "LR"))
  
  summary_table <- summary_table[order(summary_table$Time,
                                       summary_table$dim_obs,
                                       summary_table$Method,
                                       summary_table$Threshold,
                                       summary_table$LRV), ]
  
  table_output <- kable(summary_table, "html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
    add_header_above(c(" " = 5, "$\\widehat{m} - m$" = 5, "Accuracy" = 3, "Runtime (s)" = 2), escape = FALSE) %>%
    add_footnote(paste0("Summary of change point estimation over ", nrep, " realisations"),
                 notation = "none")
  
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
  
  cp_est_df <- cp_est_df %>%
    mutate(Group = case_when(
      method == "TFMseg" & threshold == "fixed"  & V_shap == "Diagonal" ~ "TFMseg: fixed, diag",
      method == "TFMseg" & threshold == "fixed"  & V_shap == "Full" ~ "TFMseg: fixed, full",
      method == "TFMseg" & threshold == "oracle" & V_shap == "Diagonal" ~ "TFMseg: oracle, diag",
      method == "TFMseg" & threshold == "oracle" & V_shap == "Full" ~ "TFMseg: oracle, full",
      method == "TFMseg.vec" & threshold == "oracle" & V_shap == "Diagonal" ~ "TFMseg.vec: oracle, diag",
      method == "TFMseg.vec" & threshold == "oracle" & V_shap == "Full" ~ "TFMseg.vec: oracle, full",
      method == "FMseg" ~ "FMseg: oracle m, diag",
      method == "LR" ~ "LR",
      TRUE ~ paste(method, threshold, V_shap, sep = ", ")
    ))
  
  list(
    summary = summary_table,
    results = results,
    table = table_output,
    cp_est_scaled = all_cp_est_scaled,
    cp_est_df = cp_est_df
  )
}
