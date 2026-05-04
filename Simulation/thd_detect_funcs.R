info_record <- function(num_simu,
                        seed_start = 888,
                        Time,
                        dim_obs = c(10,10,10),
                        dim_latent  = c(3,3,3),
                        r_hat = NULL,
                        dist = "Gaussian",
                        setting = c("s0", "s1", "s2", "s2_222", "s2_232", "s2_233", "s7", "s8_C2", "s8_C3"),
                        theta_coef = NULL,
                        trim = NULL,
                        threshold = NULL,
                        lbd = NULL,
                        method = c("fixed", "oracle"),
                        V.diag = TRUE,
                        lrv = TRUE,
                        dep = TRUE) {
  
  method  <- match.arg(method)
  setting <- match.arg(setting)
  
  if (is.null(trim)) trim <- floor(0.25 * Time / log(Time))
  if (is.null(lbd))  lbd  <- floor(0.5  * Time / log(Time))
  
  # Only define theta/m when not null
  if (setting == "s0") {
    theta <- integer(0)
    m <- 0
  } else {
    if (is.null(theta_coef)) theta_coef <- c(0.25, 0.5, 0.75)
    theta <- c(floor(Time * theta_coef[1]),
               floor(Time * theta_coef[2]),
               floor(Time * theta_coef[3]))
    m <- length(theta)
  }
  
  # -----------------------------
  # DGP settings
  # -----------------------------
  if (setting == "s1") {
    type_mode   <- list(c(2,3), c(1), c(3,1))
    type_change <- list(c("f","l"), "f", c("l","l"))
    
    make_3I <- function(r) 3 * diag(r)
    make_C2 <- function(r) { M <- diag(r); M[r,r] <- 0; M }
    make_Orth <- function(r, angle = pi/5) {
      if (r <= 1) return(diag(r))
      R2 <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
      M  <- diag(r); M[1:2, 1:2] <- R2; M
    }
    
    transform_list <- list(
      list(NULL, NULL, make_3I),
      NULL,
      list(make_Orth, NULL, make_C2)
    )
    
    add_factors <- list(c(0,1,0), c(2,0,0), c(0,0,0))
    add_factors_coeff <- list(c(0,0.3,0), c(0.3,0,0), c(0,0,0))
    shift_ind  <- replicate(3, list(NULL, NULL, NULL), simplify = FALSE)
    shift_mean <- c(0,0,0)
    shift_var  <- c(0,0,0)
    true_cp <- theta
    
  } else if (setting == "s2") {
    r0 <- 3
    C1 <- matrix(rnorm(r0^2, sd = 1/sqrt(r0)), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
    C2 <- matrix(c(1,0,0, 0,1,0, 0,0,0), 3, 3, byrow = TRUE)
    C3 <- matrix(0, 3, 3)
    C3[1,1] <- 0.5; C3[2,1] <- rnorm(1); C3[2,2] <- 1
    C3[3,1] <- rnorm(1); C3[3,2] <- rnorm(1); C3[3,3] <- 1.5
    
    transform_list <- list(
      list(C3,  diag(3), diag(3)),
      list(diag(3), C2,  diag(3)),
      list(diag(3), diag(3), C1)
    )
    type_mode   <- list(1,2,3)
    type_change <- list("l","l","l")
    shift_ind <- shift_mean <- shift_var <- NULL
    add_factors <- add_factors_coeff <- NULL
    true_cp <- theta
    
  } else if (setting == "s0") {
    true_cp <- integer(0)
    type_mode <- list()
    type_change <- list()
    transform_list <- NULL
    shift_ind <- NULL
    shift_mean <- NULL
    shift_var <- NULL
    add_factors <- NULL
    add_factors_coeff <- NULL
    
  } else if (setting == "s2_222") {
    r0 <- 2
    C1 <- matrix(rnorm(r0^2, sd = 1/sqrt(r0)), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
    C2 <- matrix(c(1,0, 0,0), 2, 2, byrow = TRUE)
    C3 <- matrix(0, 2, 2)
    C3[1,1] <- 0.5; C3[2,1] <- rnorm(1); C3[2,2] <- 1.5
    
    transform_list <- list(
      list(C3,  diag(2), diag(2)),
      list(diag(2), C2,  diag(2)),
      list(diag(2), diag(2), C1)
    )
    type_mode <- list(1,2,3)
    type_change <- list("l","l","l")
    shift_ind <- shift_mean <- shift_var <- NULL
    add_factors <- add_factors_coeff <- NULL
    true_cp <- theta
    
  } else if (setting == "s2_232") {
    r0 <- 2
    C1 <- matrix(rnorm(r0^2, sd = 1/sqrt(r0)), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
    C2 <- matrix(c(1,0,0, 0,1,0, 0,0,0), 3, 3, byrow = TRUE)
    C3 <- matrix(0, 2, 2)
    C3[1,1] <- 0.5; C3[2,1] <- rnorm(1); C3[2,2] <- 1.5
    
    transform_list <- list(
      list(C3,  diag(2), diag(2)),
      list(diag(3), C2,  diag(3)),
      list(diag(2), diag(2), C1)
    )
    type_mode   <- list(1,2,3)
    type_change <- list("l","l","l")
    shift_ind <- shift_mean <- shift_var <- NULL
    add_factors <- add_factors_coeff <- NULL
    true_cp <- theta
    
  } else if (setting == "s2_233") {
    r0 <- 3
    C1 <- matrix(rnorm(r0^2, sd = 1/sqrt(r0)), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
    C2 <- matrix(c(1,0,0, 0,1,0, 0,0,0), 3, 3, byrow = TRUE)
    C3 <- matrix(0, 2, 2)
    C3[1,1] <- 0.5; C3[2,1] <- rnorm(1); C3[2,2] <- 1.5
    
    transform_list <- list(
      list(C3,  diag(2), diag(2)),
      list(diag(3), C2,  diag(3)),
      list(diag(3), diag(3), C1)
    )
    type_mode   <- list(1,2,3)
    type_change <- list("l","l","l")
    shift_ind <- shift_mean <- shift_var <- NULL
    add_factors <- add_factors_coeff <- NULL
    true_cp <- theta
  }
  
  detected_cp_list <- vector("list", num_simu)
  cp_info <- vector("list", num_simu)
  dr_values <- numeric(num_simu)
  
  for (sim in seq_len(num_simu)) {
    set.seed(sim + seed_start)
    
    data_sim <- dgp_general(
      model = "tensor",
      dim_obs = dim_obs,
      dim_latent  = dim_latent,
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
      dist = dist
    )
    
    X_t <- data_sim$X
    dim_X <- dim(X_t)[1:3]
    
    est_load <- global_pca(X = X_t, dim_X = dim_X, r_hat = r_hat, centre = TRUE, proj = TRUE)
    G <- est_load$G_proj
    G_dim <- as.vector(est_load$r_hat)
    
    dr_values[sim] <- sum(G_dim * (G_dim + 1) / 2)
    
    if (setting == "s0") {
      out <- TFMseg(G = G, G_dim = G_dim,
                    method = "fixed",
                    threshold = Inf,
                    V.diag = V.diag, lrv = lrv)
      
      est_cp <- integer(0)
      aligned_deviations <- character(0)
      selected_thd <- NA_real_
      next_highest <- NA_real_
      
    } else {
      out <- TFMseg(G = G, G_dim = G_dim,
                    method = "oracle",
                    m = m,
                    V.diag = V.diag, lrv = lrv)
      
      est_cp <- sort(out$est.cp)
      aligned_deviations <- align_change_points(est_cp, theta, Time, trim)
      selected_thd <- out$selected_threshold[1]
      next_highest <- out$next_highest_threshold[1]
    }
    
    detected_cp_list[[sim]] <- est_cp
    
    cp_info[[sim]] <- data.frame(
      seed = seed_start + sim,
      True_CP = if (length(theta) == 0) "" else paste(theta, collapse = ", "),
      Detected_CP = if (length(est_cp) == 0) "" else paste(est_cp, collapse = ", "),
      Est_Error = if (length(aligned_deviations) == 0) "" else paste(aligned_deviations, collapse = ", "),
      Trim = trim,
      Min_Interval_Length = lbd,
      dr = dr_values[sim],
      Selected_Threshold = selected_thd,
      Next_Highest_CUSUM = next_highest,
      
      NullMax = out$null_max[1],
      NullMax_excl = out$null_max_excl[1],
      
      NullMax_len = out$null_max_len[1],
      NullMax_bound = out$null_max_bound[1],
      NullMax_st = out$null_max_st[1],
      NullMax_ed = out$null_max_ed[1],
      NullMax_cp = out$null_max_cp[1],
      
      Next_len = if (!is.null(out$next_len)) out$next_len[1] else NA_integer_,
      Next_bound = if (!is.null(out$next_bound)) out$next_bound[1] else NA_real_,
      Next_st = if (!is.null(out$next_st)) out$next_st[1] else NA_integer_,
      Next_ed = if (!is.null(out$next_ed)) out$next_ed[1] else NA_integer_,
      Next_cp = if (!is.null(out$next_cp)) out$next_cp[1] else NA_integer_
    )
  }
  
  cp_info_df <- do.call(rbind, cp_info)
  
  return(list(
    detected_cp_list = detected_cp_list,
    cp_info = cp_info_df,
    theta = theta
  ))
}





quantile_record <- function(num_simu,
                            seed_start = 888,
                            Time,
                            dim_obs = c(10,10,10),
                            dim_latent  = c(3,3,3),
                            r_hat = NULL,
                            dist = "Gaussian",
                            setting = c("s0", "s1", "s2", "s2_222", "s2_232", "s2_233", "s7", "s8_C2", "s8_C3"),
                            theta_coef,
                            trim = NULL,
                            threshold = NULL,
                            lbd = NULL,
                            method = c("fixed", "oracle"),
                            V.diag = TRUE,
                            lrv = TRUE,
                            dep = TRUE,
                            quantile_per = c(0.9, 0.95, 0.99)) {
  
  simu_results <- info_record(num_simu,
                              seed_start,
                              Time,
                              dim_obs,
                              dim_latent,
                              r_hat,
                              dist,
                              setting,
                              theta_coef,
                              trim,
                              threshold,
                              lbd,
                              method,
                              V.diag,
                              lrv,
                              dep)
  
  cp_info_df <- simu_results$cp_info
  
  # quantiles for null maxima
  quant_NullMax <- quantile(cp_info_df$NullMax, probs = quantile_per, na.rm = TRUE)
  quant_NullMax_excl <- quantile(cp_info_df$NullMax_excl, probs = quantile_per, na.rm = TRUE)
  
  quant_Next_Highest_CUSUM <- quantile(cp_info_df$Next_Highest_CUSUM, probs = quantile_per, na.rm = TRUE)
  quant_Selected_Threshold <- quantile(cp_info_df$Selected_Threshold, probs = quantile_per, na.rm = TRUE)
  
  # lower quantiles
  quant_Next_Highest_CUSUM_re <- quantile(cp_info_df$Next_Highest_CUSUM, probs = 1 - quantile_per, na.rm = TRUE)
  quant_Selected_Threshold_re <- quantile(cp_info_df$Selected_Threshold, probs = 1 - quantile_per, na.rm = TRUE)
  quant_NullMax_re <- quantile(cp_info_df$NullMax, probs = 1 - quantile_per, na.rm = TRUE)
  quant_NullMax_excl_re <- quantile(cp_info_df$NullMax_excl, probs = 1 - quantile_per, na.rm = TRUE)
  
  # null-max interval diagnostics
  quant_NullMax_len   <- quantile(cp_info_df$NullMax_len,   probs = quantile_per, na.rm = TRUE)
  quant_NullMax_bound <- quantile(cp_info_df$NullMax_bound, probs = quantile_per, na.rm = TRUE)
  
  quant_interval_length <- quantile(cp_info_df$Min_Interval_Length, probs = quantile_per, na.rm = TRUE)
  quant_dr <- quantile(cp_info_df$dr, probs = quantile_per, na.rm = TRUE)
  
  return(list(
    quant_Next_Highest_CUSUM = quant_Next_Highest_CUSUM,
    quant_Selected_Threshold = quant_Selected_Threshold,
    
    quant_NullMax = quant_NullMax,
    quant_NullMax_excl = quant_NullMax_excl,
    
    quant_Next_Highest_CUSUM_re = quant_Next_Highest_CUSUM_re,
    quant_Selected_Threshold_re = quant_Selected_Threshold_re,
    quant_NullMax_re = quant_NullMax_re,
    quant_NullMax_excl_re = quant_NullMax_excl_re,
    
    quant_interval_length = quant_interval_length,
    quant_dr = quant_dr,
    quant_NullMax_len = quant_NullMax_len,
    quant_NullMax_bound = quant_NullMax_bound,
    
    cp_info = cp_info_df
  ))
}




quantile_set <- function(
    simulation_settings,
    dist = "Gaussian",
    method = c("fixed", "oracle"),
    V_shap_values = c("diag", "full"),
    setting = c("s0", "s1", "s2", "s2_222", "s2_232", "s2_233", "s7", "s8_C2", "s8_C3"),
    num_simu,
    r_hat = NULL,
    lrv = TRUE,
    dep = TRUE,
    theta_coef = NULL,
    quantile_per = c(0.9, 0.95, 0.99)
) {
  results <- list()
  
  for (thd.type in method) {
    for (V_shap in V_shap_values) {
      for (simu_set in simulation_settings) {
        
        Time <- simu_set$Time
        dim_obs <- simu_set$dim_obs
        dim_latent <- simu_set$dim_latent
        
        cat(sprintf(
          "\nRunning simulation: T = %d, dim_obs = (%s), dim_latent = (%s), V_shap = %s, thd.type = %s\n",
          Time, paste(dim_obs, collapse = ","), paste(dim_latent, collapse = ","), V_shap, thd.type
        ))
        
        quantiles <- quantile_record(
          num_simu,
          seed_start = 888,
          Time = Time,
          dim_obs = dim_obs,
          dim_latent = dim_latent,
          r_hat = r_hat,
          dist = dist,
          setting = setting,
          theta_coef = theta_coef,
          method = thd.type,
          V.diag = (V_shap == "diag"),
          lrv = lrv,
          dep = dep,
          quantile_per = quantile_per
        )
        
        results <- append(results, list(list(
          Time = Time,
          dim_obs = dim_obs,
          dim_latent = dim_latent,
          
          quant_Next_Highest_CUSUM = quantiles$quant_Next_Highest_CUSUM,
          quant_Selected_Threshold = quantiles$quant_Selected_Threshold,
          
          # null maxima
          quant_NullMax = quantiles$quant_NullMax,
          quant_NullMax_excl = quantiles$quant_NullMax_excl,
          
          quant_Next_Highest_CUSUM_re = quantiles$quant_Next_Highest_CUSUM_re,
          quant_Selected_Threshold_re = quantiles$quant_Selected_Threshold_re,
          quant_NullMax_re = quantiles$quant_NullMax_re,
          quant_NullMax_excl_re = quantiles$quant_NullMax_excl_re,
          
          quant_interval_length = quantiles$quant_interval_length,
          quant_dr = quantiles$quant_dr,
          quant_NullMax_len = quantiles$quant_NullMax_len,
          quant_NullMax_bound = quantiles$quant_NullMax_bound,
          
          V_shap = V_shap,
          Threshold_type = thd.type
        )))
      }
    }
  }
  
  results
}
