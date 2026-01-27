info_record <- function(num_simu, 
                        seed_start = 888, 
                        Time, 
                        dim_obs = c(10,10,10),
                        dim_latent  = c(3,3,3),
                        dist = "Gaussian",
                        setting = c("s0", "s1", "s2"), # s0 for no change
                        theta_coef = NULL, 
                        trim = NULL, 
                        threshold = NULL, 
                        lbd = NULL, 
                        method = c("fixed", "oracle"), 
                        V.diag = TRUE, 
                        lrv = TRUE,
                        dep = TRUE) {
  
  method <- match.arg(method)
  setting <- match.arg(setting)
  
  if (is.null(theta_coef)) {theta_coef <- c(0.25,0.5,0.75)}
  theta <- c(floor(Time * theta_coef[1]), floor(Time * theta_coef[2]), floor(Time * theta_coef[3]))
  m <- length(theta_coef)
  if (is.null(trim)) trim <- floor(0.25 * Time / log(Time))
  if (is.null(lbd)) lbd <-  floor(0.5 * Time / log(Time))
  
  # Define settings
  if (setting == "s1") {
    type_mode <- list(
      c(2, 3),       # cp2: modes 2 & 3
      c(1),          # cp3: mode 1
      c(3, 1)        # cp4: modes 3 & 1
    )
    type_change <- list(
      c("f","l"),          # cp2: factor-number + loading change
      "f",                 # cp3: factor-number change
      c("l","l")           # cp4: loading changes
    )
    
    make_I  <- function(r) diag(r)
    make_3I <- function(r) 3 * diag(r)
    make_C2 <- function(r) { M <- diag(r); M[r,r] <- 0; M }
    make_Orth <- function(r, angle = pi/5) {
      if (r <= 1) return(diag(r))
      R2 <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
      M  <- diag(r); M[1:2, 1:2] <- R2; M
    }
    make_C3 <- function(r) {
      vals <- pmax(1 - 0.2 * (0:(r - 1)), 0.1)
      diag(vals)
    }
    
    transform_list <- list(
      list(NULL,     NULL,   make_3I),# cp2: 3I on mode 3
      NULL,                           # cp3: f only
      list(make_Orth,NULL,   make_C2) # cp4: orth on mode 1, C2 on mode 3
    )
    add_factors <- list(
      c(0,1,0),  # cp2: +1 on mode 2
      c(2,0,0),  # cp3: +2 on mode 1
      c(0,0,0)   # cp4 unchanged
    )
    add_factors_coeff <- list(
      c(0,  0.3, 0), # cp2: AR 0.3 for new mode-2 factors
      c(0.3,0,  0),  # cp3: AR 0.3 for new mode-1 factors
      c(0,  0,  0)   # cp4 unchanged
    )
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
      dist = dist)
    
      X_t <- data_sim$X
      dim_X <- dim(X_t)[1:3]
      
      #est_load <- global_pca(X = X_t, dim_X = dim_X, centre = T, proj = TRUE)
      est_load <- global_pca(X = X_t, dim_X = dim_X, r_hat = c(3,3,3), centre = T, proj = TRUE)
      G <- est_load$G_proj
      G_dim <- as.vector(est_load$r_hat)
      
      out <- TFMseg(G = G, G_dim = G_dim, method = c("oracle"), m = m, V.diag = V.diag, lrv = lrv)
      est_cp <- sort(out$est.cp)
      detected_cp_list[[sim]] <- est_cp
    
      dr_values[sim] <- sum(G_dim * (G_dim + 1) / 2)
    
      aligned_deviations <- align_change_points(est_cp, theta, Time, trim)
    
    cp_info[[sim]] <- data.frame(
      seed = seed_start + sim,
      True_CP = paste(theta, collapse = ", "),
      Detected_CP = if (length(est_cp) == 0) "" else paste(est_cp, collapse = ", "),
      Est_Error = paste(aligned_deviations, collapse = ", "),
      Trim = trim,
      Min_Interval_Length = lbd,
      dr = dr_values[sim], 
      Selected_Threshold = out$selected_threshold[1],
      Next_Highest_CUSUM = out$next_highest_threshold[1],
      Max_CUSUM = out$no_change_threshold[1]
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
                            dist = "Gaussian",
                            setting = c("s0", "s1", "s2"), # s0 for no change
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
  
  # Upper quantiles
  quant_Next_Highest_CUSUM <- quantile(cp_info_df$Next_Highest_CUSUM, probs = quantile_per, na.rm = TRUE)
  quant_Selected_Threshold <- quantile(cp_info_df$Selected_Threshold, probs = quantile_per, na.rm = TRUE)
  quant_Max_CUSUM <- quantile(cp_info_df$Max_CUSUM, probs = quantile_per, na.rm = TRUE)
  quant_interval_length <- quantile(cp_info_df$Min_Interval_Length, probs = quantile_per, na.rm = TRUE)
  quant_dr <- quantile(cp_info_df$dr, probs = quantile_per, na.rm = TRUE)
  
  # Lower quantiles 
  quant_Next_Highest_CUSUM_re <- quantile(cp_info_df$Next_Highest_CUSUM, probs = 1 - quantile_per, na.rm = TRUE)
  quant_Selected_Threshold_re <- quantile(cp_info_df$Selected_Threshold, probs = 1 - quantile_per, na.rm = TRUE)
  quant_Max_CUSUM_re <- quantile(cp_info_df$Max_CUSUM, probs = 1 - quantile_per, na.rm = TRUE)
  
  return(list(
    quant_Next_Highest_CUSUM = quant_Next_Highest_CUSUM,
    quant_Selected_Threshold = quant_Selected_Threshold,
    quant_Max_CUSUM = quant_Max_CUSUM,
    quant_Next_Highest_CUSUM_re = quant_Next_Highest_CUSUM_re,
    quant_Selected_Threshold_re = quant_Selected_Threshold_re,
    quant_Max_CUSUM_re = quant_Max_CUSUM_re,
    quant_interval_length = quant_interval_length,
    quant_dr = quant_dr,
    cp_info = cp_info_df
  ))
}


quantile_set <- function(
    simulation_settings, 
    dist = "Gaussian",
    method = c("fixed", "oracle"), 
    V_shap_values = c("diag", "full"), 
    setting = c("s0", "s1", "s2"),
    num_simu, 
    lrv = TRUE,
    dep = TRUE,
    theta_coef = NULL, 
    quantile_per = c(0.9, 0.95, 0.99)
) {
  library(dplyr)
  library(knitr)
  library(kableExtra)

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
        
        quantiles <- quantile_record(num_simu, 
                                     seed_start = 888, 
                                     Time = Time, 
                                     dim_obs = dim_obs,
                                     dim_latent = dim_latent,
                                     dist = dist,
                                     setting = setting, # s0 for no change
                                     theta_coef = theta_coef, 
                                     method = thd.type, 
                                     V.diag = (V_shap == "diag"), 
                                     lrv = lrv,
                                     dep = dep,
                                     quantile_per = c(0.9, 0.95, 0.99)) 
      
        results <- append(results, list(list(
          Time = Time,
          dim_obs = dim_obs,
          dim_latent = dim_latent,
          quant_Next_Highest_CUSUM = quantiles$quant_Next_Highest_CUSUM,
          quant_Selected_Threshold = quantiles$quant_Selected_Threshold,
          quant_Max_CUSUM = quantiles$quant_Max_CUSUM,
          quant_Next_Highest_CUSUM_re = quantiles$quant_Next_Highest_CUSUM_re,
          quant_Selected_Threshold_re = quantiles$quant_Selected_Threshold_re,
          quant_Max_CUSUM_re = quantiles$quant_Max_CUSUM_re,
          quant_interval_length = quantiles$quant_interval_length,  
          quant_dr = quantiles$quant_dr, 
          V_shap = V_shap,
          Threshold_type = thd.type
        )))
      }
    }
  }
  
  return(results)
}

