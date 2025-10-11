TNotSBS_simu <- function(
    nrep = 50,
    data_setting = c("s1","s2"),
    Time_list = c(400, 600, 800, 1000),
    dim_obs_list = list(c(20, 30, 50), c(10, 20, 30)),
    dim_latent = c(2, 2, 2),
    theta_coef = c(0.25, 0.5, 0.75),
    threshold_coef = NULL,
    dep = TRUE,
    method = c("fixed", "oracle")
) {
  # Parallel package
  #if (!requireNamespace("future.apply", quietly = TRUE)) {install.packages("future.apply")}
  library(future.apply)
  plan(multisession)
  
  method <- match.arg(method)
  results <- list()
  
  for (Time in Time_list) {
    true_cp <- c(round(Time * theta_coef[1]), round(Time * theta_coef[2]), round(Time * theta_coef[3]))
    for (dim_obs in dim_obs_list) {
      detected_cp_list <- future_lapply(seq_len(nrep), function(rep) {
        set.seed(rep+888)
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
          true_cp <- true_cp
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
            list(C1, diag(3), diag(3)),
            list(diag(3), C3, diag(3)), 
            list(diag(3), diag(3), C2)
          ) #231#132
          type_mode <- list(1, 2, 3)
          type_change <- list("l", "l", "l")
          shift_ind <- NULL
          shift_mean <- NULL
          shift_var <- NULL
          add_factors <- NULL
          add_factors_coeff <- NULL
          true_cp <- true_cp
        }
        
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
          add_factors_coeff = add_factors_coeff)
        
        X <- data_sim$X
        dim_X <- dim(X)[1:3]
        est_load <- global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
        #est_load <- global_pca(X = X, dim_X = dim_X, r_hat = c(3,3,3), proj = TRUE)
        G <- est_load$G_proj
        G_dim <- as.vector(est_load$r_hat)
        #print(G_dim)
        if (method == "oracle") {
          out <- TNotSBS(G, G_dim, method = c("oracle"), m = 3, V.diag = TRUE, lrv = TRUE)
        } else if (method == "fixed"){
          dr <- sum(G_dim * (G_dim + 1) / 2)
          int_len <- round(6*log(Time))
          thd <- pmax(exp(threshold_coef[1] * log(log(Time/int_len)) + threshold_coef[2] * log(dr)),
                      threshold_coef[3] * log(Time))
          out <- TNotSBS(G, G_dim, method = c("fixed"), threshold = thd, V.diag = TRUE, lrv = TRUE)
        }
        sort(out$est.cp)
      }, future.seed = TRUE)
      
      key <- paste("Time", Time, "dim", paste(dim_obs, collapse = "_"), sep = "_")
      
      grDevices::dev.new(noRStudioGD = TRUE)
      visualise_stats(detected_cp_list, true_cp, Time, method_name = paste0("TNotSBS-oracle (", key, ")"))
      fig <- grDevices::recordPlot()
      grDevices::dev.off()
      
      stats <- compute_stats(Time, true_cp, detected_cp_list, method_name = paste0("TNotSBS-oracle (", key, ")"))
      
      results[[key]] <- list(
        Time = Time,
        dim_obs = dim_obs,
        true_cp = true_cp,
        detected_cp_list = detected_cp_list,
        stats = stats,
        figure = fig
      )
      cat(sprintf("Finished: %s\n", key))
    }
  }
  return(results)
}
