single_identify_detection <- function(M,
                            dim_obs = c(100, 100),
                            dim_latent = c(3, 3),
                            TT = 500,
                            true_cp = TT / 2,
                            true_mode = 1,
                            reps = 100,
                            seed = 1,
                            thd_D = 0.10,
                            thd_sd = 1,
                            detection = FALSE) {
  
  K <- length(dim_obs)
  stopifnot(true_mode %in% seq_len(K))
  
  id_mode <- integer(reps)
  id_method <- character(reps)
  Dmat <- matrix(NA_real_, reps, K)   # projector distances
  SDmat <- matrix(NA_real_, reps, K)   # ratio SDs
  cp_hat <- numeric(reps) 
  
  for (rep in seq_len(reps)) {
    set.seed(seed + rep - 1)
    
    tr_list <- vector("list", K)
    tr_list[[true_mode]] <- list(M, NULL)
    
    dat <- dgp_general(
      model = ifelse(K == 1, "vector", ifelse(K == 2, "matrix", "tensor")),
      dim_obs = dim_obs,
      dim_latent = dim_latent,
      Time = TT,
      coeff = 0.7,
      true_cp = true_cp,
      type_mode = list(c(true_mode)),
      type_change = list("l"),
      transform_list = tr_list)
    
    X <- dat$X
    Time <- dim(X)[3]
    dim_X <- dim(X)[1:2]
    res_load <- global_pca(X, dim_X)
    G <- res_load$G
    G_dim <- as.vector(res_load$r_hat)
    est_cp <- cand_sbs(G, G_dim, single = TRUE)$est.cp
    cp_use <- if (!detection) true_cp else est_cp
    cp_hat[rep] <- cp_use 
    
    out <- mode_identify(dat$X, dim_obs,
                         cp = cp_use,
                         st = 1,
                         ed = TT,
                         thd_D = thd_D,
                         thd_sd = thd_sd)
    
    id_mode[rep] <- if (length(out$change_mode) == 0) NA else out$change_mode
    id_method[rep] <- out$method
    Dmat[rep, ] <- out$D
    SDmat[rep, ] <- out$ratio_sd
  }
  
  Dmat <- round(Dmat,  4)
  SDmat <- round(SDmat, 4)
  
  colnames(Dmat) <- paste0("D_mode",  seq_len(K))
  colnames(SDmat) <- paste0("SD_mode", seq_len(K))
  
  results <- data.frame(
    rep = seq_len(reps),
    cp = cp_hat,
    identified_mode = id_mode,
    method = id_method,
    Dmat,
    SDmat
  )
  
  accuracy <- round(mean(id_mode == true_mode, na.rm = TRUE), 4)
  
  list(results = results, accuracy = accuracy)
}
