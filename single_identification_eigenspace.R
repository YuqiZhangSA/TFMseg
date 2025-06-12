# mode-i covariance
mode_cov <- function(G, mode_i, a, b) {
  dims <- dim(G)
  K <- length(dims) - 1
  Time <- dims[K+1]
  ri <- dims[mode_i]
  C <- matrix(0, ri, ri)
  for (t in (a+1):b) {
    if (K == 1) {
      x <- G[, t]
      C <- C + x %*% t(x)
    } else {
      idx <- c(as.list(rep(TRUE, K)), list(t))
      slice <- do.call(`[`, c(list(G), idx, list(drop = FALSE)))
      slice <- array(slice, dim = dims[1:K])
      perm <- c(mode_i, setdiff(seq_len(K), mode_i))
      Ai <- aperm(slice, perm)
      mat_i <- matrix(Ai, nrow = ri)
      C <- C + mat_i %*% t(mat_i)
    }
  }
  C / (b - a)
}



# ============================================================
# single change point identification - repitition

single_identify <- function(M1,
                            dim_obs = c(100, 100),
                            dim_latent = c(3, 3),
                            TT = 500,
                            true_cp = TT / 2,
                            reps = 100,
                            split = TRUE,   # TRUE = two PCAs, before and after
                            seed = 1) {
  
  K <- length(dim_obs)
  Dmat <- matrix(0, reps, K)
  r_pre_table <- matrix(0, reps, K)
  r_post_table <- matrix(0, reps, K)
  
  for (rep in seq_len(reps)) {
    set.seed(seed + rep - 1)
    
    res <- dgp_general(
      model = ifelse(K == 1, "vector",
                              ifelse(K == 2, "matrix", "tensor")),
      dim_obs = dim_obs,
      dim_latent = dim_latent,
      Time = TT,
      coeff = 0.7,
      true_cp = true_cp,
      type_mode = list(c(1)),
      type_change = list("l"),
      transform_list = list(list(M1, NULL))
    )
    
    if (!split) {
      ## [0, Time]
      pca_full <- global_pca(res$X, dim_obs)
      r_pre <- r_post <- pca_full$r_hat
      cc_full <- pca_full$cc 
      
      for (k in seq_len(K)) {
        r_pre_table[rep, k] <- r_post_table[rep, k] <- r_pre[k]
        
        C0 <- mode_cov(cc_full, k, 0, true_cp)
        C1 <- mode_cov(cc_full, k, true_cp, TT)
        
        U0 <- eigen(C0, symmetric = TRUE)$vectors
        U1 <- eigen(C1, symmetric = TRUE)$vectors
        
        Dmat[rep, k] <- norm(U1 %*% t(U1) - U0 %*% t(U0), "2")
      }
      
    } else {
      ## [1,true_cp] & (true_cp, TT]
      pca_pre <- global_pca(res$X, dim_obs, st = 1, ed = true_cp)
      pca_post <- global_pca(res$X, dim_obs, st = true_cp + 1, ed = TT)
      
      r_pre <- pca_pre$r_hat
      r_post <- pca_post$r_hat
      cc_pre <- pca_pre$cc 
      cc_post <- pca_post$cc
      
      for (k in seq_len(K)) {
        r_pre_table[rep,  k] <- r_pre[k]
        r_post_table[rep, k] <- r_post[k]
        
        C0 <- mode_cov(cc_pre, k, 0, true_cp)
        C1 <- mode_cov(cc_post, k, 0, TT - true_cp)
        
        U0 <- eigen(C0, symmetric = TRUE)$vectors[, seq_len(r_pre[k]), drop = FALSE]
        U1 <- eigen(C1, symmetric = TRUE)$vectors[, seq_len(r_post[k]), drop = FALSE]
        
        Dmat[rep, k] <- norm(U1 %*% t(U1) - U0 %*% t(U0), "2")
      }
    }  
  }    

  library(tidyr);  library(dplyr);  library(ggplot2)
  
  dist_long <- as.data.frame(Dmat) |>
    setNames(paste0("Mode", seq_len(K))) |>
    mutate(rep = seq_len(reps)) |>
    pivot_longer(starts_with("Mode"),
                 names_to = "Mode", values_to = "Distance")
  
  rank_long <- data.frame(
    rep = rep(seq_len(reps), each = K),
    Mode = rep(paste0("Mode", seq_len(K)), reps),
    r_pre = as.vector(r_pre_table),
    r_post = as.vector(r_post_table)
  )
  
  list(
    hist_plot = ggplot(dist_long, aes(Distance, fill = Mode)) +
      geom_histogram(alpha = .6, bins = 30, position = "identity") +
      theme_minimal(base_size = 12) +
      labs(title = "Histogram of sub-space distances", y = "count"),
    
    scatter_plot = ggplot(dist_long, aes(rep, Distance, colour = Mode)) +
      geom_point(alpha = .55, size = 2) +
      geom_smooth(se = FALSE, span = .5) +
      theme_minimal(base_size = 12) +
      labs(x = "replication", y = "sub-space distance"),
    
    rank_table = rank_long,
    distance_mat = Dmat,
    count_maxMode = sum(apply(Dmat, 1, \(v) which.max(v) == 1))
  )
}