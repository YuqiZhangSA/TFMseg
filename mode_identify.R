mode_identify <- function(X, dim_obs,
                          cp,
                          st = 1,
                          ed = dim(X)[length(dim_obs) + 1], # by default, single case, whole interval
                          thd_D = 0.10, 
                          thd_sd = 1) { 
  
  if (cp <= st || cp >= ed)
    stop("cp must satisfy st < cp < ed.")
  
  pca_pre <- global_pca(X, dim_obs, st = st, ed = cp)
  pca_post <- global_pca(X, dim_obs, st = cp + 1, ed = ed)
  
  cc_pre <- pca_pre$cc; cc_post <- pca_post$cc
  G_pre <- pca_pre$G; G_post <- pca_post$G
  r_pre <- pca_pre$r_hat; r_post <- pca_post$r_hat
  
  K <- length(dim_obs)
  D <- numeric(K) # eigenspace distance
  ratio_sd  <- numeric(K) 
  ratio_all <- vector("list", K)
  
  n_pre  <- cp - st + 1L
  n_post <- ed - cp
  
  for (k in seq_len(K)) {
    
    ## (A) eigenspace distance from common components
    C0 <- mode_cov(cc_pre, k, 0, n_pre)
    C1 <- mode_cov(cc_post, k, 0, n_post)
    
    U0 <- eigen(C0, TRUE)$vectors[, seq_len(r_pre[k]),  drop = FALSE]
    U1 <- eigen(C1, TRUE)$vectors[, seq_len(r_post[k]), drop = FALSE]
    D[k] <- norm(U1 %*% t(U1) - U0 %*% t(U0), "2")
    
    ## (B) full eigenvalue-ratio spectrum from pseudo factor
    GG0 <- mode_cov(G_pre,  k, 0, n_pre)
    GG1 <- mode_cov(G_post, k, 0, n_post)
    
    eig_pre <- eigen(GG0, TRUE, only.values = TRUE)$values
    eig_post <- eigen(GG1, TRUE, only.values = TRUE)$values
    eig_ratio <- eig_post / eig_pre # element-wise ratios
    
    ratio_all[[k]] <- eig_ratio
    ratio_sd[k] <- sd(eig_ratio) # standard deviation for ratios
  }
  
  ## decision rule
  if (max(D) > thd_D) {
    changed <- which(D == max(D))
    method  <- "rule: eigenspace"
  } else {
    if (max(ratio_sd) < thd_sd) {
      changed <- integer(0)
      method <- "cannot identify"
    } else {
      changed <- which(ratio_sd == max(ratio_sd))
      method <- "rule: eigenvalue-ratio"
    }
  }
  
  list(
    D = D,
    thd_D = thd_D,
    ratio_all = ratio_all,  
    ratio_sd  = ratio_sd,  
    thd_sd = thd_sd,
    method = method,
    change_mode = changed,
    r_pre = r_pre,
    r_post = r_post
  )
}


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


