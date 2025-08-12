global_pca <- function(X, dim_X, r_hat = NULL,
                       st = NULL, ed = NULL,
                       centre = FALSE, proj = TRUE) {
  K <- length(dim_X)
  all_dims <- dim(X)
  TT <- all_dims[K + 1]
  if (centre) X <- centre_along_axis(X, K+1)
  if (is.null(st)) st <- 1
  if (is.null(ed)) ed <- TT
  if (st < 1 || ed > TT || st > ed) stop("`st` and `ed` must satisfy 1 ≤ st ≤ ed ≤ Time.")
  idx_time <- st:ed
  T_win <- length(idx_time)
  X_win <- slice_time(X, idx_time)
  p <- prod(dim_X)
  #if (is.null(r_hat)) {
    #X_perm <- aperm(X_win, c(length(all_dims), seq_len(length(all_dims) - 1)))
    #r_hat <- TFM_FN(X_perm, method = "PE")$factor.num
    #r_hat[1] <- r_hat[1] + 4 # manually add 3
    #r_hat[2] <- r_hat[2] + 3 
    #r_hat[3] <- r_hat[3] + 2 
   #}
  if (is.null(r_hat)) {
    st_r <- max(1L, ceiling(0.65 * TT)) 
    X_for_r <- slice_time(X, st_r:TT)   
    X_perm <- aperm(X_for_r, c(length(dim(X_for_r)), seq_len(length(dim(X_for_r)) - 1)))
    r_hat <- TFM_FN(X_perm, method = "PE")$factor.num
    r_hat[1] <- r_hat[1]+2 # manually add 2
    #r_hat[2] <- r_hat[2]+2 # manually add 1
    #r_hat <- c(5,4,3)
  }
  
  # Step 1: Initial PCA estimator
  Lambda_init <- vector("list", K)
  for (i in seq_len(K)) {
    perm <- c(i, setdiff(seq_len(K), i), K + 1)
    Xp <- aperm(X_win, perm)
    Xmat <- matrix(Xp, nrow = dim_X[i])
    M_i <- Xmat %*% t(Xmat) / T_win
    r_i <- r_hat[i]
    eig <- eigen(M_i, symmetric = TRUE)
    Lambda_init[[i]] <- eig$vectors[, seq_len(r_i), drop = FALSE] * sqrt(dim_X[i])
  }
  names(r_hat) <- sprintf("r_mode%d", seq_len(K))
  names(Lambda_init) <- sprintf("Lambda_init_mode%d", seq_len(K))
  
  # Initial pseudo-factors G and common component
  if (K == 1) {
    G <- t(Lambda_init[[1]]) %*% X_win / p
  } else {
    G_dims <- c(r_hat, T_win)
    G <- array(0, G_dims)
    for (t_idx in seq_len(T_win)) {
      Xt <- slice_time(X_win, t_idx)
      gt <- Xt
      for (i in seq_len(K)) gt <- mode_n_mult(gt, t(Lambda_init[[i]]), i)
      if (K == 2) G[, , t_idx] <- gt / p
      else if (K == 3) G[, , , t_idx] <- gt / p
      else stop("Only K = 1, 2, 3 are supported.")
    }
  }
  cc <- est_common_component(G, Lambda_init, T_win)
  
  if (!proj) {
    return(list(
      r_hat = r_hat,
      Lambda_init = Lambda_init,
      G = G,
      cc = cc
    ))
  }
  
  # Step 2: Projected estimator
  Lambda_proj <- vector("list", K)
  for (k in seq_len(K)) {
    if (K == 1) {
      # For vector case: projected == initial
      Lambda_proj[[k]] <- Lambda_init[[k]]
      next
    }
    Lambda_mk <- Lambda_init[-k]
    Lambda_kron <- kronecker_list(rev(Lambda_mk)) 
    p_k <- dim_X[k]
    p_mk <- prod(dim_X[-k])
    Y_list <- vector("list", T_win)
    perm <- c(k, setdiff(seq_len(K), k), K + 1)
    Xp <- aperm(X_win, perm)
    for (t in seq_len(T_win)) {
      Xkt <- matrix(slice_time(Xp, t), nrow = p_k)
      Ykt_t <- (1 / p_mk) * Xkt %*% Lambda_kron
      Y_list[[t]] <- Ykt_t
    }
    sum_Y <- matrix(0, nrow = p_k, ncol = p_k)
    for (t in seq_len(T_win)) {
      Ykt_t <- Y_list[[t]]  # (p_k by rmk)
      sum_Y <- sum_Y + Ykt_t %*% t(Ykt_t)
    }
    Gamma_Y <- sum_Y / (T_win * p_k)
    r_k <- r_hat[k]
    eigY <- eigen(Gamma_Y, symmetric = TRUE)
    Lambda_proj[[k]] <- eigY$vectors[, seq_len(r_k), drop = FALSE] * sqrt(p_k)
  }
  names(Lambda_proj) <- sprintf("Lambda_proj_mode%d", seq_len(K))
  
  # Projected pseudo-factors and common component
  if (K == 1) {
    G_proj <- t(Lambda_proj[[1]]) %*% X_win / p
  } else {
    G_dims <- c(r_hat, T_win)
    G_proj <- array(0, G_dims)
    for (t_idx in seq_len(T_win)) {
      Xt <- slice_time(X_win, t_idx)
      gt <- Xt
      for (i in seq_len(K)) gt <- mode_n_mult(gt, t(Lambda_proj[[i]]), i)
      if (K == 2) G_proj[, , t_idx] <- gt / p
      else if (K == 3) G_proj[, , , t_idx] <- gt / p
      else stop("Only K = 1, 2, 3 are supported.")
    }
  }
  cc_proj <- est_common_component(G_proj, Lambda_proj, T_win)
  
  return(list(
    r_hat = r_hat,
    Lambda_init = Lambda_init,
    G_init = G,
    cc_init = cc,
    Lambda_proj = Lambda_proj,
    G_proj = G_proj,
    cc_proj = cc_proj
  ))
}


# Helper: Centre array/matrix/tensor along axis
centre_along_axis <- function(A, axis) {
  d <- dim(A)
  perm <- c(setdiff(seq_along(d), axis), axis)
  Aperm <- aperm(A, perm)
  K <- length(d) - 1L
  mean_array <- apply(Aperm, MARGIN = seq_len(K), FUN = mean)
  Acent <- sweep(Aperm, MARGIN = seq_len(K), STATS = mean_array, FUN = "-")
  aperm(Acent, order(perm))
}

# Helper: mode-n multiplication (M %*% unfold_n(tensor))
mode_n_mult <- function(tensor, M, mode) {
  dims_t  <- dim(tensor)
  perm <- c(mode, setdiff(seq_len(length(dims_t)), mode))
  tmp <- aperm(tensor, perm)
  mat <- matrix(tmp, nrow = dims_t[mode])
  prod <- M %*% mat
  new_dims <- c(nrow(M), dims_t[-mode])
  arr <- array(prod, new_dims)
  inv_perm <- match(seq_len(length(dims_t)), perm)
  aperm(arr, inv_perm)
}

# Helper: Extract a slice along the temporal dimension, by index
slice_time <- function(A, ind) {
  K <- length(dim(A)) - 1
  idx <- c(rep(list(TRUE), K), list(ind))
  do.call(`[`, c(list(A), idx, list(drop = FALSE)))
}

# Helper: Remove (collapse) one mode by aggregation
remove_mode <- function(A, mode, FUN = mean, ...) {
  dims <- seq_along(dim(A))
  keep <- dims[-mode]
  res <- apply(A, MARGIN = keep, FUN = FUN, ...)
  if (length(keep) == 1) {
    dim(res) <- dim(A)[keep]
  }
  return(res)
}

# Helper: Kronecker product of a list of matrices (right-to-left order)
kronecker_list <- function(mat_list) {
  if (length(mat_list) == 0) return(matrix(1, 1, 1))  # Identity if empty (K = 1)
  result <- mat_list[[1]]
  if (length(mat_list) == 1) return(result)
  for (i in 2:length(mat_list)) {
    result <- kronecker(result, mat_list[[i]])
  }
  result
}

# Helper: Estimate common component
est_common_component <- function(G_hat, Lambda_list, Time) {
  K <- length(Lambda_list)
  obs_dims <- vapply(Lambda_list, nrow, integer(1))
  cc <- array(0, dim = c(obs_dims, Time))
  for (t in seq_len(Time)) {
    slice <- array(slice_time(G_hat, t), dim = dim(G_hat)[1:K])
    for (k in seq_len(K)) slice <- mode_n_mult(slice, Lambda_list[[k]], k)
    if (K == 1) cc[, t] <- slice
    else if (K == 2) cc[, , t] <- slice
    else if (K == 3) cc[, , , t] <- slice
    else stop("Only K = 1, 2, 3 are supported.")
  }
  cc
}

