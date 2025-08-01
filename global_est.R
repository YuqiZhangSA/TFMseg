# Global PCA (can "global" on a given window)
global_pca <- function(X, dim_X, r_hat = NULL,
                       st = NULL, ed = NULL,
                       centre = FALSE, proj = TRUE) {
  K <- length(dim_X)
  all_dims <- dim(X)
  TT <- all_dims[K + 1]
  
  if (centre) X <- centre_X(X)
  if (is.null(st)) st <- 1
  if (is.null(ed)) ed <- TT
  if (st < 1 || ed > TT || st > ed) stop("`st` and `ed` must satisfy 1 ≤ st ≤ ed ≤ Time.")
  idx_time <- st:ed
  T_win <- length(idx_time)
  X_win <- slice_time(X, idx_time)
  p <- prod(dim_X)
  
  # Step 1: Initial PCA estimator
  if (is.null(r_hat)) { r_hat <- integer(K) }
  Lambda_init <- vector("list", K)
  for (i in seq_len(K)) {
    perm <- c(i, setdiff(seq_len(K), i), K + 1)
    Xp <- aperm(X_win, perm)
    Xmat <- matrix(Xp, nrow = dim_X[i])
    M_i <- Xmat %*% t(Xmat) / T_win
    if (is.null(r_hat) || is.na(r_hat[i]) || r_hat[i] == 0) {
      # change point number est -- use Haeran's method (ABC)
      r_i <- median(abc.factor.number(M_i)$r[4:6])
      r_hat[i] <- r_i
    } else r_i <- r_hat[i]
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
      for (i in seq_len(K))
        gt <- mode_n_mult(gt, t(Lambda_init[[i]]), i)
      if (K == 2) G[, , t_idx] <- gt / p
      else if (K == 3) G[ , , , t_idx] <- gt / p
      else stop("Only K = 1, 2, 3 are supported.")
    }
  }
  cc <- est_common_component(G, Lambda_init, T_win)
  
  # If not projection, return initial
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
    Lambda_mk <- Lambda_init[-k]
    Lambda_kron <- kronecker_list(rev(Lambda_mk)) # reverse for index order
    p_k <- dim_X[k]
    p_mk <- prod(dim_X[-k])
    Y_list <- vector("list", T_win)
    for (t in seq_len(T_win)) {
      Xkt <- matrix(slice_time(X_win, t), nrow = p_k)
      Ykt_t <- (1 / p_mk) * Xkt %*% Lambda_kron
      Y_list[[t]] <- Ykt_t
    }
    sum_Y <- matrix(0, nrow = p_k, ncol = p_k)
    for (t in seq_len(T_win)) {
      Ykt_t <- Y_list[[t]]  # (p_k x rmk)
      sum_Y <- sum_Y + Ykt_t %*% t(Ykt_t)
    }
    Gamma_Y <- sum_Y / (T_win * p_k)
    r_k <- r_hat[k]
    eigY <- eigen(Gamma_Y, symmetric = TRUE)
    Lambda_proj[[k]] <- eigY$vectors[, seq_len(r_k), drop = FALSE] * sqrt(p_k)
  }
  names(Lambda_proj) <- sprintf("Lambda_proj_mode%d", seq_len(K))
  
  # projected pseudo-factors and common component
  if (K == 1) {
    G_proj <- t(Lambda_proj[[1]]) %*% X_win / p
  } else {
    G_dims <- c(r_hat, T_win)
    G_proj <- array(0, G_dims)
    for (t_idx in seq_len(T_win)) {
      Xt <- slice_time(X_win, t_idx)
      gt <- Xt
      for (i in seq_len(K))
        gt <- mode_n_mult(gt, t(Lambda_proj[[i]]), i)
      if (K == 2) G_proj[, , t_idx] <- gt / p
      else if (K == 3) G_proj[ , , , t_idx] <- gt / p
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



#' Centre an array/matrix/tensor along its last (time) dimension
#'
#' @param X An array of dimension c(p1, p2, …, pK, T).  For the vector case, this is a matrix of size p1×T; for the matrix case p1×p2×T; for a 3‐mode tensor p1×p2×p3×T;
#' @return  An array of the same dimensions, where for each fixed index (i1,…,iK) we have X_c[i1,…,iK, t] = X[i1,…,iK, t] − mean_{s=1..T} X[i1,…,iK, s].
#'
centre_X <- function(X) {
  dims <- dim(X)
  K <- length(dims) - 1  
  mean_array <- apply(X, seq_len(K), mean)
  Xc <- sweep(X, seq_len(K), mean_array, FUN = "-")
  return(Xc)
}


# mode-n multiplication
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


## slices the tensor on the last dimension
slice_time <- function(A, ind) {
  K <- length(dim(A)) - 1
  idx <- c(rep(list(TRUE), K), list(ind))
  do.call(`[`, c(list(A), idx, list(drop = FALSE)))
}


## common component estimation
est_common_component <- function(G_hat, Lambda_list, Time) {
  K <- length(Lambda_list)
  obs_dims <- vapply(Lambda_list, nrow, integer(1))
  cc <- array(0, dim = c(obs_dims, Time))
  
  for (t in seq_len(Time)) {
    slice <- array(slice_time(G_hat, t), dim = dim(G_hat)[1:K])
    for (k in seq_len(K))
      slice <- mode_n_mult(slice, Lambda_list[[k]], k)
    
    if (K == 1) cc[, t] <- slice
    else if (K == 2) cc[, , t] <- slice
    else if (K == 3) cc[, , , t] <- slice
    else stop("Only K = 1, 2, 3 are supported.")
  }
  cc
}


kronecker_list <- function(mat_list) {
  result <- mat_list[[1]]
  if (length(mat_list) == 1) return(result)
  for (i in 2:length(mat_list)) {
    result <- kronecker(result, mat_list[[i]])
  }
  result
}
