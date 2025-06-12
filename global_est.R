# Global PCA (can "global" on a given window)
global_pca <- function(X, dim_X,
                       st = NULL,
                       ed = NULL,
                       centre  = FALSE) {
  
  if (centre) X <- centre_X(X)
  
  K <- length(dim_X)
  all_dims <- dim(X)
  TT <- all_dims[K + 1]       
  
  if (is.null(st)) st <- 1
  if (is.null(ed)) ed <- TT
  if (st < 1 || ed > TT || st > ed)
    stop("`st` and `ed` must satisfy 1 ≤ st ≤ ed ≤ Time.")
  
  idx_time <- st:ed           
  T_win <- length(idx_time)      
  
  r_hat <- integer(K)
  Lambda <- vector("list", K)
  
  X_win <- slice_time(X, idx_time)
  
  for (i in seq_len(K)) {
    perm <- c(i, setdiff(seq_len(K), i), K + 1)
    Xp <- aperm(X_win, perm)
    Xmat <- matrix(Xp, nrow = dim_X[i])
    
    M_i <- Xmat %*% t(Xmat) / T_win
    r_i <- median(abc.factor.number(M_i)$r[4:6])
    eig <- eigen(M_i, symmetric = TRUE)
    
    Lambda[[i]] <- eig$vectors[, seq_len(r_i), drop = FALSE] * sqrt(dim_X[i])
    r_hat[i] <- r_i
  }
  names(r_hat) <- sprintf("r_mode%d", seq_len(K))
  names(Lambda) <- sprintf("Lambda_mode%d", seq_len(K))
  
  ## pseudo factors
  if (K == 1) {
    G <- t(Lambda[[1]]) %*% X_win
  } else {
    G_dims <- c(r_hat, T_win)
    G <- array(0, G_dims)
    
    for (t_idx in seq_len(T_win)) {
      Xt <- slice_time(X_win, t_idx)
      gt <- Xt
      for (i in seq_len(K))
        gt <- mode_n_mult(gt, t(Lambda[[i]]), i)
      
      if (K == 2) G[, , t_idx] <- gt
      else if (K == 3) G[ , , , t_idx] <- gt
      else stop("Only K = 1, 2, 3 are supported.")
    }
  }

  
  cc <- est_common_component(G, Lambda, T_win)
  
  list(
    r_hat  = r_hat,
    Lambda = Lambda,
    G = G,     
    cc = cc 
  )
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