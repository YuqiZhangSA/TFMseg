#' Global PCA estimation to derive time-invariant loadings and pseudo factors
#'
#' Performs a one‐shot principal component analysis on an array \code{X} of arbitrary order (vector, matrix, or tensor)
#' by:
#'   1. Calculating mode-wise covariance matrices across time;
#'   2. Selecting the number of factors in each mode via the ABC median rule (Haeran MOSUM paper);
#'   3. Computing orthonormal loading matrices scaled by \(\sqrt{p_i}\);
#'   4. Extracting pseudo factors \(G_t\) by successive mode-\(n\) projections of each time-slice.
#'
#' @param X An array of size \code{c(p1, p2, …, pK, T)}, where modes \(1\ldots K\) index the spatial dimensions and
#'   mode \(K+1\) indexes time (\(T\) observations).
#' @param dim_X Integer vector of length \(K\), giving the observation dimension \code{p_i} for each mode \(i=1,\dots,K\).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{r_hat}}{Integer vector of length \(K\) with the estimated number of factors in each mode.}
#'     \item{\code{Lambda}}{List of length \(K\); each element is a \code{p_i × r_hat[i]} matrix of estimated loadings, scaled by \(\sqrt{p_i}\).}
#'     \item{\code{G}}{An array of size \code{c(r_hat, T)} giving the pseudo factors \(G_t\) for each time \(t\).}
#'   }
#'
#' @details
#' For each mode \(i\):
#' * We permute \code{X} so that mode \(i\) is the first dimension and time is the last, then reshape to a matrix
#'   \(\mathbf{X}_{(i)}\) of size \(p_i × (T \prod_{j\neq i} p_j)\).
#' * We form the mode-wise covariance \(M_i = \frac1T\,\mathbf{X}_{(i)}\,\mathbf{X}_{(i)}^\top\).
#' * The ABC median rule \code{abc.factor.number} is applied to \code{M_i} to choose \(r_i\).
#' * We take the first \(r_i\) eigenvectors of \code{M_i}, scale each by \(\sqrt{p_i}\), and store as \(\Lambda_i\).
#'
#' The pseudo‐factor at time \(t\) is then
#' \[
#'   G_t \;=\;X_t
#'     \times_{i=1}^K \Lambda_i^\top
#' \]
#' i.e.\ successive \(n\)-mode products of each slice \(X_t\) by \(\Lambda_i^\top\).
#'
#'
#' @export


gloal_pca <- function(X, dim_X) {
  #X <- centre_X(X)
  
  K    <- length(dim_X)
  dims <- dim(X)
  Time <- dims[K + 1]
  
  # r_hat and loadings
  r_hat  <- integer(K)
  Lambda <- vector("list", K)
  for (i in seq_len(K)) {
    perm <- c(i, setdiff(seq_len(K), i), K + 1)
    Xp   <- aperm(X, perm)
    Xmat <- matrix(Xp, nrow = dim_X[i])
    M_i  <- Xmat %*% t(Xmat) / Time
    
    r_i <- median(abc.factor.number(M_i)$r[4:6])
    eig <- eigen(M_i, symmetric = TRUE)
    Vi <- eig$vectors[, seq_len(r_i), drop = FALSE]
    
    Lambda[[i]] <- Vi * sqrt(dim_X[i])
    r_hat[i] <- r_i
  }
  names(r_hat) <- paste0("r_mode", seq_len(K))
  names(Lambda) <- paste0("Lambda_mode", seq_len(K))
  
  # K==1, vector
  if (K == 1) {
    G <- t(Lambda[[1]]) %*% X
    
  } else {
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
    
    G_dims <- c(r_hat, Time)
    G <- array(0, G_dims)
    
    for (t in seq_len(Time)) {
      # extract the t-th slice
      idx <- c(rep(list(TRUE), K), list(t))
      X_t <- do.call(`[`, c(list(X), idx))
      
      g_t <- X_t
      for (i in seq_len(K)) {
        g_t <- mode_n_mult(g_t, t(Lambda[[i]]), i)
      }
      
      if (K == 2) {
        G[,, t] <- g_t
      } else if (K == 3) {
        G[,,, t] <- g_t
      } else {
        stop("Only K = 1,2,3 are supported.")
      }
    }
  }
  
  list(
    r_hat = r_hat,    
    Lambda = Lambda, 
    G = G     
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

