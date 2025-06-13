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

