## lower-triangular vectorisation
vech <- function(M) M[lower.tri(M, diag = TRUE)]

## single change point search on a standardised data matrix D
##
## D : d × T numeric; each column a d-vector
## st : left interval end point 
## ed : right interval end point
## trim: integer; first/last ‘trim’ positions are ignored
##
## value: data.frame(est.cp, val) – 1-row
find_single_cp <- function(D, st, ed, trim){
  Time <- ncol(D)
  st_i <- max(1L, as.integer(st) + 1L) 
  ed_i <- min(Time, as.integer(ed))
  len <- ed_i - st_i + 1L
  stopifnot(len > 2 * trim + 1)
  
  stat  <- numeric(len)
  total <- rowSums(D[, st_i:ed_i, drop = FALSE])
  lsum  <- numeric(nrow(D))
  
  for (k in seq_len(len - 1L)) {
    lsum  <- lsum + D[, st_i + k - 1L]
    rsum  <- total - lsum
    nL <- k
    nR <- len - k
    diff <- rsum / nR - lsum / nL
    stat[k] <- sqrt((nL * nR)/len) * base::norm(matrix(diff, 1L), "2")
  }
  stat[seq_len(trim)] <- 0
  stat[(len - trim + 1L):len] <- 0
  
  best_rel <- which.max(stat)
  data.frame(est.cp = st_i + best_rel - 1L,
             val = stat[best_rel],
             stringsAsFactors = FALSE)
}


## --------------------------------------------------------------------------
## main function: aggregated/stacked CUSUM
##
## G : array of size c(r1, …, rK, T)
## G_dim : integer vector c(r1, …, rK)
## trim : as in paper; default 2log(T)
## V.diag : use only diagonal of long-run covariance if TRUE
## lrv : apply HAC correction if TRUE
## n : HAC bandwidth; default T^{1/4}
## lbd : minimum interval length for seeded_intervals
## single : if TRUE search on [0,T] only; otherwise on seeded intervals
##
## returns: candidates of estimated change point locations
cand_sbs <- function(G, G_dim, trim = NA, V.diag = TRUE,
                      lrv = TRUE, n = NULL, lbd = NULL, single = FALSE){
  K <- length(G_dim)
  Time <- dim(G)[K + 1L]
  
  if (is.na(trim)) trim <- 2L * round(log(Time))
  trim <- as.integer(trim); stopifnot(trim >= 0L)
  
  if (is.null(n)) n <- floor(Time^0.25)
  
  ## single/multiple
  if (single) {
    intervals <- data.frame(st = 0L, ed = Time)
  } else {
    if (is.null(lbd)) lbd <- round(2 * log(Time)) # can adjust lbd
    intervals <- seeded_intervals(Time, minl = lbd)
  }
  
  D_list <- vector("list", K) 
  
  for (i in seq_len(K)) {
    p_i <- G_dim[i]
    d_i <- p_i * (p_i + 1L) / 2L
    
    FF <- matrix(0, d_i, Time)
    for (t in seq_len(Time)) {
      if (K == 1L) {
        g_t <- G[, t]
        S <- g_t %o% g_t
      } else {
        idx <- c(rep(list(TRUE), K), list(t))
        g_t <- do.call(`[`, c(list(G), idx))
        perm <- c(i, setdiff(seq_len(K), i))
        X_i <- matrix(aperm(g_t, perm), nrow = p_i)
        S <- X_i %*% t(X_i)
      }
      FF[, t] <- vech(S)
    }
    
    Z <- FF - rowMeans(FF)  
    
    ## long-run variance (HAC)
    V <- Z %*% t(Z) / Time
    if (lrv && n >= 1L) {
      for (ell in seq_len(n)) {
        C <- Z[, 1:(Time - ell)] %*% t(Z[, (1:(Time - ell)) + ell]) / Time
        w <- 1 - ell/(n + 1L)
        V <- V + w * (C + t(C))
      }
    }
    if (V.diag) {
      d <- diag(V); d[d <= 0] <- min(d[d > 0])
      Vinvhalf <- diag(1 / sqrt(d))
    } else {
      eg <- eigen(V, symmetric = TRUE)
      vals <- pmax(eg$values, .Machine$double.eps)
      Vinvhalf <- diag(1 / sqrt(vals)) %*% t(eg$vectors)
    }
    
    D_list[[i]] <- Vinvhalf %*% Z # d_i × T
  }
  
  D_stack <- do.call(rbind, D_list) # \sum d_i × T
  
  ## candidate(s)
  cps <- unlist(lapply(seq_len(nrow(intervals)), function(j) {
    iv <- intervals[j, ]
    find_single_cp(D_stack, iv$st, iv$ed, trim)$est.cp
  }))
  sort(unique(cps))
}

