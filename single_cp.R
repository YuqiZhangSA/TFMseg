#' Single‐CP via CUSUM + HAC
#'
#' @param D    Numeric matrix of size d \times T.
#' @param st   Start index, 0 <= st < T.
#' @param ed   End index,   st < ed <= T.
#' @param trim Positive integer; zeros CUSUM in the first/last \code{trim} points.
#' @return A 1‐row data.frame with columns \code{est.cp,val,st,ed,trim}.
#' @export
#' 

find_single_cp <- function(D, st, ed, trim) {
  Time <- ncol(D)
  st_i <- max(1, as.integer(st) + 1)
  ed_i <- min(Time, as.integer(ed))
  len <- ed_i - st_i + 1
  stopifnot(len > 2 * trim + 1)
  
  stat <- numeric(len)
  total <- rowSums(D[, st_i:ed_i, drop = FALSE])
  lsum <- numeric(nrow(D))
  
  for (k in seq_len(len - 1)) {
    lsum <- lsum + D[, st_i + k - 1]
    rsum <- total - lsum
    nL <- k
    nR <- len - k
    diff <- rsum / nR - lsum / nL
    stat[k] <- sqrt((nL * nR) / len) * sqrt(sum(diff^2))
  }
  stat[seq_len(trim)] <- 0
  stat[(len - trim + 1):len] <- 0
  
  best_rel <- which.max(stat)
  best_abs <- st_i + best_rel - 1
  
  data.frame(
    est.cp = best_abs,
    val = stat[best_rel],
    st = st,
    ed = ed,
    trim = trim,
    stringsAsFactors = FALSE
  )
}




#' Tensor CUSUM candidate generation
#'
#' @param G      Array of size c(r1, …, rK, T) of pseudo‐factors.
#' @param G_dim  Integer vector c(r1, …, rK).
#' @param trim   Non‐negatintervale integer; if NA, defaults to round(2*log(T)).
#' @param V.diag Logical; if TRUE, only diagonal of long‐run covariance is used.
#' @param lrv    Logical; if TRUE, apply HAC up to lag n.
#' @param n      HAC bandwidth; defaults to floor(T^{1/4}).
#' @param lbd    Min interval length.
#'
#' @return A data.frame with columns mode,est.cp,val,st,ed,trim — one row per mode×interval.
#' @export
#' 

cand_sbs <- function(G, G_dim, trim   = NA, V.diag = TRUE, lrv = TRUE, n  = NULL, lbd = NULL, single = FALSE) {
  K <- length(G_dim)
  Time <- dim(G)[K+1]
  
  if (is.na(trim)) trim <- 2 * round(log(Time))
  trim <- as.integer(trim)
  stopifnot(trim >= 0)
  
  if (is.null(n)) n <- floor(Time^0.25)
  
  if (single) {
    interval <- data.frame(st = 0, ed = Time)
  } else {
    if (is.null(lbd)) lbd <- round(4*log(Time))
    interval <- seeded_intervals(Time, minl = lbd)
  }
  
  cusum_paths <- list()

  # drop too‐short intervals
  #interval <- interval[(interval$ed - interval$st) > 2*trim, , drop = FALSE]
  #if (!nrow(interval)) stop("No intervals longer than 2*trim; adjust trim/lbd.")
  
  vech <- function(M) M[lower.tri(M, diag = TRUE)]
  
  all_cands <- vector("list", K)
  
  for (i in seq_len(K)) {
    p_i <- G_dim[i]
    d_i <- p_i*(p_i+1)/2
    
    FF <- matrix(0, d_i, Time)
    for (t in seq_len(Time)) {
      if (K == 1) {
        # vector 
        g_t <- G[, t]
        S <- g_t %o% g_t
      } else {
        # matrix/tensor 
        idx <- c(rep(list(TRUE), K), list(t))
        g_t <- do.call(`[`, c(list(G), idx))
        perm<- c(i, setdiff(seq_len(K), i))
        X_i <- matrix(aperm(g_t, perm), nrow = p_i)
        S <- X_i %*% t(X_i)
      }
      FF[,t] <- vech(S)
    }
    mu_FF <- rowMeans(FF)
    Z <- FF - mu_FF
    
    # HAC long‐run variance
    Vmat <- Z %*% t(Z) / Time
    if (lrv && n >= 1) {
      for (ell in seq_len(n)) {
        C <- Z[, 1:(Time-ell)] %*% t(Z[,(1:(Time-ell))+ell]) / Time
        w <- 1 - ell/(n+1)
        Vmat <- Vmat + w*(C + t(C))
      }
    }
    if (V.diag) {
      dd <- diag(Vmat); dd[dd<=0] <- min(dd[dd>0])
      Vhalf_inv <- diag(1/sqrt(dd))
    } else {
      eg <- eigen(Vmat, symmetric=TRUE)
      vals <- pmax(eg$values, .Machine$double.eps)
      Vhalf_inv <- diag(1/sqrt(vals)) %*% t(eg$vectors)
    }
    Dstd <- Vhalf_inv %*% Z
    
    # Compute full CUSUM path across [1..Time]
    stat <- numeric(Time)
    total <- rowSums(Dstd)
    lsum <- numeric(nrow(Dstd))
    
    for (k in 1:(Time - 1)) {
      lsum <- lsum + Dstd[, k]
      rsum <- total - lsum
      nL <- k
      nR <- Time - k
      diff <- rsum / nR - lsum / nL
      stat[k] <- sqrt((nL * nR) / Time) * sqrt(sum(diff^2))
    }
    
    # trimming
    stat[1:trim] <- 0
    stat[(Time - trim + 1):Time] <- 0
    cusum_paths[[paste0("mode_", i)]] <- stat
    
    
    cands_i <- lapply(seq_len(nrow(interval)), function(j) {
      s <- interval$st[j]
      e <- interval$ed[j]
      df <- find_single_cp(Dstd, s, e, trim)
      cbind(mode = i, df)
    })
    cdf <- do.call(rbind, cands_i)
    
    cdf <- cdf[!duplicated(cdf$est.cp), , drop = FALSE]
    rownames(cdf) <- NULL
    all_cands[[i]] <- cdf
  }
  
  result <- do.call(rbind, all_cands)
  rownames(result) <- NULL
  
  attr(result, "cusum_paths") <- cusum_paths
  class(result) <- c("cand_sbs_result", class(result))
  return(result)
  
  
}


