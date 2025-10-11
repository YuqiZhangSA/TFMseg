op_norm <- function(M) {
  d <- svd(M, nu = 0L, nv = 0L)$d
  if (length(d)) d[1L] else 0
}

mat_norm <- function(M, type = c("op","fro","max")) {
  type <- match.arg(type)
  switch(type,
         op  = op_norm(M),
         fro = sqrt(sum(M * M)),
         max = max(abs(M))
  )
}

scalar_from_cov <- function(C, type = c("trace","op")) {
  type <- match.arg(type)
  if (type == "trace") sum(diag(C)) else op_norm(C)
}

`%||%` <- function(x, y) if (is.null(x)) y else x


mode_cov_from_tensor <- function(Z, mode_k, a, b, centre = FALSE) {
  dims <- dim(Z); K <- length(dims) - 1L
  Ttot <- dims[K + 1L]
  if (a < 0L || b > Ttot || a >= b) stop("Invalid (a,b] window.")
  pk <- dims[mode_k]
  C  <- matrix(0, pk, pk)
  
  for (t in (a + 1L):b) {
    idx   <- c(rep(list(TRUE), K), list(t))
    slice <- do.call(`[`, c(list(Z), idx, list(drop = FALSE)))
    slice <- array(slice, dim = dims[1:K])
    Mk    <- matrix(aperm(slice, c(mode_k, setdiff(seq_len(K), mode_k))), nrow = pk)
    if (centre) Mk <- sweep(Mk, 1L, rowMeans(Mk), "-")
    C <- C + Mk %*% t(Mk)
  }
  C / (b - a)
}


m11_pair_general <- function(C_pre, C_post,
                             xi_norm = c("op","fro","max"),
                             scalar_type = c("trace","op"),
                             extra_div = c("none","pre","post","both","jump")) {
  xi_norm <- match.arg(xi_norm)
  scalar_type <- match.arg(scalar_type)
  extra_div <- match.arg(extra_div)
  
  s0 <- scalar_from_cov(C_pre, scalar_type)
  s1 <- scalar_from_cov(C_post, scalar_type)
  
  if (!is.finite(s0) || !is.finite(s1) || s0 <= 0 || s1 <= 0) {
    return(list(zeta = NA_real_, Xi = NA,
                tr_pre = sum(diag(C_pre)), tr_post = sum(diag(C_post))))
  }
  
  Xi <- (C_post / s1) - (C_pre / s0)
  z  <- mat_norm(Xi, xi_norm)
  
  if (extra_div != "none") {
    denom <- switch(extra_div,
                    pre = op_norm(C_pre),
                    post = op_norm(C_post),
                    both = sqrt(op_norm(C_pre) * op_norm(C_post)),
                    jump = op_norm(C_post - C_pre)
    )
    if (!is.finite(denom) || denom <= 0) {
      return(list(zeta = NA_real_, Xi = Xi,
                  tr_pre = sum(diag(C_pre)), tr_post = sum(diag(C_post))))
    }
    z <- z / denom
  }
  
  list(zeta = z, Xi = Xi, tr_pre = sum(diag(C_pre)), tr_post = sum(diag(C_post)))
}


#====================================================

m11_zeta_from_X <- function(
    X, dim_obs, r_hat = NULL,
    cp_vec,
    st = 1L,
    ed = dim(X)[length(dim_obs) + 1L],
    centre_cov = FALSE,
    xi_norm = c("op","fro","max"),
    scalar_type = c("trace","op"),
    extra_div = c("none","pre","post","both","jump"),
    return_Xi = FALSE,
    centre = TRUE,
    G_hac = FALSE,    
    hac_n = NULL
) {
  xi_norm <- match.arg(xi_norm)
  scalar_type <- match.arg(scalar_type)
  extra_div <- match.arg(extra_div)
  
  # 1) estimate pseudo-factors G 
  dim_X <- as.integer(dim_obs)
  est <- if (is.null(r_hat)) {
    global_pca(X = X, dim_X = dim_X, centre = centre, proj = TRUE)
  } else {
    global_pca(X = X, dim_X = dim_X, r_hat = r_hat, centre = centre, proj = TRUE)
  }
  G <- est$G_proj                              
  G_dim <- as.vector(est$r_hat)
  K <- length(G_dim)
  Ttot <- dim(G)[K + 1L]
  
  if (ed > Ttot) stop("`ed` exceeds time dimension of G.")
  
  cp_vec <- sort(unique(as.integer(cp_vec)))
  cp_vec <- cp_vec[cp_vec > st & cp_vec < ed]
  if (!length(cp_vec)) stop("No interior change points between st and ed.")
  
  bnds <- c(st - 1L, cp_vec, ed) 
  Mseg <- length(bnds) - 1L     
  
  # 3) covariances of G
  Cov <- vector("list", Mseg)
  
  W_list <- vector("list", K)
  if (isTRUE(G_hac)) {
    for (k in seq_len(K)) {
      Sigma_k <- .hac_longrun_for_mode(G, mode_k = k, n = hac_n %||% NULL, centre_rows = TRUE)
      W_list[[k]] <- .psd_inv_sqrt(Sigma_k)
    }
  }
  
  for (m in seq_len(Mseg)) {
    Cm <- vector("list", K)
    a_m <- bnds[m]; b_m <- bnds[m + 1L]
    for (k in seq_len(K)) {
      pk <- dim(G)[k]
      Ck <- matrix(0, pk, pk)
      
      for (t in (a_m + 1L):b_m) {
        idx <- c(rep(list(TRUE), K), list(t))
        Gt <- do.call(`[`, c(list(G), idx, list(drop = FALSE)))
        Gt <- array(Gt, dim = dim(G)[1:K])
        
        perm <- c(k, setdiff(seq_len(K), k))
        Mk <- matrix(aperm(Gt, perm), nrow = pk)
        
        if (centre_cov) Mk <- sweep(Mk, 1L, rowMeans(Mk), "-")
        
        # HAC standardisation (left-multiply by W_k)
        if (isTRUE(G_hac)) Mk <- W_list[[k]] %*% Mk
        
        Ck <- Ck + Mk %*% t(Mk)
      }
      
      Cm[[k]] <- Ck / (b_m - a_m)
    }
    Cov[[m]] <- Cm
  }
  
  
  qhat <- length(cp_vec)
  rows <- vector("list", qhat * K)
  Xi_keep <- if (return_Xi) vector("list", qhat * K) else NULL
  zeta_mat <- matrix(NA_real_, nrow = K, ncol = qhat,
                     dimnames = list(paste0("mode", seq_len(K)),
                                     paste0("cp", seq_len(qhat))))
  idx <- 0L
  
  for (j in seq_len(qhat)) {
    a_pre <- bnds[j]; b_pre <- bnds[j + 1L]
    a_post <- bnds[j + 1L]; b_post <- bnds[j + 2L]
    for (k in seq_len(K)) {
      idx <- idx + 1L
      res <- m11_pair_general(
        Cov[[j]][[k]], Cov[[j + 1L]][[k]],
        xi_norm = xi_norm, scalar_type = scalar_type, extra_div = extra_div
      )
      zeta_mat[k, j] <- res$zeta
      rows[[idx]] <- data.frame(
        cp_index = j,
        cp_time  = cp_vec[j],
        mode = k,
        seg_pre = sprintf("(%d,%d]", a_pre,  b_pre),
        seg_post = sprintf("(%d,%d]", a_post, b_post),
        len_pre = b_pre - a_pre,
        len_post = b_post - a_post,
        zeta = res$zeta,
        xi_norm = xi_norm,
        scalar = scalar_type,
        extra_div = extra_div,
        stringsAsFactors = FALSE
      )
      if (return_Xi) Xi_keep[[idx]] <- res$Xi
    }
  }
  
  tab <- do.call(rbind, rows)
  
  out <- list(
    table = tab,          
    zeta_mat = zeta_mat,    
    cp_vec = cp_vec,
    G_dim = G_dim,
    xi_norm = xi_norm,
    scalar = scalar_type,
    extra_div = extra_div
  )
  if (return_Xi) out$Xi <- Xi_keep
  out
}



.psd_inv_sqrt <- function(S, eps = .Machine$double.eps) {
  E <- eigen((S + t(S))/2, symmetric = TRUE)
  vals <- pmax(E$values, eps)
  E$vectors %*% diag(1/sqrt(vals), length(vals)) %*% t(E$vectors)
}

.hac_longrun_for_mode <- function(G, mode_k, n = NULL, centre_rows = TRUE) {
  dims <- dim(G); K <- length(dims) - 1L; Ttot <- dims[K + 1L]
  rk <- dims[mode_k]; q <- prod(dims[setdiff(seq_len(K), mode_k)])
  if (is.null(n)) n <- max(1L, floor(Ttot^0.25))
  
  M <- vector("list", Ttot); row_mean <- matrix(0, rk, q)
  for (t in seq_len(Ttot)) {
    idx <- c(rep(list(TRUE), K), list(t))
    Gt  <- do.call(`[`, c(list(G), idx, list(drop = FALSE)))
    Gt  <- array(Gt, dim = dims[1:K])
    perm <- c(mode_k, setdiff(seq_len(K), mode_k))
    Mk <- matrix(aperm(Gt, perm), nrow = rk)
    M[[t]] <- Mk; row_mean <- row_mean + Mk
  }
  if (centre_rows) {
    row_mean <- row_mean / Ttot
    for (t in seq_len(Ttot)) M[[t]] <- M[[t]] - row_mean
  }
  
  Sigma <- Reduce(`+`, lapply(M, function(A) A %*% t(A))) / Ttot

  if (n >= 1L) {
    for (ell in seq_len(min(n, Ttot - 1L))) {
      w <- 1 - ell / (n + 1)
      acc <- matrix(0, rk, rk)
      for (t in (ell + 1L):Ttot) acc <- acc + M[[t]] %*% t(M[[t - ell]])
      acc <- acc / (Ttot - ell)
      Sigma <- Sigma + w * (acc + t(acc))
    }
  }
  Sigma
}

