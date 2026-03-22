library(lubridate)
library(dplyr)
library(kableExtra)

# factor Number
NbFactors2 <- function(X, kmax = NULL) {
  Tn <- nrow(X); N <- ncol(X)
  if (is.null(kmax)) kmax <- min(Tn, N)
  kmax_use <- min(kmax, Tn, N)
  if (kmax_use < 1) stop("Segment too small for factor number estimation.")
  
  if (Tn < N) {
    eig_all <- eigen(X %*% t(X), symmetric = TRUE)
    factors_kmax <- sqrt(Tn) * eig_all$vectors[, 1:kmax_use, drop = FALSE] 
    loadings_kmax <- t(X) %*% factors_kmax / Tn
    betahat_kmax <- loadings_kmax %*% chol(t(loadings_kmax) %*% loadings_kmax / N)
  } else {
    eig_all <- eigen(t(X) %*% X, symmetric = TRUE)
    loadings_kmax <- sqrt(N) * eig_all$vectors[, 1:kmax_use, drop = FALSE] 
    factors_kmax <- X %*% loadings_kmax / N
    betahat_kmax <- (t(X) %*% X) %*% loadings_kmax / (N * Tn)
  }
  
  Z_kmax <- X - X %*% betahat_kmax %*% solve(t(betahat_kmax) %*% betahat_kmax) %*% t(betahat_kmax)
  var_Z_kmax <- sum(Z_kmax^2) / (N * Tn)
  
  V <- numeric(kmax_use)
  for (k in seq_len(kmax_use)) {
    if (Tn < N) {
      eig_k <- eigen(X %*% t(X), symmetric = TRUE)
      factors_k <- sqrt(Tn) * eig_k$vectors[, 1:k, drop = FALSE] 
      loadings_k <- t(X) %*% factors_k / Tn
      betahat_k <- loadings_k %*% chol(t(loadings_k) %*% loadings_k / N)
    } else {
      eig_k <- eigen(t(X) %*% X, symmetric = TRUE)
      loadings_k <- sqrt(N) * eig_k$vectors[, 1:k, drop = FALSE] 
      factors_k <- X %*% loadings_k / N
      betahat_k <- (t(X) %*% X) %*% loadings_k / (N * Tn)
    }
    Zk <- X - X %*% betahat_k %*% solve(t(betahat_k) %*% betahat_k) %*% t(betahat_k)
    V[k] <- sum(Zk^2) / (N * Tn)
  }
  
  CNT <- min(N, Tn)
  Penalty <- c((N + Tn)/(N * Tn) * log((N * Tn)/(N + Tn)),
               (N + Tn)/(N * Tn) * log(CNT),
               log(CNT)/CNT)
  Penalty <- matrix(rep(Penalty, each = kmax_use), nrow = kmax_use)
  
  IC <- matrix(NA, nrow = kmax_use + 1, ncol = 4)
  PC <- matrix(NA, nrow = kmax_use + 1, ncol = 4)
  AIC3 <- matrix(NA, nrow = kmax_use + 1, ncol = 2)
  BIC3 <- matrix(NA, nrow = kmax_use + 1, ncol = 2)
  
  mean_X2 <- sum(X^2) / (N * Tn)
  IC[1, 1] <- 0; IC[1, 2:4] <- log(mean_X2)
  PC[1, 1] <- 0; PC[1, 2:4] <- mean_X2
  AIC3[1, 1] <- 0; AIC3[1, 2] <- mean_X2
  BIC3[1, 1] <- 0; BIC3[1, 2] <- mean_X2
  
  kk <- matrix(1:kmax_use, nrow = kmax_use, ncol = 3)
  PC[2:(kmax_use + 1), 1] <- 1:kmax_use
  IC[2:(kmax_use + 1), 1] <- 1:kmax_use
  PC[2:(kmax_use + 1), 2:4] <- matrix(rep(V, 3), ncol = 3) + var_Z_kmax * kk * Penalty
  IC[2:(kmax_use + 1), 2:4] <- log(matrix(rep(V, 3), ncol = 3)) + kk * Penalty
  BIC3[2:(kmax_use + 1), 1] <- 1:kmax_use
  BIC3[2:(kmax_use + 1), 2] <- V + (1:kmax_use) * var_Z_kmax * (N + Tn - (1:kmax_use)) * log(N * Tn) / (N * Tn)
  AIC3[2:(kmax_use + 1), 1] <- 1:kmax_use
  AIC3[2:(kmax_use + 1), 2] <- V + (1:kmax_use) * var_Z_kmax * (N + Tn - (1:kmax_use)) * 2 / (N * Tn)
  
  khat_IC <- apply(IC[, 2:4, drop = FALSE], 2, function(col) which.min(col) - 1)
  khat_PC <- apply(PC[, 2:4, drop = FALSE], 2, function(col) which.min(col) - 1)
  khat_AIC3 <- which.min(AIC3[, 2]) - 1
  khat_BIC3 <- which.min(BIC3[, 2]) - 1
  khat <- c(khat_IC, khat_PC, khat_AIC3, khat_BIC3)
  list(
    khat = khat,
    khat_IC = khat_IC,
    khat_PC = khat_PC,
    khat_AIC3 = khat_AIC3,
    khat_BIC3 = khat_BIC3,
    kmax = kmax_use,
    IC = IC,
    PC = PC,
    AIC3 = AIC3,
    BIC3 = BIC3,
    Vkmax = var_Z_kmax
  )
}


panelFactorNew_R <- function(X, r) {
  Tn <- nrow(X)
  N  <- ncol(X)
  
  XX <- X %*% t(X) / (N * Tn)
  sv <- svd(XX)
  
  factor <- sv$u[, 1:r, drop = FALSE] * sqrt(Tn)
  lambda <- t(X) %*% factor / Tn
  VNT <- diag(sv$d[1:r], nrow = r)
  
  list(factor = factor, lambda = lambda, VNT = VNT)
}

# VAR/HAC
fitVAR <- function(Y, p_min = 1, p_max = 1) {
  p <- p_min
  n <- nrow(Y)
  k <- ncol(Y)
  if (p < 1) stop("VAR order must be >= 1.")
  
  Ymat <- NULL
  for (i in 1:p) {
    Ymat <- cbind(Ymat, Y[(p - i + 1):(n - i), ])
  }
  Y_response <- Y[(p + 1):n, ]
  X_reg <- cbind(1, Ymat)
  
  Betahat <- matrix(NA, ncol = k, nrow = ncol(X_reg))
  ehat <- matrix(NA, nrow = n - p, ncol = k)
  for (j in 1:k) {
    fit <- lm(Y_response[, j] ~ X_reg - 1)
    Betahat[, j] <- coef(fit)
    ehat[, j] <- residuals(fit)
  }
  list(Betahat = Betahat, ehat = ehat, X_reg = X_reg)
}

HAC_NW94_new1 <- function(X, e, cons = 0, kernel = 1, a = 4, pw = 0, p_min = 1, p_max = 1, m = 1.5) {
  N <- nrow(X); k <- ncol(X)
  w <- if (cons == 0) rep(1, k) else c(0, rep(1, k - 1))
  Y <- X * matrix(e, nrow = N, ncol = k)
  
  if (pw == 0) {
    Z <- Y
    Betahat <- NULL
  } else if (pw == 1) {
    varfit <- fitVAR(Y, p_min, p_max)
    ehat <- varfit$ehat
    Z <- ehat
    Betahat <- varfit$Betahat[-1, , drop = FALSE]
  }
  
  Tn <- nrow(Z)
  
  if (kernel == 0) {
    V <- crossprod(Z) / Tn
  } else if (kernel == 1) {
    q <- 1
    cr <- 1.1447
    nlag <- floor(a * (Tn / 100)^(2 / 9))
    nlag <- min(nlag, Tn - 1)  
    if (nlag < 0) nlag <- 0      
    s <- numeric(nlag + 1)
    s0 <- 0
    sq <- 0
    for (v in 0:nlag) {
      if ((Tn - v) < 1 || (v + 1) > Tn) next  
      Gamma_v <- t(Z[1:(Tn - v), , drop = FALSE]) %*% Z[(v + 1):Tn, , drop = FALSE] / Tn
      s[v + 1] <- as.numeric(t(w) %*% Gamma_v %*% w)
      s0 <- s0 + 2 * s[v + 1]
      sq <- sq + 2 * v^q * s[v + 1]
    }
    s0 <- s0 - s[1]
    rhat <- m * cr * ((sq / s0)^2)^(1 / (2 * q + 1))
    ST <- floor(rhat * Tn^(1 / (2 * q + 1)))
    ST <- min(ST, Tn - 1) 
    if (ST < 0) ST <- 0   
    Gama <- vector("list", ST + 1)
    Phihat <- matrix(0, k, k)
    for (j in 0:ST) {
      if ((Tn - j) < 1 || (j + 1) > Tn) next  
      Gama[[j + 1]] <- t(Z[1:(Tn - j), , drop = FALSE]) %*% Z[(j + 1):Tn, , drop = FALSE] / Tn
      Phihat <- Phihat + (1 - j / (ST + 1)) * (Gama[[j + 1]] + t(Gama[[j + 1]]))
    }
    V <- Phihat - Gama[[1]]
  } else {
    stop("Only kernel = 0 (White) or kernel = 1 (Bartlett) are supported.")
  }
  
  if (pw == 0) {
    output <- V
  } else if (pw == 1) {
    D <- diag(k)
    for (i in seq_len(nrow(Betahat) / k)) {
      idx <- ((i - 1) * k + 1):(i * k)
      D <- D - Betahat[idx, ]
    }
    output <- solve(D) %*% V %*% t(solve(D))
  }
  output
}

simbrownian_v <- function(t, n, dim) {
  steps <- round(t * n)
  w <- matrix(0, nrow = dim, ncol = steps)
  w[,1] <- rnorm(dim) / sqrt(n)
  for (j in 2:steps) {
    w[,j] <- w[,j-1] + rnorm(dim) / sqrt(n)
  }
  w
}


critical_value4_lr_HAC_m <- function(F_hat, tau1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  Tn <- nrow(F_hat)
  r0 <- ncol(F_hat)
  
  FF_I <- matrix(0, nrow = Tn, ncol = r0^2)
  for (t in seq_len(Tn)) {
    M <- tcrossprod(F_hat[t, ]) - diag(r0)
    FF_I[t, ] <- as.vector(M)
  }
  
  Omega_hat <- HAC_NW94_new1(
    FF_I, rep(1, Tn),
    cons = 0, kernel = 1, a = 4, pw = 0,
    p_min = 1, p_max = 1, m = 1.5
  )
  
  n <- 2000
  niter <- 1000
  result <- numeric(niter)
  j_range <- round(tau1 * n):round((1 - tau1) * n)
  
  for (i in seq_len(niter)) {
    W <- simbrownian_v(1, n, r0^2)
    LR <- rep(NA_real_, n)
    for (j in j_range) {
      diffW <- W[, j] - (j / n) * W[, n]
      LR[j] <- crossprod(diffW, Omega_hat %*% diffW) / (2 * (j / n) * (1 - j / n))
    }
    result[i] <- max(LR[j_range], na.rm = TRUE)
  }
  
  result_sort <- sort(result)
  idx <- c(floor(0.90 * niter), floor(0.95 * niter), floor(0.99 * niter))
  as.numeric(result_sort[idx])
}

critical_value4_lr_HAC_m_cached <- function(F_hat, tau1, seed = NULL, cv_cache) {
  r0 <- ncol(F_hat)
  key <- paste0("r", r0, "_tau", tau1, "_seed", if (is.null(seed)) "NULL" else seed)
  
  if (!is.null(cv_cache[[key]])) {
    return(list(cv = cv_cache[[key]], cv_cache = cv_cache))
  }
  
  cv <- critical_value4_lr_HAC_m(F_hat, tau1, seed = seed)
  cv_cache[[key]] <- cv
  
  list(cv = cv, cv_cache = cv_cache)
}

# New: global-trim BS for LR
estimate_r_global <- function(X, kmax = 10, ic_col = 1) {
  nb <- NbFactors2(X, kmax = kmax)
  r_est <- nb$khat[ic_col]
  max(1, as.integer(r_est))
}

logdet_safe <- function(S, eps = 1e-6) {
  k <- nrow(S)
  val <- tryCatch({
    cholS <- chol(S + diag(eps, k))
    2 * sum(log(diag(cholS)))
  }, error = function(e) NA_real_)
  if (is.na(val)) {
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    ev <- pmax(ev, eps)
    sum(log(ev))
  } else val
}

LR_globaltrim_once <- function(X, s, e, min_gap, r_est,
                               tau_for_cv = 0.1, eps_ridge = 1e-6,
                               seed = NULL, cv_cache) {
  Tn <- nrow(X)
  N <- ncol(X)
  
  if ((e - s + 1) < 2 * min_gap + 1) {
    return(list(reject = 0L, k_hat = NA_integer_, lr = -Inf, cv_cache = cv_cache))
  }
  
  k_lo <- s + min_gap
  k_hi <- e - min_gap
  
  Xseg <- X[s:e, , drop = FALSE]
  Tseg <- nrow(Xseg)
  
  pf <- panelFactorNew_R(Xseg, r_est)
  F_hat <- pf$factor
  
  idx_lo <- k_lo - s + 1L
  idx_hi <- k_hi - s + 1L
  if (idx_lo < 1L || idx_hi > Tseg || idx_lo >= idx_hi) {
    return(list(reject = 0L, k_hat = NA_integer_, lr = -Inf, cv_cache = cv_cache))
  }
  
  best_val <- Inf
  best_k_local <- NA_integer_
  for (k_local in idx_lo:idx_hi) {
    S1 <- crossprod(F_hat[1:k_local, , drop = FALSE]) / k_local
    S2 <- crossprod(F_hat[(k_local + 1):Tseg, , drop = FALSE]) / (Tseg - k_local)
    v <- k_local * logdet_safe(S1, eps_ridge) + (Tseg - k_local) * logdet_safe(S2, eps_ridge)
    if (v < best_val) {
      best_val <- v
      best_k_local <- k_local
    }
  }
  
  k_hat_global <- s - 1L + best_k_local
  
  S1 <- crossprod(F_hat[1:best_k_local, , drop = FALSE]) / best_k_local
  S2 <- crossprod(F_hat[(best_k_local + 1):Tseg, , drop = FALSE]) / (Tseg - best_k_local)
  S0 <- crossprod(F_hat) / Tseg
  
  LR_stat <- Tseg * logdet_safe(S0, eps_ridge) -
    best_k_local * logdet_safe(S1, eps_ridge) -
    (Tseg - best_k_local) * logdet_safe(S2, eps_ridge)
  
  tmp <- critical_value4_lr_HAC_m_cached(F_hat, tau_for_cv, seed = seed, cv_cache = cv_cache)
  cv <- tmp$cv
  cv_cache <- tmp$cv_cache
  reject <- as.integer(LR_stat > cv[3]) 
  
  list(reject = reject, k_hat = k_hat_global, lr = LR_stat, cv_cache = cv_cache)
}

bs_LR_globaltrim <- function(X, tau_global = 0.1, min_size = 20, r_est = NULL,
                             tau_for_cv = 0.1, max_cps = NULL, eps_ridge = 1e-6,
                             seed = NULL) {
  Tn <- nrow(X)
  min_gap <- max(1L, floor(tau_global * Tn))
  
  if (is.null(r_est)) r_est <- estimate_r_global(X, kmax = 10, ic_col = 1)
  if (is.null(max_cps)) {
    max_cps <- max(0L, floor((Tn - 1) / min_gap) - 1L)
  }
  
  cv_cache <- list()
  
  seg_stack <- list(c(1L, Tn))
  cps <- integer(0)
  lr_at_cp <- numeric(0)
  
  while (length(seg_stack) > 0 && length(cps) < max_cps) {
    seg <- seg_stack[[length(seg_stack)]]
    seg_stack <- seg_stack[-length(seg_stack)]
    
    s <- seg[1]
    e <- seg[2]
    
    if ((e - s + 1) < max(min_size, 2 * min_gap + 1)) next
    
    res <- LR_globaltrim_once(
      X, s, e, min_gap, r_est,
      tau_for_cv = tau_for_cv,
      eps_ridge = eps_ridge,
      seed = seed,
      cv_cache = cv_cache
    )
    cv_cache <- res$cv_cache
    
    if (res$reject == 1L && !is.na(res$k_hat)) {
      k <- res$k_hat
      if (length(cps) == 0 || all(abs(k - cps) >= min_gap)) {
        cps <- c(cps, k)
        lr_at_cp <- c(lr_at_cp, res$lr)
        
        if ((k - s + 1) >= max(min_size, 2 * min_gap + 1)) {
          seg_stack[[length(seg_stack) + 1]] <- c(s, k)
        }
        if ((e - k) >= max(min_size, 2 * min_gap + 1)) {
          seg_stack[[length(seg_stack) + 1]] <- c(k + 1L, e)
        }
      }
    }
  }
  
  if (length(cps) > 1) {
    o <- order(cps)
    cps <- cps[o]
    lr_at_cp <- lr_at_cp[o]
    
    keep <- rep(TRUE, length(cps))
    i <- 1
    while (i < length(cps)) {
      j <- i + 1
      while (j <= length(cps) && (cps[j] - cps[i]) < min_gap) j <- j + 1
      if (j - i > 1) {
        block <- i:(j - 1)
        best <- block[which.max(lr_at_cp[block])]
        keep[block] <- FALSE
        keep[best] <- TRUE
      }
      i <- j
    }
    cps <- cps[keep]
  }
  
  list(
    n = length(cps),
    cps = sort(unique(as.integer(cps))),
    min_gap = min_gap,
    r_est = r_est
  )
}


