# duanbai.R
# Attention: In duan&bai's method, the input dataset is the transpose of ours

library(parallel)

# Number of Factors Estimation
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



# Newey-West HAC estimator (NW94)
# VAR fit function (order = p)
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
  
  # OLS fit for each column
  Betahat <- matrix(NA, ncol = k, nrow = ncol(X_reg))
  ehat <- matrix(NA, nrow = n - p, ncol = k)
  for (j in 1:k) {
    fit <- lm(Y_response[, j] ~ X_reg - 1)
    Betahat[, j] <- coef(fit)
    ehat[, j] <- residuals(fit)
  }
  list(Betahat = Betahat, ehat = ehat, X_reg = X_reg)
}

# HAC
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



# Brownian motion simulation
simbrownian_v <- function(t, n, dim) {
  steps <- round(t * n)
  w <- matrix(0, nrow = dim, ncol = steps)
  w[,1] <- rnorm(dim) / sqrt(n)
  for (j in 2:steps) {
    w[,j] <- w[,j-1] + rnorm(dim) / sqrt(n)
  }
  w
}

# Critical values for sup-LR test (10%, 5%, 1%)
cv_cache <- list()

critical_value4_lr_HAC_m <- function(F_hat, tau1) {
  Tn <- nrow(F_hat); r0 <- ncol(F_hat)

  FF_I <- matrix(0, nrow = Tn, ncol = r0^2)
  for (t in 1:Tn) {
    M <- tcrossprod(F_hat[t,]) - diag(r0)
    FF_I[t, ] <- as.vector(M)
  }
  Omega_hat <- HAC_NW94_new1(FF_I, rep(1, Tn), cons = 0, kernel = 1, a = 4, pw = 0, p_min = 1, p_max = 1, m = 1.5)
  
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
  cv <- result_sort[idx]
  as.numeric(cv)
}

critical_value4_lr_HAC_m_cached <- function(F_hat, tau1) {
  r0 <- ncol(F_hat)
  key <- paste0("r", r0)
  if (!is.null(cv_cache[[key]])) {
    return(cv_cache[[key]])
  } else {
    cv <- critical_value4_lr_HAC_m(F_hat, tau1)
    cv_cache[[key]] <<- cv
    return(cv)
  }
}

# LR test - single change point detection
# tau: truncation, like a trimming parameter
LR <- function(X, tau1, tau2, r_est = NULL) {
  Tn <- nrow(X);  N <- ncol(X)
  # estimate r by Bai & Ng
  if (is.null(r_est)) {
    nb <- NbFactors2(X, kmax = 10)
    r_est <- nb$khat[1]
  } else{
    r_est <- min(r_est, Tn - 1, N)
  }
  
  eig_all <- if (Tn < N) {
    eigen(X %*% t(X), symmetric = TRUE)
  } else {
    eigen(t(X) %*% X, symmetric = TRUE)
  }
  if (Tn < N) {
    F_hat <- sqrt(Tn) * eig_all$vectors[, 1:r_est, drop = FALSE]
  } else {
    loadings <- sqrt(N) * eig_all$vectors[, 1:r_est, drop = FALSE]
    F_hat <- X %*% loadings / N
  }
  
  # QML location
  T_min <- round(Tn * tau1);  T_max <- round(Tn * tau2)
  SSR <- numeric(T_max - T_min + 1)
  for (k in T_min:T_max) {
    S1 <- crossprod(F_hat[1:k, , drop = FALSE]) / k
    S2 <- crossprod(F_hat[(k+1):Tn, , drop = FALSE]) / (Tn - k)
    SSR[k - T_min + 1] <- k * log(det(S1)) + (Tn - k) * log(det(S2))
  }
  hatk <- which.min(SSR) - 1 + T_min
  
  # LR statistic
  S1 <- crossprod(F_hat[1:hatk, , drop = FALSE]) / hatk
  S2 <- crossprod(F_hat[(hatk+1):Tn, , drop = FALSE]) / (Tn - hatk)
  S0 <- crossprod(F_hat) / Tn
  LR_stat <- Tn * log(det(S0)) - hatk * log(det(S1)) - (Tn - hatk) * log(det(S2))
  
  # critical values and 5% decision
  # cv <- critical_value4_lr_HAC_m(F_hat, tau1)
  cv <- critical_value4_lr_HAC_m_cached(F_hat, tau1)
  Result_rej1 <- as.integer(LR_stat > cv[2])
  
  list(Result_rej1 = Result_rej1, hatk = hatk)
}




# LR test - multiple change points detection
recurse_LR <- function(Xsub, tau1 = 0.3, offset = 0, min_size = 20, r_est = NULL) {
  if (nrow(Xsub) < min_size) return(list(n = 0, cps = integer(0)))
  out <- LR(Xsub, tau1, 1 - tau1, r_est)
  if (out$Result_rej1 == 0) return(list(n = 0, cps = integer(0)))
  cp <- out$hatk
  left <- Xsub[1:cp, , drop = FALSE]
  right <- Xsub[(cp+1):nrow(Xsub), , drop = FALSE]
  me <- 1
  my_cps <- offset + cp
  a1 <- recurse_LR(left, tau1, offset, min_size, r_est)
  a2 <- recurse_LR(right, tau1, offset + cp, min_size, r_est)
  return(list(n = me + a1$n + a2$n,
              cps = sort(c(my_cps, a1$cps, a2$cps))))
}




# source("duanbai_update_v2.R")



