simulate_m11_zeta <- function(
    nrep,
    data_setting = c("s1","s2","s3"),
    dim_obs, dim_latent,
    Time,
    dep = TRUE,
    dist = c("Gaussian","heavy"),
    theta_coef = c(0.2, 0.4, 0.6, 0.8),
    r_hat = NULL,                  
    cp_source = c("true","estimate"),
    threshold_coef = NULL, 
    V.diag = TRUE, lrv = TRUE,
    xi_norm = c("op","fro","max"),
    scalar_type = c("trace","op"),
    extra_div = c("none","pre","post","both","jump"),
    G_hac = FALSE,
    seed_start = 888,
    reverse = FALSE 
) {
  data_setting <- match.arg(data_setting)
  cp_source    <- match.arg(cp_source)
  dist         <- match.arg(dist)
  xi_norm      <- match.arg(xi_norm)
  scalar_type  <- match.arg(scalar_type)
  extra_div    <- match.arg(extra_div)
  
  if (!all(theta_coef == 0)) {
    theta <- floor(Time * theta_coef)
    theta <- sort(unique(pmax(1L, pmin(Time - 1L, theta))))
  } else {
    theta <- integer(0)
  }
  m <- length(theta)
  
  res_list <- vector("list", nrep)
  
  for (sim in seq_len(nrep)) {
    set.seed(seed_start + sim)
    
    if (data_setting == "s1") {
      type_mode <- list(
        c(1),          # cp1 unchanged
        c(2, 3),       # cp2: modes 2 & 3
        c(1),          # cp3: mode 1
        c(3, 1)        # cp4: modes 3 & 1
      )
      type_change <- list(
        "l",                 # cp1 unchanged
        c("f","l"),          # cp2: factor-number + loading change
        "f",                 # cp3: factor-number change
        c("l","l")           # cp4: loading changes
      )
      
      make_I  <- function(r) diag(r)
      make_3I <- function(r) 3 * diag(r)
      make_C2 <- function(r) { M <- diag(r); M[r,r] <- 0; M }
      make_Orth <- function(r, angle = pi/5) {
        if (r <= 1) return(diag(r))
        R2 <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
        M  <- diag(r); M[1:2, 1:2] <- R2; M
      }
      make_C3 <- function(r) {
        vals <- pmax(1 - 0.2 * (0:(r - 1)), 0.1)
        diag(vals)
      }
      
      transform_list <- list(
        list(make_I,   NULL,   NULL),   # cp1: I on mode 1
        list(NULL,     NULL,   make_3I),# cp2: 3I on mode 3
        NULL,                           # cp3: f only
        list(make_Orth,NULL,   make_C2) # cp4: orth on mode 1, C2 on mode 3
      )
      add_factors <- list(
        c(0,0,0),  # cp1 unchanged
        c(0,1,0),  # cp2: +1 on mode 2
        c(2,0,0),  # cp3: +2 on mode 1
        c(0,0,0)   # cp4 unchanged
      )
      add_factors_coeff <- list(
        c(0,  0,  0),  # cp1 unchanged
        c(0,  0.3, 0), # cp2: AR 0.3 for new mode-2 factors
        c(0.3,0,  0),  # cp3: AR 0.3 for new mode-1 factors
        c(0,  0,  0)   # cp4 unchanged
      )
      shift_ind  <- replicate(4, list(NULL, NULL, NULL), simplify = FALSE)
      shift_mean <- c(0,0,0,0)
      shift_var  <- c(0,0,0,0)
      true_cp <- theta
    } else if (data_setting == "s2") {
      r0 <- 3
      C0 <- diag(r0)
      C01 <- 5 * diag(r0)
      C1 <- matrix(rnorm(r0^2, sd = 1/sqrt(r0)), nrow = r0); C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
      C2 <- matrix(c(1,0,0, 0,1,0, 0,0,0), 3, 3, byrow = TRUE)
      C3 <- matrix(0, 3, 3)
      C3[1,1] <- 0.5; C3[2,1] <- rnorm(1); C3[2,2] <- 1
      C3[3,1] <- rnorm(1); C3[3,2] <- rnorm(1); C3[3,3] <- 1.5
      
      transform_list <- list(
        list(C0,  diag(3), diag(3)),
        list(C01, diag(3), diag(3)),
        list(C1,  diag(3), diag(3)),
        list(diag(3), C3,  diag(3)),
        list(diag(3), diag(3), C2)
      )
      type_mode   <- list(1,1,1,2,3)
      type_change <- list("l","l","l","l","l")
      shift_ind <- shift_mean <- shift_var <- NULL
      add_factors <- add_factors_coeff <- NULL
      true_cp <- theta
      
    } else if (data_setting == "s3") {
      type_mode <- list(
        c(1),        # cp1: no change
        c(2, 3),     # cp2: f on mode 2, l on mode 3
        c(1),        # cp3: f on mode 1
        c(1)         # cp4: no change
      )
      type_change <- list(
        character(0),        # cp1: no change
        c("f","l"),          # cp2
        "f",                 # cp3
        character(0)         # cp4: no change
      )
      transform_list <- list(
        list(make_I,   NULL,     NULL),  # cp1: I on mode 1
        list(NULL,     NULL,     make_3I), # cp2: 3I on mode 3
        NULL,                               # cp3: f only
        list(NULL,     NULL,     make_I)      # cp4: I on mode 3
      )
      add_factors <- list(
        c(0, 0, 0),    # cp1
        c(0, 1, 0),    # cp2: +1 on mode 2
        c(2, 0, 0),    # cp3: +2 on mode 1
        c(0, 0, 0)     # cp4
      )
      add_factors_coeff <- list(
        c(0.0, 0.0, 0.0),  # cp1
        c(0.0, 0.3, 0.0),  # cp2: AR 0.3 for new mode-2 factors
        c(0.3, 0.0, 0.0),  # cp3: AR 0.3 for new mode-1 factors
        c(0.0, 0.0, 0.0)   # cp4
      )
      true_cp <- theta 
    }
    
    data_sim <- dgp_general(
      model = "tensor",
      dim_obs = dim_obs,
      dim_latent  = dim_latent,
      Time = Time,
      dep = dep,
      coeff = 0.7,
      true_cp = true_cp,
      type_mode = type_mode,
      type_change = type_change,
      shift_ind = shift_ind,
      shift_mean = shift_mean,
      shift_var = shift_var,
      transform_list = transform_list,
      add_factors = add_factors,
      add_factors_coeff = add_factors_coeff,
      dist = dist
    )
    
    X <- data_sim$X
    dim_X <- dim(X)[1:3]
    Ttot  <- dim(X)[4L]
    
    if (isTRUE(reverse)) {
      X <- X[,,, Ttot:1, drop = FALSE]
    }
    
    if (cp_source == "true") {
      cp_vec <- true_cp %||% integer(0)
      if (length(cp_vec) && isTRUE(reverse)) {
        cp_vec <- sort(Ttot - cp_vec)
      }
    } else {
      est <- if (is.null(r_hat)) {
        global_pca(X = X, dim_X = dim_X, centre = TRUE, proj = TRUE)
      } else {
        global_pca(X = X, dim_X = dim_X, r_hat = r_hat, centre = TRUE, proj = TRUE)
      }
      G     <- est$G_proj
      G_dim <- as.vector(est$r_hat)
      dr <- sum(G_dim * (G_dim + 1L) / 2L)
      int_len <- max(2L, round(6 * log(Ttot)))
      thd <- pmax(
        exp(threshold_coef[1] * log(log(Ttot / int_len)) + threshold_coef[2] * log(dr)),
        threshold_coef[3] * log(Ttot)
      )
      out <- TNotSBS(G, G_dim, method = "fixed", threshold = thd, V.diag = V.diag, lrv = lrv)
      cp_vec <- out$est.cp %||% integer(0)
    }
    
    if (!length(cp_vec)) {
      res_list[[sim]] <- list(
        table = data.frame(),
        zeta_mat = matrix(numeric(0), nrow = 0, ncol = 0),
        cp_vec = integer(0),
        info = list(data_setting = data_setting, cp_source = cp_source, reversed = reverse)
      )
      next
    }
    
    zeta_res <- m11_zeta_from_X(
      X = X,
      dim_obs = dim_X,
      r_hat = r_hat,
      cp_vec = cp_vec,
      st = 1L, ed = Ttot,
      centre_cov = FALSE,
      xi_norm = xi_norm,
      scalar_type = scalar_type,
      extra_div = extra_div,
      return_Xi = FALSE,
      centre = TRUE,
      G_hac = G_hac
    )
    
    res_list[[sim]] <- c(zeta_res,
                         list(info = list(data_setting = data_setting,
                                          cp_source = cp_source,
                                          seed = seed_start + sim,
                                          reversed = reverse)))
  }
  
  # align - for est not true
  same_q <- all(vapply(res_list, function(x) ncol(x$zeta_mat), integer(1L)) ==
                  ncol(res_list[[which.max(vapply(res_list, function(x) ncol(x$zeta_mat), integer(1L)) )]]$zeta_mat))
  agg <- NULL
  if (same_q && length(res_list) && ncol(res_list[[1]]$zeta_mat) > 0) {
    K <- nrow(res_list[[1]]$zeta_mat)
    q <- ncol(res_list[[1]]$zeta_mat)
    agg <- array(NA_real_, dim = c(K, q, nrep),
                 dimnames = list(rownames(res_list[[1]]$zeta_mat),
                                 colnames(res_list[[1]]$zeta_mat),
                                 paste0("rep", seq_len(nrep))))
    for (i in seq_len(nrep)) {
      zm <- res_list[[i]]$zeta_mat
      if (all(dim(zm) == c(K, q))) agg[,,i] <- zm
    }
  }
  
  list(reps = res_list, zeta_array = agg)
}




library(ggplot2)

plot_zeta_array <- function(zeta_array, bins = 20) {
  stopifnot(length(dim(zeta_array)) == 3)
  
  K <- dim(zeta_array)[1]; q <- dim(zeta_array)[2]; nrep <- dim(zeta_array)[3]
  dn <- dimnames(zeta_array)
  
  if (is.null(dn)) {
    dimnames(zeta_array) <- list(
      paste0("mode", seq_len(K)),
      paste0("cp", seq_len(q)),
      paste0("rep", seq_len(nrep))
    )
    dn <- dimnames(zeta_array)
  }
  
  df <- as.data.frame.table(zeta_array, responseName = "xi", stringsAsFactors = FALSE)
  names(df) <- c("Mode", "CP", "Rep", "xi")
  
  df$Mode <- factor(df$Mode, levels = dn[[1]])
  df$CP <- factor(df$CP, levels = dn[[2]])
  
  df <- df[is.finite(df$xi), , drop = FALSE]
  
  ggplot(df, aes(x = xi)) +
    geom_histogram(bins = bins, colour = "grey20") +
    facet_grid(rows = vars(Mode), cols = vars(CP), scales = "free_y") +
    labs(x = expression(xi), y = "Count") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90", colour = NA),
          panel.grid.minor = element_blank())
}
