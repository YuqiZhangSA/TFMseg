#' General data‐generating function for vector, matrix and tensor factor models with multiple change points
#'
#' @param model Character; one of "vector", "matrix", or "tensor".
#' @param dim_obs Integer vector (length = model order) of observation dimensions per mode.
#' @param dim_latent Integer vector (length = model order) of initial factor dimensions per mode.
#' @param Time Integer; total number of time points.
#' @param dep Logical; if FALSE, factors are i.i.d. over time (sets coeff = 0).
#' @param coeff Numeric in (-1,1); AR(1) coefficient for all initial factors.
#' @param true_cp Integer vector of sorted change points (< Time).
#' @param type_mode List length q_true; at cp j, integer vector of modes that change.
#' @param type_change List length q_true; each contains any of "l" (loading shift/transform) and/or "f" (factor-number change).
#' @param shift_ind Optional list length q_true; each a list length = M; per change & mode: c(n_rows, n_cols) submatrix to shift; NULL = skip.
#' @param shift_mean Numeric vector length q_true; mean for Gaussian loading perturbations.
#' @param shift_var Numeric vector length q_true; variance for loading perturbations.
#' @param transform_list Optional list length q_true; each a list length = M; matrices/functions to apply to loadings when "l" present.
#' @param add_factors Optional list length q_true; each an integer vector length = M of new factor counts per mode when "f" present.
#' @param add_factors_coeff Optional list length q_true; numeric vector length = M of AR(1) coefficients for newly added factors.
#' @param idio Logical; include idiosyncratic term E_t.
#' @param dist Character; "Gaussian" or "heavy" (t_7 rescaled to var 1).
#'
#' @return list(X, Ft_hist, loadings_hist)

dgp_general <- function(
    model = c("vector", "matrix", "tensor"),
    dim_obs,
    dim_latent,
    Time,
    dep = TRUE,
    coeff,
    true_cp,
    type_mode,
    type_change,
    shift_ind = NULL,
    shift_mean = NULL,
    shift_var = NULL,
    transform_list = NULL,
    add_factors = NULL,
    add_factors_coeff = NULL,
    idio = TRUE,
    dist = c("Gaussian","heavy")
) {
  model <- match.arg(model)
  dist <- match.arg(dist)
  M <- switch(model, vector = 1L, matrix = 2L, tensor = 3L)
  stopifnot(length(dim_obs) == M, length(dim_latent) == M)
  q_true <- length(true_cp)
  stopifnot(q_true == length(type_mode), q_true == length(type_change))
  cps <- sort(true_cp)
  if (!dep) coeff <- 0
  
  # -------- helpers --------
  runit_noise <- function(n) {
    if (dist == "Gaussian") rnorm(n) else {
      df <- 7
      rt(n, df = df) * sqrt((df - 2) / df)  # unit variance
    }
  }
  matrixvariate_draw <- function(shape, size)
    array(runit_noise(prod(shape) * size), c(shape, size))
  
  VAR1_vec <- function(n, size, coeff) {
    stopifnot(abs(coeff) < 1)
    A <- coeff * diag(n)
    innov_sd <- sqrt(1 - coeff^2)  # stationary var 1
    ft <- runit_noise(n)
    out <- array(0, c(n, size))
    out[, 1] <- ft
    for (t in 2:size) {
      ft <- A %*% ft + innov_sd * runit_noise(n)
      out[, t] <- ft
    }
    out
  }
  
  mode_n_prod <- function(Tarr, Mmat, mode) {
    dims <- dim(Tarr)
    if (mode == 1) {
      pm <- Mmat %*% matrix(Tarr, nrow = dims[1])
      return(array(pm, c(nrow(Mmat), dims[2], dims[3])))
    }
    if (mode == 2) {
      pm <- Mmat %*% matrix(aperm(Tarr, c(2,1,3)), nrow = dims[2])
      tmp <- array(pm, c(nrow(Mmat), dims[1], dims[3]))
      return(aperm(tmp, c(2,1,3)))
    }
    if (mode == 3) {
      pm <- Mmat %*% matrix(aperm(Tarr, c(3,1,2)), nrow = dims[3])
      tmp <- array(pm, c(nrow(Mmat), dims[1], dims[2]))
      return(aperm(tmp, c(2,3,1)))
    }
    stop("mode must be 1/2/3.")
  }
  
  apply_shift <- function(mat, num_idx, mu, sd) {
    if (is.null(num_idx) || length(num_idx) != 2) return(mat)
    rows <- sample(nrow(mat), num_idx[1])
    cols <- sample(ncol(mat), num_idx[2])
    mat[rows, cols] <- mat[rows, cols] + rnorm(prod(num_idx), mu, sd)
    mat
  }
  
  # -------- initial loadings + factors --------
  loadings_current <- vector("list", M)
  for (m in seq_len(M)) {
    loadings_current[[m]] <- matrix(
      runif(dim_obs[m] * dim_latent[m], -1, 1),
      nrow = dim_obs[m]
    )
  }
  cur_latent <- dim_latent
  final_latent <- dim_latent
  if (!is.null(add_factors)) {
    for (j in seq_len(q_true)) final_latent <- final_latent + add_factors[[j]]
  }
  
  if (M == 1) {
    Ft_full <- array(0, c(final_latent, Time))
  } else if (M == 2) {
    Ft_full <- array(0, c(final_latent[1], final_latent[2], Time))
  } else {
    Ft_full <- array(0, c(final_latent, Time))
  }
  
  base <- VAR1_vec(prod(dim_latent), Time, coeff)
  Ft_base <- array(base, c(dim_latent, Time))
  if (M == 1) {
    Ft_full[1:dim_latent, ] <- Ft_base
  } else if (M == 2) {
    Ft_full[1:dim_latent[1], 1:dim_latent[2], ] <- Ft_base
  } else {
    Ft_full[1:dim_latent[1], 1:dim_latent[2], 1:dim_latent[3], ] <- Ft_base
  }
  
  loadings_hist <- vector("list", Time)
  Ft_hist <- vector("list", Time)
  
  E <- switch(model,
              vector = array(matrixvariate_draw(c(dim_obs,1), Time), c(dim_obs,1,Time)),
              matrix = matrixvariate_draw(dim_obs, Time),
              tensor = matrixvariate_draw(dim_obs, Time))
  
  record_state <- function(t) {
    loadings_hist[[t]] <<- lapply(loadings_current, function(L) L)
    if (M == 1) {
      Ft_hist[[t]] <<- Ft_full[1:cur_latent, t]
    } else if (M == 2) {
      Ft_hist[[t]] <<- Ft_full[1:cur_latent[1], 1:cur_latent[2], t]
    } else {
      slice <- Ft_full[1:cur_latent[1], 1:cur_latent[2], 1:cur_latent[3], t]
      Ft_hist[[t]] <<- array(slice, dim = cur_latent)
    }
  }
  
  # -------- iterate segments and apply changes --------
  start_t <- 1L
  for (j in seq_len(q_true + 1L)) {
    end_t <- if (j <= q_true) cps[j] else Time
    for (t in start_t:end_t) record_state(t)
    if (j > q_true) break
    
    modes_j <- type_mode[[j]]
    ops_j   <- type_change[[j]]
    
    # (1) loading shifts / transforms
    if ("l" %in% ops_j) {
      for (m in modes_j) {
        if (!is.null(shift_ind)) {
          loadings_current[[m]] <- apply_shift(
            loadings_current[[m]],
            shift_ind[[j]][[m]],
            shift_mean[j],
            sqrt(shift_var[j])
          )
        }
        if (!is.null(transform_list) && length(transform_list) >= j) {
          Tj <- transform_list[[j]]
          if (!is.null(Tj) && length(Tj) >= m) {
            Tjm <- Tj[[m]]
            if (!is.null(Tjm)) {
              Lm <- loadings_current[[m]]
              if (is.function(Tjm)) Tjm <- Tjm(ncol(Lm))
              Tjm <- as.matrix(Tjm)
              if (!is.numeric(Tjm)) stop("Transform must be numeric.")
              if (ncol(Lm) != nrow(Tjm)) {
                stop(sprintf("Transform dim mismatch at cp %d, mode %d: ncol(L)=%d, nrow(T)=%d",
                             j, m, ncol(Lm), nrow(Tjm)))
              }
              loadings_current[[m]] <- Lm %*% Tjm
            }
          }
        }
      }
    }
    
    # (2) factor-number changes — add new factors
    if ("f" %in% ops_j) {
      for (m in modes_j) {
        k <- add_factors[[j]][m]
        if (k > 0) {
          t0 <- cps[j] + 1L
          t1 <- Time                 
          L  <- t1 - t0 + 1L
          
          prod_other <- if (M == 1) 1L else prod(cur_latent[-m])
          raw <- VAR1_vec(k * prod_other, L, add_factors_coeff[[j]][m])
          
          for (ii in seq_len(L)) {
            idx <- t0 + ii - 1L
            if (M == 2) {
              if (m == 1) {
                arr_n <- array(raw[, ii], c(k, cur_latent[2]))
                Ft_full[(cur_latent[1] + 1):(cur_latent[1] + k), seq_len(cur_latent[2]), idx] <- arr_n
              }
              if (m == 2) {
                arr_n <- array(raw[, ii], c(cur_latent[1], k))
                Ft_full[seq_len(cur_latent[1]), (cur_latent[2] + 1):(cur_latent[2] + k), idx] <- arr_n
              }
            }
            if (M == 3) {
              sel1 <- seq_len(cur_latent[1]); sel2 <- seq_len(cur_latent[2]); sel3 <- seq_len(cur_latent[3])
              if (m == 1) {
                arr_n <- array(raw[, ii], c(k, cur_latent[2], cur_latent[3]))
                Ft_full[(cur_latent[1] + 1):(cur_latent[1] + k), sel2, sel3, idx] <- arr_n
              }
              if (m == 2) {
                arr_n <- array(raw[, ii], c(cur_latent[1], k, cur_latent[3]))
                Ft_full[sel1, (cur_latent[2] + 1):(cur_latent[2] + k), sel3, idx] <- arr_n
              }
              if (m == 3) {
                arr_n <- array(raw[, ii], c(cur_latent[1], cur_latent[2], k))
                Ft_full[sel1, sel2, (cur_latent[3] + 1):(cur_latent[3] + k), idx] <- arr_n
              }
            }
          }
          
          # extend loadings to match the new factors and update rank
          loadings_current[[m]] <- cbind(
            loadings_current[[m]],
            matrix(runif(dim_obs[m] * k, -2, 2), nrow = dim_obs[m])
          )
          cur_latent[m] <- cur_latent[m] + k
        }
      }
    }
    
    start_t <- cps[j] + 1L
  }
  
  # -------- observations --------
  if (idio) {
    X <- switch(model,
                vector = {
                  mat <- matrix(0, dim_obs, Time)
                  for (t in 1:Time) mat[, t] <- loadings_hist[[t]][[1]] %*% Ft_hist[[t]] + E[,1,t]
                  mat
                },
                matrix = {
                  arr <- array(0, c(dim_obs, Time))
                  for (t in 1:Time) arr[,,t] <- loadings_hist[[t]][[1]] %*% Ft_hist[[t]] %*% t(loadings_hist[[t]][[2]]) + E[,,t]
                  arr
                },
                tensor = {
                  arr <- array(0, c(dim_obs, Time))
                  for (t in 1:Time) {
                    A1 <- mode_n_prod(Ft_hist[[t]], loadings_hist[[t]][[1]], 1)
                    A2 <- mode_n_prod(A1, loadings_hist[[t]][[2]], 2)
                    A3 <- mode_n_prod(A2, loadings_hist[[t]][[3]], 3)
                    arr[,,,t] <- A3 + E[,,,t]
                  }
                  arr
                }
    )
  } else {
    X <- switch(model,
                vector = {
                  mat <- matrix(0, dim_obs, Time)
                  for (t in 1:Time) mat[, t] <- loadings_hist[[t]][[1]] %*% Ft_hist[[t]]
                  mat
                },
                matrix = {
                  arr <- array(0, c(dim_obs, Time))
                  for (t in 1:Time) arr[,,t] <- loadings_hist[[t]][[1]] %*% Ft_hist[[t]] %*% t(loadings_hist[[t]][[2]])
                  arr
                },
                tensor = {
                  arr <- array(0, c(dim_obs, Time))
                  for (t in 1:Time) {
                    A1 <- mode_n_prod(Ft_hist[[t]], loadings_hist[[t]][[1]], 1)
                    A2 <- mode_n_prod(A1, loadings_hist[[t]][[2]], 2)
                    A3 <- mode_n_prod(A2, loadings_hist[[t]][[3]], 3)
                    arr[,,,t] <- A3
                  }
                  arr
                }
    )
  }
  
  list(X = X, Ft_hist = Ft_hist, loadings_hist = loadings_hist)
}
