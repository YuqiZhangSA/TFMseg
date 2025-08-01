#' General data‐generating function for vector, matrix and tensor factor models with multiple change points
#'
#' This function simulates observations \eqn{X_t} for a latent factor model of type "vector", "matrix", or "tensor",
#' allowing for multiple change points where loadings and/or factor dimensions may shift.
#'
#' @param model Character; one of "vector", "matrix", or "tensor" indicating factor‐structure type.
#' @param dim_obs Integer vector of length equal to model order (1, 2, or 3) giving observation dimensions per mode.
#' @param dim_latent Integer vector same length as \code{dim_obs}, initial factor dimensions per mode (before any changes).
#' @param Time Integer; total number of time points to simulate.
#' @param coeff Numeric in (-1,1); AR(1) coefficient applied to all initial factors.
#' @param true_cp Integer vector of sorted change points (< Time) where structural changes occur.
#' @param type_mode List of length \eqn{q_true=length(true_cp)}; each element is an integer vector of mode indices that change.
#' @param type_change List of same length as \code{true_cp}; each element is a character vector containing any of "l" (loading shift) and/or "f" (factor‐number change).
#' @param shift_ind Optional list of length \eqn{q_true}, each a list of length=M; for each change and mode: a length‐2 integer vector \code{c(n_rows, n_cols)} specifying submatrix size to shift. NULL skips shift.
#' @param shift_mean Numeric vector length \eqn{q_true}; Gaussian shift mean for loading perturbations at each change.
#' @param shift_var Numeric vector length \eqn{q_true}; variance for loading shifts at each change.
#' @param transform_list Optional list of length \eqn{q_true}, each a list length=M of transformation matrices to apply to loadings when type contains "l".
#' @param add_factors Optional list of length \eqn{q_true}, each an integer vector length=M of new factor counts per mode when type contains "f".
#' @param add_factors_coeff Optional list of length \eqn{q_true}, each numeric vector length=M of AR coefficients for newly added factors.
#'
#' @return A 3‐way array \code{X} (or matrix for the vector model) of simulated observations with dimensions \code{c(dim_obs, Time)}.
#' @export
dgp_general <- function(
    model = c("vector", "matrix", "tensor"),
    dim_obs,
    dim_latent,
    Time,
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
    idio = TRUE
) {
  model <- match.arg(model)
  M <- switch(model, vector = 1, matrix = 2, tensor = 3)
  stopifnot(length(dim_obs) == M, length(dim_latent) == M)
  q_true <- length(true_cp)
  stopifnot(q_true == length(type_mode), q_true == length(type_change))
  cps <- sort(true_cp)
  
  # Preparation
  matrixvariate_normal <- function(shape, size) array(rnorm(prod(shape) * size), c(shape, size))
  tensorvariate_normal <- function(shape, size) array(rnorm(prod(shape) * size), c(shape, size))
  VAR1_vec <- function(n, size, coeff) {
    stopifnot(abs(coeff) < 1)
    A <- coeff * diag(n)
    sigma <- sqrt(1 - coeff^2)
    ft <- rnorm(n)
    out <- array(0, c(n, size))
    out[, 1] <- ft
    for (t in 2:size) {
      ft <- A %*% ft + rnorm(n, 0, sigma)
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
  }
  apply_shift <- function(mat, num_idx, mu, sd) {
    if (is.null(num_idx) || length(num_idx) != 2) return(mat)
    rows <- sample(nrow(mat), num_idx[1])
    cols <- sample(ncol(mat), num_idx[2])
    mat[rows, cols] <- mat[rows, cols] + rnorm(prod(num_idx), mu, sd)
    mat
  }
  
  # Initialisation
  loadings_current <- vector("list", M)
  for (m in seq_len(M)) {
    loadings_current[[m]] <- matrix(runif(dim_obs[m] * dim_latent[m], -1, 1),
                                    nrow = dim_obs[m])
  }
  cur_latent <- dim_latent
  
  final_latent <- dim_latent
  if (!is.null(add_factors)) {
    for (j in seq_len(q_true)) {
      final_latent <- final_latent + add_factors[[j]]
    }
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
  
  # Histories & noise
  loadings_hist <- vector("list", Time)
  Ft_hist <- vector("list", Time)
  # Gaussian noise
  E <- switch(model,
              vector = array(matrixvariate_normal(c(dim_obs,1), Time), c(dim_obs,1,Time)),
              matrix = matrixvariate_normal(dim_obs, Time),
              tensor = tensorvariate_normal(dim_obs, Time)
  )
  
  record_state <- function(t) {
    loadings_hist[[t]] <<- lapply(loadings_current, function(L) L)
    if (M == 1) {
      Ft_hist[[t]] <<- Ft_full[1:cur_latent, t]
    } else if (M == 2) {
      Ft_hist[[t]] <<- Ft_full[1:cur_latent[1], 1:cur_latent[2], t]
    } else {
      Ft_hist[[t]] <<- Ft_full[1:cur_latent[1],
                               1:cur_latent[2],
                               1:cur_latent[3],
                               t]
    }
  }
  
  start_t <- 1
  for (j in seq_len(q_true + 1)) {
    end_t <- if (j <= q_true) cps[j] else Time
    for (t in start_t:end_t) record_state(t)
    if (j > q_true) break
    
    modes_j <- type_mode[[j]]
    ops_j   <- type_change[[j]]
    
    # Loading shifts / transforms
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
        if (!is.null(transform_list)) {
          loadings_current[[m]] <- loadings_current[[m]] %*%
            transform_list[[j]][[m]]
        }
      }
    }
    
    # Factor number changes
    if ("f" %in% ops_j) {
      for (m in modes_j) {
        k <- add_factors[[j]][m]
        if (k > 0) {
          t0 <- cps[j] + 1
          t1 <- if (j < q_true) cps[j+1] else Time
          L  <- t1 - t0 + 1
          raw <- VAR1_vec(k * prod(cur_latent[-m]), L,
                          add_factors_coeff[[j]][m])
          
          for (ii in seq_len(L)) {
            idx <- t0 + ii - 1
            if (M == 2) {
              arr_n <- array(raw[, ii],
                             c(if (m==1) k else cur_latent[1],
                               if (m==2) k else cur_latent[2]))
              if (m==1) Ft_full[(cur_latent[1]+1):(cur_latent[1]+k), , idx] <- arr_n
              if (m==2) Ft_full[, (cur_latent[2]+1):(cur_latent[2]+k), idx] <- arr_n
            }
            if (M == 3) {
              dims <- cur_latent
              arr_n <- array(raw[, ii],
                             c(if (m==1) k else dims[1],
                               if (m==2) k else dims[2],
                               if (m==3) k else dims[3]))
              if (m==1) Ft_full[(dims[1]+1):(dims[1]+k), , , idx] <- arr_n
              if (m==2) Ft_full[, (dims[2]+1):(dims[2]+k), , idx] <- arr_n
              if (m==3) Ft_full[, , (dims[3]+1):(dims[3]+k), idx] <- arr_n
            }
          }
          
          loadings_current[[m]] <- cbind(
            loadings_current[[m]],
            matrix(runif(dim_obs[m] * k, -1, 1), nrow = dim_obs[m])
          )
          cur_latent[m] <- cur_latent[m] + k
        }
      }
    }
    
    start_t <- cps[j] + 1
  }
  
  # Observations
  if (idio) {
    X <- switch(model,
                vector = {
                  mat <- matrix(0, dim_obs, Time)
                  for (t in 1:Time) {
                    mat[, t] <- loadings_hist[[t]][[1]] %*% Ft_hist[[t]] + E[,1,t]
                  }
                  mat
                },
                matrix = {
                  arr <- array(0, c(dim_obs, Time))
                  for (t in 1:Time) {
                    arr[,,t] <- loadings_hist[[t]][[1]] %*%
                      Ft_hist[[t]] %*%
                      t(loadings_hist[[t]][[2]]) +
                      E[,,t]
                  }
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
  } else{
    X <- switch(model,
                vector = {
                  mat <- matrix(0, dim_obs, Time)
                  for (t in 1:Time) {
                    mat[, t] <- loadings_hist[[t]][[1]] %*% Ft_hist[[t]]
                  }
                  mat
                },
                matrix = {
                  arr <- array(0, c(dim_obs, Time))
                  for (t in 1:Time) {
                    arr[,,t] <- loadings_hist[[t]][[1]] %*%
                      Ft_hist[[t]] %*%
                      t(loadings_hist[[t]][[2]])
                  }
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
  
  
  list(
    X = X,
    Ft_hist = Ft_hist,
    loadings_hist = loadings_hist
  )
}
