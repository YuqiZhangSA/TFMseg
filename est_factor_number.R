#' Estimate the Number of Factors (Barigozzi, Cho and Trapani (2024))
#'
#' Approach proposed by Alessi et al.(2010) in combination with the three information criteria of Bai and Ng (2002) procedure
#' to estimate the number of factors in high-dimensional approximate factor models.
#' This method is designed to be stable across a range of penalty constants and is suitable for both stationary and change point models.
#'
#' @param x A matrix of dimension \code{p × n}, where \code{p} is the number of variables and \code{n} the number of time points.
#' @param r.max Maximum number of factors to consider. If \code{NULL}, defaults to \code{min(50, floor(sqrt(min(n - 1, p)))}.
#' @param center Logical. Should the columns of \code{x} be mean-centred? Default is \code{TRUE}.
#' @param p.seq A numeric vector of row (dimension) subsample sizes used for computing the cross-validation error. If \code{NULL}, default values are used.
#' @param n.seq A numeric vector of time subsample sizes used for computing the cross-validation error. If \code{NULL}, default values are used.
#' @param do.plot Logical. If \code{TRUE}, produces diagnostic plots of estimated rank and variability across constants. Default is \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{r.hat}{A numeric vector of estimated ranks (length 6), corresponding to six different IC variants.}
#' }
#' With attributes:
#' \describe{
#'   \item{Sc}{Matrix of variances of selected ranks across constants.}
#'   \item{const.seq}{The sequence of penalty constants used.}
#'   \item{r.mat}{Array of selected ranks across constants and subsamples.}
#' }
#'
#' @details
#' This function implements the ABC criterion for factor number estimation as proposed in Barigozzi, Cho, and Trapani (2024).
#' It uses a set of IC-type (information criterion) rules indexed by penalty constants and subsample sizes, selecting
#' the rank with the smallest variability across constants. The function evaluates both standard and log-transformed ICs.
#'
#' The output \code{r.hat} contains six values corresponding to the combinations of three penalties (IC1, IC2, IC3)
#' and two transformation types (raw and log).
#'
#' @references
#' Barigozzi, M., Cho, H., & Trapani, L. (2024). *Moving sum procedure for multiple change point detection in large factor models.*
#' arXiv preprint arXiv:2410.02918.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 200), nrow = 100)
#' result <- abc.factor.number(X)
#' result$r.hat
#' }

abc.factor.number <- function(x, r.max = NULL, center = TRUE,
                              p.seq = NULL, n.seq = NULL, do.plot = FALSE) {

  p <- dim(x)[1]
  n <- dim(x)[2]
  ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
  xx <- x - mean.x

  if(is.null(r.max)) r.max <- min(50, floor(sqrt(min(n - 1, p))))

  if(is.null(p.seq)) p.seq <- floor(4 * p / 5 + (1:10) * p / 50)
  if(is.null(n.seq)) n.seq <- floor(4 * n / 5 + (1:10) * n / 50)
  const.seq <- seq(.01, 3, by = 0.01)
  IC <- array(Inf, dim = c(r.max + 1, length(const.seq), 10, 6))

  for(kk in 1:min(length(n.seq), length(p.seq))) {
    nn <- n.seq[kk]
    pp <- p.seq[kk]

    int <- sort(sample(n, nn, replace = FALSE))

    pen <- c((nn + pp) / (nn * pp) * log(nn * pp / (nn + pp)),
             (nn + pp) / (nn * pp) * log(min(nn, pp)),
             log(min(nn, pp)) / min(nn, pp))

    covx <- xx[, int] %*% t(xx[, int]) / nn

    sv <- svd(covx[1:pp, 1:pp], nu = 0, nv = 0)
    tmp <- rev(cumsum(rev(sv$d))) / pp
    if(pp > r.max) tmp <- tmp[1:(r.max + 1)]
    for (jj in 1:length(const.seq)) {
      for (ic.op in 1:3) {
        IC[1:length(tmp), jj, kk, ic.op] <-
          tmp + (1:length(tmp) - 1) * const.seq[jj] * pen[ic.op]
        IC[1:length(tmp), jj, kk, 3 * 1 + ic.op] <-
          log(tmp) + (1:length(tmp) - 1) * const.seq[jj] * pen[ic.op]
      }
    }
  }

  r.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(r.mat, c(1, 3), var)
  r.hat <- rep(0, 6)
  for(ii in 1:6){
    ss <- Sc[, ii]
    if(min(ss) > 0){
      r.hat[ii] <- min(r.mat[max(which(ss == min(ss))),, ii]) - 1
    } else{
      if(sum(ss[-length(const.seq)] != 0 & ss[-1] == 0)) {
        r.hat[ii] <-
          r.mat[which(ss[-length(const.seq)] != 0 &
                        ss[-1] == 0)[1] + 1, dim(r.mat)[2], ii] - 1
      }else{
        r.hat[ii] <- min(r.mat[max(which(ss == 0)),, ii]) - 1
      }
    }
  }

  out <- list(r.hat = r.hat)
  attr(out, "data") <- list(Sc = Sc, const.seq = const.seq, r.mat = r.mat)

  if(do.plot){
    data <- attr(out, "data")
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow = c(2, 3))
    Sc <- data$Sc
    const.seq <- data$const.seq
    q.mat <- data$q.mat
    for(ii in 1:6){
      plot(const.seq, r.mat[, dim(r.mat)[2], ii] - 1, type = "b", pch = 1, col = 2,
           bty = "n", axes = FALSE, xlab = "constant", ylab = "", main = paste("IC ", ii))
      box()
      axis(1, at = pretty(range(const.seq)))
      axis(2, at = pretty(range(r.mat[, dim(r.mat)[2], ii] - 1)),
           col = 2, col.ticks = 2, col.axis = 2)
      par(new = TRUE)
      plot(const.seq, Sc[, ii], col = 4, pch = 2, type = "b",
           bty = "n", axes = FALSE, xlab = "", ylab = "")
      axis(4, at = pretty(range(Sc[, ii])), col = 4, col.ticks = 4, col.axis = 4)
      legend("topright", legend = c("r", "Sc"), col = c(2, 4), lty = c(1, 1), pch = c(1, 2), bty = "n")
    }
  }

  return(out)

}
