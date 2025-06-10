#' Generate Seeded Intervals
#'
#' Constructs a collection of overlapping intervals of varying lengths across a time series,
#' used to scan for potential change points. The intervals decay in size geometrically,
#' ensuring coverage at multiple resolutions.
#'
#' @param n An integer specifying the total length of the time series. Must be at least 3.
#' @param a A decay parameter controlling the geometric scaling of intervals. Must be in \eqn{[1/2, 1)}.
#'          Default is \eqn{(1/2)^1 = 0.5}.
#' @param minl An integer specifying the minimum allowed length of any interval. Default is 2.
#'
#' @return A data frame with columns \code{st} (start index) and \code{ed} (end index) for each interval.
#'         The intervals are sorted in increasing order of \code{st} and \code{ed}, and duplicates are removed.
#'
#' @details
#' The intervals are generated at multiple resolution levels indexed by \eqn{k}, with resolution increasing
#' (i.e., interval size decreasing) as \eqn{k} increases. For each level, \eqn{2m - 1} intervals are placed,
#' where \eqn{m = \lceil (1/a)^{k - 1} \rceil}. The width of intervals is controlled by the decay rate \code{a}.
#'
#' The first interval always spans the full series \eqn{[0, n]}.
#'
#' @examples
#' # Generate seeded intervals for a series of length 100
#' intervals <- seeded_intervals(n = 100)
#' head(intervals)
#'
#' @export

seeded_intervals <- function(n, a = (1/2)^(1), minl = 2) {
  if (n < 3) stop("n should be at least 3")
  if (a < 0.5 || a >= 1) stop("Decay parameter 'a' must be in [1/2, 1)")

  Time <- n
  intervals <- data.frame(st = 0, ed = n)

  k_max <- ceiling(log(Time) / log(1 / a))

  for (k in 2:k_max) {
    m <- ceiling((1 / a)^(k - 1))
    if (m <= 1) next

    D <- Time - Time * a^(k - 1)
    denom <- 2 * m - 2

    for (i in 1:(2 * m - 1)) {
      st <- floor((i - 1) * D / denom)
      ed <- ceiling((i + 1) * D / denom + Time * a^(k - 1))
      ed <- min(Time, ed)

      if ((ed - st) >= minl) {
        intervals <- rbind(intervals, data.frame(st = st, ed = ed))
      }
    }
  }

  intervals <- unique(intervals)
  intervals <- intervals[order(intervals$st, intervals$ed), ]

  return(intervals)
}
