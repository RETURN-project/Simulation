#' YrYr recovery function
#'
#' @param ts Vector containing the times (same size as ys)
#' @param ys Vector containing the values (same size as ts)
#' @param tpert Time of the perturbation. Default set to 0 yr
#' @param deltat Reference recovery time (in years). Default set to 5 yr
#'
#' @return The YrYr parameter for the given time series
#' @export
#'
#' @examples
#' # Generate an example time series
#' ts <- seq(-2, 10, by = 0.1) # as a vector of times
#' ys <- 3 + 2 * ts # plus a vector of values
#' yryr(ts, ys)
yryr <- function(ts, ys, tpert=0, deltat=5) {
  # Check input
  if ((tpert < min(ts)) || (tpert + deltat > max(ts))) {
    stop("Error: 'tpert' and/or 'tpert + deltat' are outside the bounds imposed by 'ts'")
  }

  # Auxiliary interpolation function. Given a time, returns the corresponding value.
  # If the time is in ts, returns the corresponding ys. If the time is not in ts,
  # returns a linearly interpolated value
  V <- approxfun(x = ts, y = ys)

  # The result is the mean slope between t = tpert and t = tpert + deltat
  return((V(tpert + deltat) - V(tpert)) / deltat)
}

#' R80p recovery function
#'
#' @param ts Vector containing the times (same size as ys). Perturbation assumed to happen at 0
#' @param ys Vector containing the values (same size as ts)
#' @param r Ratio. Default set to 0.8
#' @param ts_pre Sampling times for estimating Vpre. Default set to -1
#' @param ts_post Sampling times for estimating Vpost. Default set to c(4, 5)
#'
#' @return The R80p indicator (R_r_p if r != 0.8 is provided)
#' @export
#'
#' @examples
#' # Generate an example time series
#' ts <- seq(-2, 10, by = 0.1) # as a vector of times
#' ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.25) # plus a vector of values
#' r80p(ts, ys)
r80p <- function(ts, ys, r=0.8, ts_pre=-1, ts_post=c(4, 5)) {
  # Auxiliary interpolation function. Given a time, returns the corresponding value.
  # If the time is in ts, returns the corresponding ys. If the time is not in ts,
  # returns a linearly interpolated value
  V <- approxfun(x = ts, y = ys)

  # The typical value before perturbation is defined as the average of the values
  # sampled at the times contained in ts_pre
  Vpre  <- mean(V(ts_pre))

  # The typical value after perturbation is defined as the maximum of the values
  # sampled at the times contained in ts_post
  Vpost <- max(V(ts_post))

  return(Vpost / (Vpre * r))
}
