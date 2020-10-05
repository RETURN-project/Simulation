#' YrYr recovery function
#'
#' Year on year average
#'
#' @param ts Vector containing the times (same size as ys)
#' @param ys Vector containing the values (same size as ts)
#' @param tpert Time of the perturbation. Default set to 0 yr
#' @param deltat Reference recovery time (in years). Default set to 5 yr
#'
#' @return The YrYr parameter for the given time series
#' @export
#'
#' @references Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018).
#' Analyzing spatial and temporal variability in short-term rates of post-fire vegetation return
#' from Landsat time series. Remote Sensing of Environment, 205, 32–45.
#' \url{https://doi.org/10.1016/j.rse.2017.11.007}
#'
#' @examples
#' # Generate an example time series
#' ts <- seq(-2, 10, by = 0.1) # as a vector of times
#' ys <- 3 + 2 * ts # plus a vector of values
#' yryr(ts, ys)
yryr <- function(ts, ys, tpert=0, deltat=5) {
  # Check input
  if ((min(tpert) < min(ts, na.rm = TRUE)) || (min(tpert) + max(deltat) > max(ts, na.rm = TRUE))) {
  # if ((tpert < min(ts, na.rm = TRUE)) || (tpert + deltat > max(ts, na.rm = TRUE))) {
    stop("Error: 'tpert' and/or 'tpert + deltat' are outside the bounds imposed by 'ts'")
  }

  # Auxiliary interpolation function. Given a time, returns the corresponding value.
  # If the time is in ts, returns the corresponding ys. If the time is not in ts,
  # returns a linearly interpolated value
  V <- approxfun(x = ts, y = ys)

  dy <- (mean(V(min(tpert, na.rm = T) + deltat)) - mean(V(tpert)))
  dx <- (mean(min(tpert, na.rm = T) + deltat) - mean(tpert))

  # The result is the mean slope between t = tpert and t = tpert + deltat
  return(dy / dx)
}

#' R80p recovery function
#'
#' Ratio of Eighty Percent (R80P). The ratios can be customized.
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
#' @references Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018).
#' Analyzing spatial and temporal variability in short-term rates of post-fire vegetation return
#' from Landsat time series. Remote Sensing of Environment, 205, 32–45.
#' \url{https://doi.org/10.1016/j.rse.2017.11.007}
#'
#' @examples
#' # Generate an example time series
#' ts <- seq(-2, 10, by = 0.1) # as a vector of times
#' ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.25) # plus a vector of values
#' r80p(ts, ys)
r80p <- function(ts, ys, r=0.8, ts_pre=c(-1,-2), ts_post=c(4, 5)) {
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

  # Return the result
  return(Vpost / (Vpre * r))
}

#' RRI recovery function
#'
#' Relative Recovery Indicator
#'
#' @param ts Vector containing the times (same size as ys). Perturbation assumed to happen at 0
#' @param ys Vector containing the values (same size as ts)
#' @param ts_pre Sampling times for estimating Vpre. Default set to -1
#' @param ts_post Sampling times for estimating Vpost. Default set to c(4, 5)
#' @param tpert Time of the perturbation. Default set to 0 yr
#'
#' @return The RRI indicator
#' @export
#'
#' @references Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018).
#' Analyzing spatial and temporal variability in short-term rates of post-fire vegetation return
#' from Landsat time series. Remote Sensing of Environment, 205, 32–45.
#' \url{https://doi.org/10.1016/j.rse.2017.11.007}
#'
#' @examples
#' # Generate an example time series
#' ts <- seq(-2, 10, by = 0.1) # as a vector of times
#' ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.25) # plus a vector of values
#' rri(ts, ys)
rri <- function(ts, ys, tpert=0, ts_pre=-1, ts_post=c(4, 5)) {
  # Auxiliary interpolation function. Given a time, returns the corresponding value.
  # If the time is in ts, returns the corresponding ys. If the time is not in ts,
  # returns a linearly interpolated value
  V <- approxfun(x = ts, y = ys)

  # The disturbance is assumed to happen at t = 0
  Vdist <- mean(V(tpert))

  # The typical value before perturbation is defined as the average of the values
  # sampled at the times contained in ts_pre
  Vpre  <- mean(V(ts_pre))

  # The typical value after perturbation is defined as the maximum of the values
  # sampled at the times contained in ts_post
  Vpost <- max(V(ts_post))

  # The deltas are closely related to the typical values
  delta_dist <- abs(Vdist - Vpre)
  delta_rec <- abs(Vpost - Vdist)

  return(delta_rec / delta_dist)
}
