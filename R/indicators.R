#' YrYr recovery function
#'
#' @param ts Vector containing the times (same size as ys)
#' @param ys Vector containing the values (same size as ts)
#' @param tpert Time of the perturbation. Default set to 0 yr
#' @param deltat Reference time (in years). Default set to 5 yr
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

  # The result is the mean slope between t = 0 and t = deltat
  return( (V(tpert + deltat) - V(tpert)) / deltat )
}
