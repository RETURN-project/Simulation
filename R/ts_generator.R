#' Piecewise linear decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#'
#' @return The time series
#' @export
piecewise <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1) {
  m <- -pert / (2 * thalf) # Slope of the transitory regime
  ttrans <- 2*thalf # Duration of the transitory regime
  y <- offset                             * (t < tpert) +
       (offset + pert + m *(t - tpert))   * (t >= tpert) * (t <= tpert + ttrans) + # Transitory regime
       offset                             * (t > tpert + ttrans)
  y
}


#' Exponential decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#'
#' @return The time series
#' @export
exponential <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1) {
  c <- log(2)/thalf # Translate the half-life to a multiplicative constant
  y <- offset + pert * exp(-c*(t-tpert)) * (t >= tpert)
  y
}
