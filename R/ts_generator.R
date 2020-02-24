#' Piecewise linear decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param recover Perturbation half-life
#'
#' @return The time series
#' @export
#'
#' @examples
piecewise <- function(t, offset = 0, pert = 0, tpert = 0, recover = 1) {
  halflife <- (recover - tpert) / 2
  m <- -pert/(recover - tpert) # Slope of the transitory regime
  y <- offset                             * (t <= tpert) +
       (offset + pert + m *(t - tpert))   * (t > tpert) * (t <= recover) + # Transitory regime
       offset                             * (t > recover)
  y
}


#' Exponential decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param recover Perturbation half-life
#'
#' @return The time series
#' @export
#'
#' @examples
exponential <- function(t, offset = 0, pert = 0, tpert = 0, recover = 1) {
  y <- offset + pert * exp(-recover*(t-tpert)) * (t >= tpert)
  y
}

