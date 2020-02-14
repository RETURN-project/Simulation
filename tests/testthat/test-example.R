context("Characterise time series")

test_that("Decomposition harmonic + trend + noise", {

  source_rmd <- function(file, local = FALSE, ...){
    options(knitr.duplicate.label = 'allow')

    tempR <- tempfile(tmpdir = ".", fileext = ".R")
    on.exit(unlink(tempR))
    knitr::purl(file, output=tempR, quiet = TRUE)

    envir <- globalenv()
    source(tempR, local = envir, ...)
  }

  source_rmd('fun_simulate.Rmd')

  set.seed(197)
  # generate time stamps
  tm <- seq(as.Date('2000-01-01'),as.Date('2005-12-31'), by='1 month',)

  # seasonality represented by harmonic function
  seas <- 0.5*sin(2*pi*(1:(12*6))/12)

  # trend
  tr <- 0.3 + 0.02*(1:(12*6))

  # noise
  rn <- runif(12*6, min=0, max=100)/200
  rn <- rn - mean(rn)

  # time series equals sum of components
  tsi <- as.data.frame(t(c(1,2,t(seas + tr + rn))))
  names(tsi) <- c('lat','lon', tm)

  # decompose
  dc <- decompTSbfast(tsi, 6, 12)

  #compare simulated and derived seasonality
  difftr <- mean(as.numeric(dc$Trend[,-c(1,2)]) - tr)
  diffseas <- mean(as.numeric(dc$Seasonality[,-c(1,2)]) - seas)

  # test
  expect_equal(diffseas, 0, tolerance = 1e-4)
  expect_equal(diffseas, 0, tolerance = 1e-4)
})

