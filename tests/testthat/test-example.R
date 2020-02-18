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


test_that("ARMA coefficients - white noise",{
  source_rmd <- function(file, local = FALSE, ...){
    library('knitr')
    options(knitr.duplicate.label = 'allow')

    tempR <- tempfile(tmpdir = ".", fileext = ".R")
    on.exit(unlink(tempR))
    knitr::purl(file, output=tempR, quiet = TRUE)#

    # source(tempR)
    envir <- globalenv()
    source(tempR, local = envir, ...)
  }

  source_rmd('fun_simulate.Rmd')
  # Generate white noise
  set.seed(197)
  rn <- arima.sim(model = list(order = c(0, 0, 0) ), n = 2000)
  ma1 <- arima.sim(model = list(order = c(0, 0, 1), ma = .9 ), n = 2000)
  ar1 <- arima.sim(model = list(order = c(1, 0, 0), ar = .9 ), n = 2000)
  ma2 <- arima.sim(model = list(order = c(0, 0, 2), ma = c(.7, .2) ), n = 2000)
  getARMAcoef(rn)
  getARMAcoef(ma1)
  getARMAcoef(ar1)
  getARMAcoef(ma2)

})


# ARMA coefficients
getARMAcoef <- function(tsx){

  tsx <- ma2


  arm <- auto.arima(tsx, seasonal=F) #arm$arma #A compact form of the specification, as a vector giving
  #the number of AR, MA, seasonal AR and seasonal MA coefficients,
  #plus the period, and the number of non-seasonal, and seasonal differences.
  teller <- 0
  coefmod <- list()
  # number of AR coef
  if(arm$arma[1]>0){
    coefmod <- c(coefmod, ar = as.numeric(arm$coef[1:arm$arma[1]]))
    teller <- teller + arm$arma[1]
  }
  # number of MA coef
  if(arm$arma[2]>0){
    coefmod <- c(coefmod, ma = as.numeric(arm$coef[(1+teller):(teller+arm$arma[2])]))
    #teller <- teller + arm$arma[2]
  }
  coefmod <- c(coefmod, order = as.numeric(c(arm$arma[1], arm$arma[6], arm$arma[2])))
  coefmod
}



