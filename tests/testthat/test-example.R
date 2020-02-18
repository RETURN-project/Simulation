context("Characterise time series")

test_that("Decomposition harmonic + trend + noise", {
  source('fun_simulate.R')

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

  #compare simulated and derived seasonality and trend
  difftr <- mean(as.numeric(dc$Trend[,-c(1,2)]) - tr)
  diffseas <- mean(as.numeric(dc$Seasonality[,-c(1,2)]) - seas)

  # test
  expect_equal(diffseas, 0, tolerance = 1e-4)
  expect_equal(diffseas, 0, tolerance = 1e-4)
})


test_that("ARMA coefficients",{
  source('fun_simulate.R')
  # Generate time series with predefined ARMA coefficients
  set.seed(197) # does not work, need to look for alternative solution
  rn <- arima.sim(model = list(order = c(0, 0, 0) ), n = 50000)
  ma1 <- arima.sim(model = list(order = c(0, 0, 1), ma = .9 ), n = 50000)
  ar1 <- arima.sim(model = list(order = c(1, 0, 0), ar = .9 ), n = 50000)
  ar2 <- arima.sim(model = list(order = c(2, 0, 0), ar = c(.7, .2) ), n = 50000)

  # estimate coefficients
  coef_rn <- getARMAcoef(rn)
  coef_ma1 <- getARMAcoef(ma1)
  coef_ar1 <- getARMAcoef(ar1)
  coef_ar2 <- getARMAcoef(ar2)

  # test whether coefficients were correnctly estimated
  expect_equal(list(order = c(coef_rn$order1, coef_rn$order2, coef_rn$order3)), list(order = c(0, 0, 0)), tolerance = 1e-4)
  expect_equal(coef_ma1, list(ma = 0.9, order1 = 0, order2 = 0, order3 = 1), tolerance = 1e-1)
  expect_equal(coef_ar1, list(ar = 0.9, order1 = 1, order2 = 0, order3 = 0), tolerance = 1e-1)
  expect_equal(coef_ar2, list(ar1 = 0.7, ar2 = 0.2, order1 = 2, order2 = 0, order3 = 0), tolerance = 1e-1)
})




