context("Characterise time series")

test_that("Decomposition harmonic + trend + noise", {
  source('../../R/fun_simulate.R')

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
  source('../../R/fun_simulate.R')
  # Generate time series with predefined ARMA coefficients
  nobs <- 50000
  set.seed(200);innov <- rnorm(nobs) # fix the innovations, so code is reproducible
  set.seed(200);start.innov <- rnorm(4) # fix the innovations of the burn-in period, so code is reproducible

  rn <- arima.sim(model = list(order = c(0, 0, 0) ), n = nobs, innov =innov,
                  n.start = 4, start.innov =  start.innov)
  ma1 <- arima.sim(model = list(order = c(0, 0, 1), ma = .9 ), n = nobs, innov =innov,
                   n.start = 4, start.innov =  start.innov)
  ar1 <- arima.sim(model = list(order = c(1, 0, 0), ar = .9 ), n = nobs, innov =innov,
                   n.start = 4, start.innov =  start.innov)
  ar2 <- arima.sim(model = list(order = c(2, 0, 0), ar = c(.7, .2) ), n = nobs, innov =innov,
                   n.start = 4, start.innov =  start.innov)

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

test_that("Disturbance simulation",{
  source('../../R/fun_simulate.R')
  # disturbance and recovery before end of time series
  distT1 <- 5
  distRec1 <- 12
  distMag1 <- -10
  nobs1 <- 20

  ts1 <- simulDist(distT1, distRec1, distMag1, nobs1)

  # recovery after end of time series
  distT2 <- 5
  distRec2 <- 20
  distMag2 <- -10
  nobs2 <- 20

  ts2 <- simulDist(distT2, distRec2, distMag2, nobs2)

  # disturbance after end of time series
  distT3 <- 25
  distRec3 <- 16
  distMag3 <- -10
  nobs3 <- 20

  ts3 <- simulDist(distT3, distRec3, distMag3, nobs3)

  # tests
  expect_equal(distMag1,min(ts1))# magnitude disturbance
  expect_equal(distT1,which(ts1 == min(ts1)))# timing disturbance
  expect_equal(distRec1,length(which(ts1 < 0)), tolerance = 1)# recovery period

  expect_equal(distMag2,min(ts2))# magnitude disturbance
  expect_equal(distT2,which(ts2 == min(ts2)))# timing disturbance
  expect_equal(distMag2/(distRec2-1), (ts2[distT2] - ts2[distT2 + 1]))# recovery slope

  expect_equal(nobs3,length(which(ts3 == 0)))# no disturbance and recovery
})

test_that("Time series simulation",{
  source('../../R/fun_simulate.R')
  nyr <- 5 # number of years
  nobsyr <- 12 # number of observations per year
  tMiss <- c(1,5,11,23)# observations having missing values
  nDr <- 0 # nimber of drought years
  seasAv <- rep((c(1,2,3,4,5,6,6,5,4,3,2,1)-3.5),3)# seasonal pattern
  seasAmp <- 3 # seasonal amplitude
  trAv <- 5 # offset
  remSd <- 2 # standard deviation remainder
  remMod <- list(order1 = 0, order2 = 0, order3 = 0) # model remainder series
  distMag <- -6 # disturbance magnitude
  distT <- 27 # disturbance timing (observation number)
  distRec <- 25 # recovery period (number of observations)

  # simulate time series
  tsi <- simulTS(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remMod,
                 distMag, distT, distRec)

  expect_equal(tMiss,which(is.na(tsi[[2]][,5])))# missing values
  expect_equal((nyr*nobsyr),length(tsi[[2]][,5]))# number of observations
  expect_equal((seasAv[1:nobsyr])/(max(seasAv[1:nobsyr]))*seasAmp,tsi[[2]][1:nobsyr,1])# seasonal pattern
  expect_equal(seasAmp,max(tsi[[2]][,1]))# seasonal amplitude
  expect_equal(remSd,sd(tsi[[2]][,3]))# standard deviation remainder
  expect_equal(trAv,as.numeric(tsi[[2]][1,2]))# offset
  expect_equal(distMag,min(tsi[[2]][,4]))# magnitude disturbance
  expect_equal(distT,which(tsi[[2]][,4] == min(tsi[[2]][,4])))# timing disturbance
  expect_equal(distRec,length(which(tsi[[2]][,4] < 0)), tolerance = 1)# recovery period
})

test_that('simulation case', {
  source('../../R/fun_simulate.R')
  # simulation settings
  nrep <- 5 # number of repetitions
  nyr <- 5 # number of years
  nobsYr <- 12 # number of observations per year
  mval <- 0.2 # fraction of missing values
  mvaldist <- 'random' # random distribution of missing values
  nDr <- 0 # number of drought years
  seasAv <- rep((c(1,2,3,4,5,6,6,5,4,3,2,1)-3.5),3)# seasonal pattern
  seasAmp <- 3 # seasonal amplitude
  trAv <- 5 # offset
  remSd <- 2 # standard deviation remainder
  remcoef <- list(list(order1 = 0, order2 = 0, order3 = 0),list(order1 = 0, order2 = 0, order3 = 0)) # model remainder series
  distMaglim <- c(-6,-6) # disturbance magnitude
  distTy <- 1 # disturbance timing (year)
  distReclim <- c(25,25) # recovery period (number of observations)

  # simulate set of time series
  tsi <- simulCase(nrep, nyr, nobsYr, nDr,seasAv,seasAmp,
                   trAv, remSd,distMaglim,distTy,distReclim, remcoef, mval, mvaldist)

  # tests
  expect_equal(nrep,dim(tsi$timeSeries)[1])# number of repetitions
  expect_equal((nyr*nobsYr),dim(tsi$timeSeries)[2])# number of repetitions
  expect_equal(apply(tsi$Disturbance, 1, min), rep(distMaglim[1],nrep))# disturbance magnitude
  expect_equal(apply(tsi$Remainder,1,sd), rep(remSd,nrep))# standard deviation of remainder
  expect_equal(colSums(apply(tsi$timeSeries,1,is.na))/(nyr*nobsYr), rep(mval,nrep))# fraction of missing values
})


