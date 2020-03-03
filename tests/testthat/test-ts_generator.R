context("Simulation of benchmarking time series")

test_that("Exponential decay", {

  ts <- seq(0, 25, by = 0.1)
  pert <- 1.5
  tpert <- 0
  thalf <- 1

  ys <- exponential(ts, pert = pert, tpert = tpert, thalf = thalf)

  expect_equal(abs(min(ys)), 0, tolerance = 1e-4)
  expect_equal(abs(max(ys)), abs(pert), tolerance = 1e-4)

})

test_that("Linear decay", {

  ts <- seq(0, 5, by = 0.01)
  pert <- 1.5
  ys <- piecewise(ts, pert = pert)

  expect_equal(abs(min(ys)), 0, tolerance = 1e-4)
  expect_equal(abs(max(ys)), abs(pert), tolerance = 1e-4)
})

test_that("Disturbance simulation",{
  # source('../R/fun_simulate.R')

  # disturbance and recovery before end of time series
  distT1 <- 5
  distRec1 <- 12/2
  distMag1 <- -10
  nobs1 <- 20

  ts1 <- piecewise(1:nobs1, pert = distMag1, tpert = distT1, thalf = distRec1)

  # recovery after end of time series
  distT2 <- 5
  distRec2 <- 20/2
  distMag2 <- -10
  nobs2 <- 20

  ts2 <- piecewise(1:nobs2, pert = distMag2, tpert = distT2, thalf = distRec2)

  # disturbance after end of time series
  distT3 <- 25
  distRec3 <- 16/2
  distMag3 <- -10
  nobs3 <- 20

  ts3 <- piecewise(1:nobs3, pert = distMag3, tpert = distT3, thalf = distRec3)

  # tests
  expect_equal(distMag1,min(ts1))# magnitude disturbance
  expect_equal(distT1,which(ts1 == min(ts1)))# timing disturbance
  expect_equal(distRec1,length(which(ts1 < 0)), tolerance = 1)# recovery period

  expect_equal(distMag2,min(ts2))# magnitude disturbance
  expect_equal(distT2,which(ts2 == min(ts2)))# timing disturbance
  expect_equal(distMag2/((2*distRec2)), (ts2[distT2] - ts2[distT2 + 1]))# recovery slope

  expect_equal(nobs3,length(which(ts3 == 0)))# no disturbance and recovery
})

test_that("Time series simulation",{
  #source('../R/fun_simulate.R')

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
  distRec <- 25/2 # recovery half time (number of observations)

  # simulate time series
  tsi <- simulTS(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remMod,
                 distMag, distT, distRec, distType = 'piecewise')

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

test_that('simulation case',{
  #source('../R/fun_simulate.R')
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
  distReclim <- c(25/2,25/2) # recovery half time (number of observations)

  # simulate set of time series
  tsi <- simulCase(nrep, nyr, nobsYr, nDr,seasAv,seasAmp,
                   trAv, remSd,distMaglim,distTy,distReclim, remcoef, mval, mvaldist, distType = 'piecewise')

  # tests
  expect_equal(nrep,dim(tsi$timeSeries)[1])# number of repetitions
  expect_equal((nyr*nobsYr),dim(tsi$timeSeries)[2])# number of repetitions
  expect_equal(apply(tsi$Disturbance, 1, min), rep(distMaglim[1],nrep))# disturbance magnitude
  expect_equal(apply(tsi$Remainder,1,sd), rep(remSd,nrep))# standard deviation of remainder
  expect_equal(colSums(apply(tsi$timeSeries,1,is.na))/(nyr*nobsYr), rep(mval,nrep))# fraction of missing values
})
