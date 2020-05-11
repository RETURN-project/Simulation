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

test_that("Realistic decay (deterministic)", {

  ts <- seq(0, 25, by = 0.1)
  pert <- 1.5
  tpert <- 0
  thalf <- 1
  offset <- 2

  ys_exp <- exponential(ts, offset = offset, pert = pert, tpert = tpert, thalf = thalf)
  ys_rea <- realistic(ts, offset = offset, pert = pert, tpert = tpert, thalf = thalf)

  # In the absence of noise, realistic and exponential decays should be identical
  expect_equal(ys_rea, ys_exp, tolerance = 1e-4)
})

test_that("Realistic decay (stochastic)", {

  ts <- seq(0, 1000, by = 0.01)
  offset <- 2
  tpert <- 50
  pert <- 5
  thalf <- 1

  sd_expected <- 0.7 # Desired sd for the generated time series
  sigma <- sd_expected * sqrt(2 * log(2) / thalf)
  # The assignation above uses the relationship between infinitesimal (sigma) and asymptotic (measured) standard deviation.
  # i.e: such a sigma creates a time series with sd_expected (provided t_0 = 0 t_end = inf)
  # More info: https://math.stackexchange.com/questions/2558659/expectation-and-variance-of-stochastic-differential-equations

  ys_r <- realistic(ts, offset = offset, pert = pert, tpert = tpert, thalf = thalf, noise = 0, sigma = sigma)
  sd_measured <- sd(ys_r)

  expect_equal(sd_measured, sd_expected, tolerance = 0.1)

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

test_that('Average parameter generation',{
  # generate a settings file
  set.seed(197)
  Vqntl <- c( .05, .275, .5, .725, .95)
  seasVImax <- rnorm(100)
  tsVIMissVal <- rnorm(100)
  Rem_VIsd <- rnorm(100)
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = quantile(seasVImax, Vqntl), obs = seasVImax)
  sttngs$'missVal' <- list(type = 'dist',  vals = quantile(tsVIMissVal, Vqntl), obs = tsVIMissVal)
  sttngs$'remSd' <- list(type = 'dist',  vals = quantile(Rem_VIsd, Vqntl), obs = Rem_VIsd)
  sttngs$'nyr' <- list(type = 'range', vals = c(20,36))
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.1,0.2,0.3,0.4,0.5))
  sttngs$'distT' <- list(type = 'range', vals =10)
  sttngs$'distRec' <- list(type = 'range', vals = seq(0.5,6.5,by=1.5)*6)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0))
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'))
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = 0.7)
  sttngs$'general' <- list(
    eval = c('distMag', 'seasAmp'),
    nTS = 2,
    nobsYr = 6,
    seasAv = c(1,2,3,4,5,6),
    remcoef = list(10,20,30),
    parSetUp = 'avg')#avg dist, comb

  # derive parameters using the settings file
  pars <- setParamValues(sttngs)

  # tests
  expect_equal(length(pars),2)# 2 evaluated parameters
  expect_equal(pars$distMag$`-0.1`$nrep,c(1,1))
  expect_equal(pars$distMag$`-0.1`$nyr,rep(mean(c(20,36)),2))
  expect_equal(pars$distMag$`-0.1`$nobsYr,rep(6,2))
  expect_equal(pars$distMag$`-0.1`$nDr,rep(0,2))
  expect_equal(pars$distMag$`-0.1`$seasAv[[1]],c(1:6))
  expect_equal(pars$distMag$`-0.1`$seasAmp,rep(mean(seasVImax),2), tolerance = 1e-5)
  expect_equal(pars$distMag$`-0.1`$trAv,rep(0.7,2))
  expect_equal(pars$distMag$`-0.1`$remSd,rep(mean(Rem_VIsd),2))
  expect_equal(pars$distMag$`-0.1`$distMag,rep(-0.1,2))
  expect_equal(pars$distMag$`-0.1`$distT,rep(10,2))
  expect_equal(pars$distMag$`-0.1`$missVal,rep(mean(tsVIMissVal),2))
  expect_equal(pars$distMag$`-0.1`$DistMissVal,rep('random',2))
  expect_equal(pars$distMag$`-0.1`$distType,rep('piecewise',2))

})

test_that('Combination parameter generation',{
  # generate a settings file
  set.seed(197)
  Vqntl <- c( .5)
  seasVImax <- rnorm(100)
  tsVIMissVal <- rnorm(100)
  Rem_VIsd <- rnorm(100)
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = quantile(seasVImax, Vqntl), obs = seasVImax)
  sttngs$'missVal' <- list(type = 'dist',  vals = quantile(tsVIMissVal, Vqntl), obs = tsVIMissVal)
  sttngs$'remSd' <- list(type = 'dist',  vals = quantile(Rem_VIsd, Vqntl), obs = Rem_VIsd)
  sttngs$'nyr' <- list(type = 'range', vals = c(20,25))
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.1,0.5))
  sttngs$'distT' <- list(type = 'range', vals =10)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,6.5)*6)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0))
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'))
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = 0.7)
  sttngs$'general' <- list(
    eval = c('distMag', 'seasAmp'),
    nTS = 2,
    nobsYr = 6,
    seasAv = c(1,2,3,4,5,6),
    remcoef = list(10,20,30),
    parSetUp = 'comb')#avg dist, comb

  # derive parameters using the settings file
  pars <- setParamValues(sttngs)

  # tests
  expect_equal(length(pars),2)# 2 evaluated parameters
  expect_equal(pars$distMag$`-0.1`$nrep,rep(1,4))
  expect_equal(pars$distMag$`-0.1`$nyr,rep(c(20,25),2))
  expect_equal(pars$distMag$`-0.1`$nobsYr,rep(6,4))
  expect_equal(pars$distMag$`-0.1`$nDr,rep(0,4))
  expect_equal(pars$distMag$`-0.1`$seasAv[[1]],c(1:6))
  expect_equal(pars$distMag$`-0.1`$seasAmp,rep(quantile(seasVImax, .5),4), tolerance = 1e-5)
  expect_equal(pars$distMag$`-0.1`$trAv,rep(0.7,4))
  expect_equal(pars$distMag$`-0.1`$remSd,rep(quantile(Rem_VIsd,.5),4))
  expect_equal(pars$distMag$`-0.1`$distMag,rep(-0.1,4))
  expect_equal(pars$distMag$`-0.1`$distT,rep(10,4))
  expect_equal(pars$distMag$`-0.1`$missVal,rep(quantile(tsVIMissVal,.5),4))
  expect_equal(pars$distMag$`-0.1`$DistMissVal,as.factor(rep('random',4)))
  expect_equal(pars$distMag$`-0.1`$distType,as.factor(rep('piecewise',4)))

})
