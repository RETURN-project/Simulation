#' Piecewise linear decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the additive white noise (standard deviation)
#'
#' @return The time series
#' @export
piecewise <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {
  m <- -pert / (2 * thalf) # Slope of the transitory regime
  ttrans <- 2*thalf # Duration of the transitory regime
  y <- offset                             * (t < tpert) +
       (offset + pert + m *(t - tpert))   * (t >= tpert) * (t <= tpert + ttrans) + # Transitory regime
       offset                             * (t > tpert + ttrans)

  y <- y + rnorm(length(t), sd = noise) # Add the noise

  return(y)
}


#' Exponential decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the additive white noise (standard deviation)
#'
#' @return The time series
#' @export
exponential <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {
  r <- log(2)/thalf # Translate the half-life to a multiplicative constant
  y <- offset + pert * exp(-r*(t-tpert)) * (t >= tpert)

  y <- y + rnorm(length(t), sd = noise) # Add the noise

  return(y)
}


#' Realistic time series simulation
#'
#' Simulates a return to equilibrium after a perturbation under the influence of a stochastic differential equation where:
#'
#' 1. The deterministic dynamics are given by an exponential decay with the given half life
#' 2. The stochastic dynamics have the given infinitesimal standard deviation
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the additive white noise (standard deviation)
#'
#' @return The time series
#' @export
realistic <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {

  # Avoid wrong inputs
  # The current implementation doesn't use any tpert different than 0
  # TODO: implement this functionality
  if(tpert != 0) { stop("Currently tpert different than 0 is not supported by this method") }

  # Translate parameters to the language of differential equations
  y0 <- pert # The perturbation represents the initial condition
  r <- log(2)/thalf # Translate the half-life to a multiplicative constant

  # Scale the noise
  # We want to interpret the parameters in piecewise, exponential and realistic uniformly.
  # Parameter noise is of course no exception. In order to make sure that the standard deviation
  # of the resulting time series is (roughly) equal to the parameter noise, we need to rescale
  # the infinitesimal standard deviation term used in the stochastic differential equation
  tStep <- max(diff(t)) # Numerical time step
  isd <- noise/sqrt(tStep) # Infinitesimal standard deviation

  # Pose the differential equation dy = f(y,t) dt + g(y,t) dW
  #
  # With
  # f(y, t) = -r * y (exponential decay)
  # and
  # g(y, t) = s (white noise)
  #
  # Unfortunately, the package sde uses a very obscure syntax. Instead of functions it expects
  # expressions depending on x and t as an input. When those object contain, additionally,
  # parameters, we need to pass them via the substitute command.
  #
  # e.g: substitute(a + x, list(a = 2)) returns 2 + x
  # 2 + x is an object of class call, that must be converted to expression
  f <- as.expression(
                      substitute(-r * x,
                                 list(r = r))
                    )
  g <- as.expression(
                     substitute(s,
                                list(s = isd))
                     )

  # Solve
  sol <- sde::sde.sim(X0 = y0,
                      T = max(t),
                      N = length(t)-1,
                      drift = f,
                      sigma = g,
                      sigma.x = 0.0,
                      method = 'euler')

  # Extract the state only (so the output has the same structure as in piecewise and exponential)
  sol <- as.data.frame(sol)
  y <- as.numeric(sol$x)

  # Don't forget to add the offset
  y <- y + offset

  return(y)
}

#' Simulate one time series with disturbance
#'
#' @param nyr number of years that need to be simulated
#' @param nobsyr number of observations per year that will be simulated
#' @param tMiss timing of missing values [observation number]. If tMiss equals NA, no missing values are introduced.
#' @param nDr number of drought years that are introduced [i.e. setting seasonality of a year equal to its minimum value]. These drought years are randomly chosen.
#' @param seasAv average seasonality profile
#' @param seasAmp seasonality amplitude
#' @param trAv offset value of time series
#' @param remSd standard deviation of the remainder
#' @param remMod ARMA model of remainder
#' @param distMag magnitude of the disturbance
#' @param distT timing of the disturbance [observation number]
#' @param distRec duration of the recovery [number of observations]
#' @param distType type of disturbance-recovery process: 'piecewise' represents a step function with linear recovery, 'exponential' represents an exponential decay
#'
#' @return  a list containign the years for which a drought was introduced and a time series object, containing the simulated seasonality, trend, remainder, disturbance, and the sum of these components.
#' @export
simulTS <- function(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remMod, distMag, distT, distRec, distType){
  #-------------------------------------------------
  # simulate seasonality
  simSeas <- rep(as.numeric(seasAv[1:(nobsyr*2)]), times=ceiling(nyr/2))
  if (ceiling(nyr/2) != (nyr/2)){
    simSeas <- simSeas[1:(length(simSeas)-nobsyr)]
  }
  simSeas <- simSeas/max(seasAv)*seasAmp

  # introduce drought years (seasonality equal to minimum value)
  if (nDr==0){ydr <- 0}
  if (nDr > 0){
    ydr <- sample(1:(nyr-1), nDr)
    #print(ydr)
    for (i in 1:length(ydr)){
      offs <- which(simSeas[1:nobsyr]==min(simSeas[1:nobsyr], na.rm=T))
      strt <- offs[1] + ((ydr[i]-1) * nobsyr)
      endt <- offs[1] + ydr[i]*nobsyr
      simSeas[strt : endt] <- -seasAmp
    }
  }

  #-------------------------------------------------
  # simulate offset
  simTr <- rep(trAv, times=(nyr*nobsyr))

  #-------------------------------------------------
  # simulate remainder
  simRem <- arima.sim(model = remMod, n = nobsyr*nyr, sd = remSd)
  simRem <- (simRem - mean(simRem))# zero mean
  simRem <- simRem/sd(simRem)*remSd # set standard deviation

  #-------------------------------------------------
  # simulate disturbance

  if(distType == 'piecewise'){
    simDist <- piecewise(1:(nobsyr*nyr), pert=distMag, tpert = distT, thalf = distRec)
  }

  if(distType == 'exponential'){
    simDist <- exponential(1:(nobsyr*nyr), pert=distMag, tpert = distT, thalf = distRec)
  }

  #-------------------------------------------------
  # Sum components
  simTS <-simSeas+simTr+simRem+simDist
  # set missing values
  if (is.na(tMiss[1])==F){
    simTS[tMiss] <- NA
  }
  list(ydr, ts(t(rbind(simSeas, simTr, simRem, simDist, simTS)), frequency = nobsyr))
}


#' Simulation of nrep disturbance time series.
#'
#' @param nrep number of time series to simulate
#' @param nyr number of years that need to be simulated
#' @param nobsYr number of observations per year that will be simulated
#' @param nDr number of drought years that are introduced [i.e. setting seasonality of a year equal to its minimum value]. These drought years are randomly chosen for each of the simulated time series.
#' @param seasAv average seasonality profile
#' @param seasAmp seasonality amplitude
#' @param trAv offset value of time series
#' @param remSd standard deviation of the remainder
#' @param distMaglim limits of the disturbance magnitude, should be a vector with the minimum and maximum value. If the minimum equals the maximum value, the disturbance magnitude is fixed for each simulated time series (and equal to the minimum value). When the minimum value does not equal the maximum value, a disturbance magnitude is randomly chosen in the given interval for each simulated time series.
#' @param distTy year of the disturbance. If distTy equals one, the disturbance will take place in the first year. The exact disturbance date (day or year) is randomly chosen per time series.
#' @param distReclim limits of the halftime period of the recovery [number of observations], should be a vector with the minimum and maximum value. If the minimum equals the maximum value, the recovery period is fixed for each simulated time series (and equal to the minimum value). When the minimum value does not equal the maximum value, a recovery period is randomly chosen in the given interval for each simulated time series.
#' @param remcoef list of ARMA models of remainder component. For each simulated time series, a remainder model is randomly chosen from the list.
#' @param mval number of missing values to be introduced. If mval equals NA, no missing values are introduced. For missing values with a random interval (see mvaldist), this should equal the fraction of missing values (mval equal to 0.1 will result in an NA value for 10 percent of the time series). For missing values having a regular interval, every mval observations one value is kept (eg for a daily time series, a mval equal to 5 will result in one observation every 5 days).
#' @param mvaldist the distribution of the missing values. Should equal 'random'or 'interval'.
#' @param distType the type of disturbance. piecewise refers to a linear decay function, while exponential refers to an exponential decay
#'
#' @return a list with the simulated time series, offest, seasonality, remainder, disturbance component and parameters used for the simulation. The time series (components) are stored as matrix where each row is a time series and the columns are associated with the observation numbers.
#' @export
simulCase <- function(nrep, nyr, nobsYr, nDr, seasAv, seasAmp,
                      trAv, remSd, distMaglim, distTy, distReclim, remcoef, mval, mvaldist, distType){
  # matrices to store the time series
  TSsim <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # time series
  TSsimTr <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # trend
  TSsimSeas <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # seasonality
  TSsimRem <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # remainder
  TSsimDist <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # disturbance
  # variables to store settings
  m_remcoef <- list() # ARMA model of remainders
  m_year_dr <- list() # drought year(s)
  m_nyr <- matrix(NA, nrow = nrep, ncol = 1)
  m_nobsYr <- matrix(NA, nrow = nrep, ncol = 1)
  m_nDr <- matrix(NA, nrow = nrep, ncol = 1)
  m_seasAmp <- matrix(NA, nrow = nrep, ncol = 1)
  m_trAv <- matrix(NA, nrow = nrep, ncol = 1)
  m_remSd <- matrix(NA, nrow = nrep, ncol = 1)
  m_distMag <- matrix(NA, nrow = nrep, ncol = 1)
  m_distT <- matrix(NA, nrow = nrep, ncol = 1)
  m_distRec <- matrix(NA, nrow = nrep, ncol = 1)
  m_distType<- matrix(NA, nrow = nrep, ncol = 1)
  m_mVal<- matrix(NA, nrow = nrep, ncol = 1)


  for (ii in 1:nrep){
    # randomly select a remainder model
    modi <- sample((1:length(remcoef)), 1)
    # randomly select a disturbance magnitude within the given limits
    if(distMaglim[1]==distMaglim[2]){
      distMag <- distMaglim[1]
    }else{distMag <- sample((10*distMaglim[1]):(10*distMaglim[2]), 1)/10}
    # ranomly select a recovery period within the given limits
    if(distReclim[1]==distReclim[2]){
      distRec <- distReclim[1]
    }else{distRec <- sample(distReclim[1]:distReclim[2],1)}

    # randomly select a disturbance day within the defined disturbance year
    distT <- ((distTy-1)*nobsYr) + sample((1:nobsYr), 1)

    # missing values to be introduced
    if (is.na(mval)){tMiss = NA}else{
      if(mvaldist == 'interval'){
        tMiss <- 1:(nyr*nobsYr) # missing values
        tMiss <- tMiss[-seq(1,(nyr*nobsYr),mval)]
      }
      if(mvaldist == 'random'){
        tMiss <- sample(1:(nyr*nobsYr), round(nyr*nobsYr*mval))
      }
    }
    # simulate time series
    sts <- simulTS(nyr, nobsYr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remcoef[[modi]], distMag, distT, distRec, distType)
    # store simulated time series
    TSsim[ii,] <- sts[[2]][,5]
    TSsimTr[ii,] <- sts[[2]][,2]
    TSsimSeas[ii,] <- sts[[2]][,1]
    TSsimRem[ii,] <- sts[[2]][,3]
    TSsimDist[ii,] <- sts[[2]][,4]
    # store selected parameters
    m_remcoef[[ii]] <-  remcoef[[modi]]
    m_year_dr[[ii]] <- sts[[1]]
    m_nyr[ii] <- nyr
    m_nobsYr[ii] <- nobsYr
    m_nDr[ii] <- nDr
    m_seasAmp[ii] <- seasAmp
    m_trAv[ii] <-trAv
    m_remSd[ii] <- remSd
    m_distMag[ii] <- distMag
    m_distT[ii] <-distT
    m_distRec[ii] <- distRec
    m_distType[ii] <- distType
    m_mVal[ii] <- mval
    rm(modi, distT, sts)
  }

  # Dataframe selected parameters
  TSsimParam <- data.frame(m_nyr,
                           m_nobsYr,
                           m_nDr,
                           m_seasAmp,
                           m_trAv,
                           m_remSd,
                           m_distMag,
                           m_distT,
                           m_distRec,
                           m_distType,
                           m_mVal)
  names(TSsimParam) <- c("number_yrs", "obs_per_year", "number_droughts",
                         "seas_amp", "trend_av", 'rem_sd',
                         'dist_magn', 'dist_time', 'dist_rec', 'dist_type', 'miss_val')
  TSsimParam$year_drought <- m_year_dr
  TSsimParam$rem_coef <- m_remcoef

  out <- list( TSsim, TSsimTr, TSsimSeas, TSsimRem, TSsimDist, TSsimParam)
  names(out) <- c('timeSeries', 'Trend', 'Seasonality', 'Remainder', 'Disturbance', 'Parameters')
  out
}

