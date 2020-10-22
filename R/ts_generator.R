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
      (offset + pert + m *(t - tpert))    * (t >= tpert) * (t <= tpert + ttrans) + # Transitory regime
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
#' 3. Gaussian noise is added to the final time series to simulate measurement errors
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the stochastic term in the differential equation (standard deviation of the integrated time series)
#'
#' @return The time series
#' @export
realistic <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {

  # Translate parameters to the language of differential equations
  y0 <- pert # The perturbation represents the initial condition
  r <- log(2)/thalf # Translate the half-life to a multiplicative constant
  sigma <- noise * sqrt(2 * log(2) / thalf) # Infinitesimal standard deviation
  # The infinitesimal standard deviation `sigma` yields a standard deviation of magnitude `noise` after integration
  # Reference: https://math.stackexchange.com/questions/2558659/expectation-and-variance-of-stochastic-differential-equations

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
                                list(s = sigma))
                     )


  # Solve
  sol <- sde::sde.sim(X0 = y0,
                      T = max(t),
                      N = length(t)-1,
                      drift = f,
                      sigma = g,
                      sigma.x = 0.0,
                      method = 'euler')


  # Shift the time series (only if tpert is not zero)
  if(tpert != 0) {
    # Shift the original time series
    ts <- time(sol) # Store the original times
    i <- min(which(ts >= tpert)) # Find the index corresponding to tpert
    sol <- lag(sol, -i + 1) # Displace the time series, so it begins at tpert

    # Create the time series before tpert (only dynamic noise around equilibrium)
    tfill <- seq(ts[1], ts[i-1], by = 1 / frequency(sol))
    fill <- ts(realistic(tfill, noise = noise), start = ts[1], end = ts[i-1], frequency = frequency(sol))

    # Paste both time series together
    sol <- ts(c(fill, sol), start = start(fill), frequency = frequency(fill))

    # Trim the tail, so the displaced time series is equal in size to the original
    sol <- head(sol, length(sol) - length(fill))
  }

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
simulTS <- function(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, distMag, distT, distRec, distType){
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
    for (i in 1:length(ydr)){
      offs <- which(simSeas[1:nobsyr]==min(simSeas[1:nobsyr], na.rm=T))
      strt <- offs[1] + ((ydr[i]-1) * nobsyr)
      endt <- offs[1] + ydr[i]*nobsyr
      simSeas[strt : endt] <- -seasAmp
    }
  }

    #-------------------------------------------------
  # simulate disturbance-recovery

  if(distType == 'piecewise'){
    simDist <- piecewise(1:(nobsyr*nyr), offset = trAv, pert=distMag, tpert = distT, thalf = distRec, noise = remSd)#disturbance with noise component
    simTruth <- piecewise(1:(nobsyr*nyr), offset = trAv, pert=distMag, tpert = distT, thalf = distRec, noise = 0)# disturbance without noise compontent = truth
  }

  if(distType == 'exponential'){
    simDist <- exponential(1:(nobsyr*nyr), offset = trAv, pert=distMag, tpert = distT, thalf = distRec, noise = remSd)
    simTruth <- exponential(1:(nobsyr*nyr), offset = trAv, pert=distMag, tpert = distT, thalf = distRec, noise = 0)
  }

  if(distType == 'diffEq'){
    simDist <- realistic(1:(nobsyr*nyr), offset = trAv, pert=distMag, tpert = distT, thalf = distRec, noise = remSd)
    simTruth <- realistic(1:(nobsyr*nyr), offset = trAv, pert=distMag, tpert = distT, thalf = distRec, noise = 0)
  }

  #-------------------------------------------------
  # Sum components
  simTS <-simSeas+simDist
  # set missing values
  if (is.na(tMiss[1])==F){
    simTS[tMiss] <- NA
  }
  list(ydr, ts(t(rbind(simSeas, simDist, simTruth, simTS)), frequency = nobsyr))
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
#' @param mval number of missing values to be introduced. If mval equals NA, no missing values are introduced. For missing values with a random interval (see mvaldist), this should equal the fraction of missing values (mval equal to 0.1 will result in an NA value for 10 percent of the time series). For missing values having a regular interval, every mval observations one value is kept (eg for a daily time series, a mval equal to 5 will result in one observation every 5 days).
#' @param mvaldist the distribution of the missing values. Should equal 'random'or 'interval'.
#' @param distType the type of disturbance. piecewise refers to a linear decay function, while exponential refers to an exponential decay
#'
#' @return a list with the simulated time series, offest, seasonality, remainder, disturbance component and parameters used for the simulation. The time series (components) are stored as matrix where each row is a time series and the columns are associated with the observation numbers.
#' @export
simulCase <- function(nrep, nyr, nobsYr, nDr, seasAv, seasAmp,
                      trAv, remSd, distMaglim, distTy, distReclim, mval, mvaldist, distType){
  # matrices to store the time series
  TSsim <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # time series
  TSsimSeas <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # seasonality
  TSsimTruth <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # disturbance without noise
  TSsimDist <- matrix(NA, nrow=nrep, ncol=(nyr*nobsYr)) # disturbance with noise
  # variables to store settings
  # m_remcoef <- list() # ARMA model of remainders
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
    # modi <- sample((1:length(remcoef)), 1)
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
    sts <- simulTS(nyr, nobsYr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, distMag, distT, distRec, distType)

    # store simulated time series
    TSsim[ii,] <- sts[[2]][,4]
    TSsimSeas[ii,] <- sts[[2]][,1]
    TSsimTruth[ii,] <- sts[[2]][,3]
    TSsimDist[ii,] <- sts[[2]][,2]
    # store selected parameters
    # m_remcoef[[ii]] <-  remcoef[[modi]]
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
    rm(distT, sts)
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
  # TSsimParam$rem_coef <- m_remcoef

  out <- list( TSsim, TSsimSeas, TSsimTruth, TSsimDist, TSsimParam)
  names(out) <- c('timeSeries', 'Seasonality', 'Truth', 'Disturbance', 'Parameters')
  out
}


#' Create list of parameter values.
#'
#' @param sttngs: settings list defining the parameter set-up
#'
#' @return list of parameter settings
#' @export
#'
setParamValues <- function(sttngs){
  if(sttngs$general$parSetUp == 'avg'){
    out <- set_avg_par(sttngs)
  } else if(sttngs$general$parSetUp == 'dist'){
    out <- set_dist_par(sttngs)
  } else if(sttngs$general$parSetUp == 'comb'){
    out <- set_combine_par(sttngs)
  }else if(sttngs$general$parSetUp == 'int'){
    out <- set_interval_par(sttngs)
  }
  out
}

#' Get average parameter values. Helper function for the set_avg_par() function
#'
#' @param param parameter list containing the (i) 'type' = type of parameter, (ii) 'vals' = evaluated parameter values and (iii) 'obs' = the observed parameter values.
#' @param nval number of values to be provided
#'
#' @return a vector containing nval parameter values
#' @export
#'
get_avg_val <- function(param, nval){
  # if observed values of the parameter are available (type = dist), return the mean of these observed values
  if(param$type == 'dist'){
    out <- rep(mean(param$obs, na.rm = T), nval)
    # if only a range of reasonable values for the parameter are known (type = range), get the mean value of the range
  }else if(param$type == 'range'){
    out <- rep(mean(c(min(param$vals, na.rm = T), max(param$vals, na.rm = T))), nval)
    # if the parameter only contains categorical values, randomly sample the categorical values
  }else if(param$type == 'cat'){
    out <- sample(param$vals,nval, replace = T)
  }
  out
}

#' Sample parameter values based on their distribution
#'
#' @param param parameter list containing the (i) 'type' = type of parameter, (ii) 'vals' = evaluated parameter values and (iii) 'obs' = the observed parameter values.
#' @param nval number of values to be provided
#'
#' @return a vector containing nval parameter values
#' @export
#'
get_dist_val <- function(param, nval){
  # if observed values of the parameter are available (type = dist), sample values using the frequencies of the observed values
  if(param$type == 'dist'){
    h <- hist(param$obs,breaks = (seq(min(param$obs),max(param$obs), l = 100)))
    out <- sample(h$mids, nval, replace = T, prob = h$counts)
    # if only a range of reasonable values for the parameter are known (type = range), randomly sample observations within the range
  }else if(param$type == 'range'){
    if(min(param$vals, na.rm = T) == max(param$vals, na.rm = T)){
      out <- rep(min(param$vals, na.rm = T), nval)
    }else{
      out <- sample(seq(min(param$vals, na.rm = T), max(param$vals, na.rm = T), l = 100000),nval, replace = T)
    }
    # if the parameter only contains categorical values, randomly sample these categorical values
  }else if(param$type == 'cat'){
    out <- sample(param$vals,nval, replace = T)
  }
  out
}

#' Sample parameter values in a particular interval
#'
#' @param param parameter list containing the (i) 'type' = type of parameter, (ii) 'vals' = evaluated parameter values and (iii) 'obs' = the observed parameter values and (iv) 'fix' = interval of the fixed observations.
#' @param nval number of values to be provided
#'
#' @return a vector containing nval parameter values
#' @export
#'
get_int_val <- function(param, nval){
  # if observed values of the parameter are available (type = dist), sample values using the frequencies of the observed values
  if((param$type == 'range') | (param$type == 'dist')){
    if(min(param$fix, na.rm = T) == max(param$fix, na.rm = T)){
      out <- rep(min(param$fix, na.rm = T), nval)
    }else{
      out <- sample(seq(min(param$fix, na.rm = T), max(param$fix, na.rm = T), l = 100000),nval, replace = T)
    }
    # if the parameter only contains categorical values, randomly sample these categorical values
  }else if(param$type == 'cat'){
    out <- sample(param$fix,nval, replace = T)
  }
  out
}

#' Set list of parameter values. Each evaluated parameter is stepwise altered over a predefined set of values. While evaluating parameter x, the other parameters are kept constant and equal to their mean value.
#'
#' @param sttngs settings file for the time series simulation set-up
#'
#' @return a list containing the parameter values for each evaluated parameter
#' @export
#'
set_avg_par <- function(sttngs){
  pars <- list()

  for (i in 1:length(sttngs$general$eval)){#iterate over the parameters that need to be evaluated
    vari <- sttngs$general$eval[i]# parameter that needs to be evaluated
    varvali <- sttngs[[vari]]$vals# values of the parameter that needs to be evaluated


    for(ii in 1:length(varvali)){# iterate over the values of the evaluated parameter
      nTS <- sttngs$general$nTS# number of time series to be simulated per evaluated parameter value
      pari <- list(nrep = rep(1,nTS),#number of time series to be simulated per parameter combination
                   nyr = round(get_avg_val(sttngs$nyr,nTS)),# number of years to be simulated
                   nobsYr = rep(sttngs$general$nobsYr,nTS),# number of observations per year to be simulated
                   nDr = floor(get_avg_val(sttngs$nDr,nTS)),# number of drought years
                   seasAv = list(sttngs$general$seasAv),# seasonal average values
                   seasAmp = get_avg_val(sttngs$seasAmp,nTS),# seasonal amplitude
                   trAv = get_avg_val(sttngs$trAv,nTS),# offset
                   remSd = get_avg_val(sttngs$remSd,nTS),# standard deviation of the remainder
                   distMag = get_dist_val(sttngs$distMag,nTS),# magnitude of the disturbance
                   distT = get_avg_val(sttngs$distT,nTS),# timing disturbance
                   distRec = get_dist_val(sttngs$distRec,nTS),# recovery period after disturbance
                   missVal = get_avg_val(sttngs$missVal,nTS),# fraction of missing values
                   DistMissVal = get_avg_val(sttngs$DistMissVal,nTS),# distribution of missing values
                   distType = get_avg_val(sttngs$distType,nTS)# type of recovery
      )
      pari[[vari]] <- rep(varvali[ii],nTS)
      pars[[vari]][[as.character(varvali[ii])]] <- pari
    }
  }
  pars
}


#' Set list of parameter values. Each evaluated parameter is stepwise altered over a predefined set of values. While evaluating parameter x, all combinations of predifined parameter values for all other parameters are made.
#'
#' @param sttngs settings file for the time series simulation set-up
#'
#' @return a list containing the parameter values for each evaluated parameter
#' @export
#'
set_combine_par <- function(sttngs){
  pars <- list()

  for (i in 1:length(sttngs$general$eval)){#iterate over the parameters that need to be evaluated
    vari <- sttngs$general$eval[i] # parameter to be evaluated
    varvali <- sttngs[[vari]]$vals  # values of the parameter to be evaluated
    vars <- names(sttngs)[names(sttngs) != 'general' & names(sttngs) != vari]# other parameters that should be combined

    for(ii in 1:length(varvali)){# iterate over the values of the evaluated parameter
      comb <- expand.grid(lapply(sttngs[vars], function(x){x$vals}))# combine all parameter values, except the evaluated one
      comb[[vari]] <- rep(varvali[ii],length(comb[[1]]))# add value of parameter to be evaluated to the list
      comb$nrep <- rep(1,length(comb[[1]])) # set the number of repetitions
      comb$nobsYr <- rep(sttngs$general$nobsYr,length(comb[[1]])) # set the number of observations per year
      comb$seasAv <- list(sttngs$general$seasAv)# set the seasonality

      pars[[vari]][[as.character(varvali[ii])]] <- comb # add the parameter values to the parameter list
    }
  }
  pars
}

#' Set list of parameter values. Each evaluated parameter is stepwise altered over a predefined set of values. While evaluating parameter x, the other parameters are sampled using their observed distribution.
#'
#' @param sttngs settings file for the time series simulation set-up
#'
#' @return a list containing the parameter values for each evaluated parameter
#' @export
#'
set_dist_par <- function(sttngs){
  pars <- list()

  for (i in 1:length(sttngs$general$eval)){#iterate over the parameters that need to be evaluated
    vari <- sttngs$general$eval[i]# parameter to be evaluated
    varvali <- sttngs[[vari]]$vals# values of the parameter to be evaluated


    for(ii in 1:length(varvali)){# iterate over the values of the evaluated parameter
      nTS <- sttngs$general$nTS
      pari <- list(nrep = rep(1,nTS),
                   nyr = round(get_dist_val(sttngs$nyr,nTS)),
                   nobsYr = rep(sttngs$general$nobsYr,nTS),
                   nDr = floor(get_dist_val(sttngs$nDr,nTS)),
                   seasAv = list(sttngs$general$seasAv),
                   seasAmp = get_dist_val(sttngs$seasAmp,nTS),
                   trAv = get_dist_val(sttngs$trAv,nTS),
                   remSd = get_dist_val(sttngs$remSd,nTS),
                   distMag = get_dist_val(sttngs$distMag,nTS),
                   distT = get_dist_val(sttngs$distT,nTS),
                   distRec = get_dist_val(sttngs$distRec,nTS),
                   missVal = get_dist_val(sttngs$missVal,nTS),
                   DistMissVal = get_dist_val(sttngs$DistMissVal,nTS),
                   distType = get_dist_val(sttngs$distType,nTS)
      )
      pari[[vari]] <- rep(varvali[ii],nTS)
      pars[[vari]][[as.character(varvali[ii])]] <- pari
    }
  }
  pars
}

#' Set list of parameter values. Each evaluated parameter is stepwise altered over multiple intervals. While evaluating parameter x, the other parameters are sampled.
#'
#' @param sttngs settings file for the time series simulation set-up
#'
#' @return a list containing the parameter values for each evaluated parameter
#' @export
#'
set_interval_par <- function(sttngs){
  pars <- list()

  for (i in 1:length(sttngs$general$eval)){#iterate over the parameters that need to be evaluated
    vari <- sttngs$general$eval[i]# parameter that needs to be evaluated
    varvali <- sttngs[[vari]]$vals# values of the parameter that needs to be evaluated
    # if(sttngs[[vari]]$type == 'dist'){
      varvali <- matrix(varvali,nrow=2,ncol=length(varvali)/2)
      len <- length(varvali)/2
    # } else{len <- length(varvali)}

    for(ii in 1:len){# iterate over the values of the evaluated parameter
      nTS <- sttngs$general$nTS# number of time series to be simulated per evaluated parameter value
      pari <- list(nrep = rep(1,nTS),#number of time series to be simulated per parameter combination
                   nyr = round(get_int_val(sttngs$nyr,nTS)),# number of years to be simulated
                   nobsYr = rep(sttngs$general$nobsYr,nTS),# number of observations per year to be simulated
                   nDr = floor(get_int_val(sttngs$nDr,nTS)),# number of drought years
                   seasAv = list(sttngs$general$seasAv),# seasonal average values
                   seasAmp = get_int_val(sttngs$seasAmp,nTS),# seasonal amplitude
                   trAv = get_int_val(sttngs$trAv,nTS),# offset
                   remSd = get_int_val(sttngs$remSd,nTS),# standard deviation of the remainder
                   distMag = get_int_val(sttngs$distMag,nTS),# magnitude of the disturbance
                   distT = get_int_val(sttngs$distT,nTS),# timing disturbance
                   distRec = get_int_val(sttngs$distRec,nTS),# recovery period after disturbance
                   missVal = get_int_val(sttngs$missVal,nTS),# fraction of missing values
                   DistMissVal = get_int_val(sttngs$DistMissVal,nTS),# distribution of missing values
                   distType = get_int_val(sttngs$distType,nTS)# type of recovery
                   )
      # if (sttngs[[vari]]$type == 'dist'){
        pari[[vari]] <- sample(seq(min(varvali[,ii], na.rm = T), max(varvali[,ii], na.rm = T), l = 100000),nTS, replace = T)
        pars[[vari]][[as.character(mean(varvali[,ii]))]] <- pari
      # }
      # else{
      #   pari[[vari]] <- rep(varvali[ii],nTS)
      #   pars[[vari]][[as.character(varvali[ii])]] <- pari
      # }
    }
  }
  pars
}

#' #' #' Sample parameter values over an interval
#' #'
#' #' @param param parameter list containing the (i) 'type' = type of parameter, (ii) 'vals' = evaluated parameter values and (iii) 'obs' = the observed parameter values.
#' #' @param nval number of values to be provided
#' #'
#' #' @return a vector containing nval parameter values
#' #' @export
#' #'
#' get_int_val <- function(param, nval){
#'   # if observed values of the parameter are available (type = dist), sample values using the frequencies of the observed values
#'   if(param$type == 'dist'){
#'     vals <- matrix(param$vals,nrow=2,ncol=length(param$vals)/2)
#'     vals <- vals[,ceiling((1+dim(vals)[2])/2)]
#'     out <- sample(seq(min(vals, na.rm = T), max(vals, na.rm = T), l = 100000),nval, replace = T)
#'     # if only a range of reasonable values for the parameter are known (type = range), randomly sample observations within the range
#'   }else if(param$type == 'range'){
#'     if(min(param$vals, na.rm = T) == max(param$vals, na.rm = T)){
#'       out <- rep(min(param$vals, na.rm = T), nval)
#'     }else{
#'       out <- sample(seq(min(param$vals, na.rm = T), max(param$vals, na.rm = T), l = 100000),nval, replace = T)
#'     }
#'     # if the parameter only contains categorical values, randomly sample these categorical values
#'   }else if(param$type == 'cat'){
#'     out <- sample(param$vals,nval, replace = T)
#'   }
#'   out
#' }
