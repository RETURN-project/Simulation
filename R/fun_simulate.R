#' ---
#' title: "Simulation functions"
#' author: "Wanda De Keersmaecker"
#' date: "2/5/2020"
#' output: html_document
#' ---
#'
## ----setup, include=FALSE------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

#'
#' # 1. Time series characterisation
#'
#' ## 1.1 Decompose time series into trend, seasonality and remainder
#' This function decomposes time series into three components using BFAST01 functionality: trend, seasonality and remainder. Trends are fitted using linear regression without breaks, seasonality is fitted using a first order harmonic function and the remainder equals the anomalies (i.e. time series - trend - seasonality).
#'
#' Function inputs:
#'
#' * df: a dataframe with time series that need to be decomposed. The dataframe needs to be structured as follows: each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the time series values for each observation date.
#' * nyr: number of years of the input time series
#' * nobsYr: number of observations per year of the input time series</font>
#'
#' Function outputs:
#'
#' * Seasonality: a dataframe with the seasonality of each pixel. Each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the seasonality values for each observation date.
#' * Remainder: a dataframe with the remainder of each pixel. Dataframe is structured in the same way as the seasonality.
#' * Trend: a dataframe with the trend of each pixel. Dataframe is structured in the same way as the seasonality.
#' * Seasonality_coefficients: a dataframe with the coeficients of the fitted harmonic function. Each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the coefficients of the fitted harmonic function.
#'
#'
#'
## ------------------------------------------------------------------------
decompTSbfast <- function(df, nyr, nobsYr){
  # Initialize each output data frame with NA values
  # seasonality
  dfSeas <- data.frame(data= matrix(NA, nrow = dim(df)[1],ncol= dim(df)[2]))
  colnames(dfSeas) <-  c('lat', 'lon', colnames(df)[-c(1,2)])
  # remainder
  dfRem <- data.frame(data= matrix(NA, nrow = dim(df)[1],ncol= dim(df)[2]))
  colnames(dfRem) <- c('lat', 'lon', colnames(df)[-c(1,2)])
  # trend
  dfTr <- data.frame(data= matrix(NA, nrow = dim(df)[1],ncol= dim(df)[2]))
  colnames(dfTr) <- c('lat', 'lon', colnames(df)[-c(1,2)])
  # coef harmonic
  dfSeasCoef <- data.frame(data= matrix(NA, nrow = dim(df)[1],ncol= (4)))
  colnames(dfSeasCoef) <-  c('lat', 'lon', 'cos', 'sin')

  library(bfast)
  # Iterate over each pixel and decompose time series
  for (i in 1:dim(df)[1]){
    err <- 0
    # Create time series object
    tss <- ts(as.numeric(df[i,-c(1,2)]), frequency =nobsYr )
    # Use BFAST01 to decompose time series
    tryCatch({
      tmp <- bfast01(tss,
                     formula = response ~ trend + harmon, order=1, level = 0.0001)#'periodic'
    }, error = function(e) {
      if(e$message == "inadmissable change points: 'from' is larger than 'to'"){
        #print(i)
        err <- 1
      }
    })
    if(err == 1){next}
    # coeficients of fitted trend and seasonality
    coef<-tmp$model[[1]]$coefficients
    # fitted seasonality and trend
    fit <- coef[1]+coef[2]* 1:length(tss) + coef[3]*cos(2*pi*time(tss))+ coef[4]*sin(2*pi*time(tss))
    # remainder
    rem <- as.numeric(df[i,-c(1,2)]) - as.numeric(fit)

    dfSeasCoef[i,] <- c(df$lat[i], df$lon[i],coef[3],coef[4])
    dfSeas[i,] <- c(df$lat[i], df$lon[i], coef[3]*cos(2*pi*time(tss))+ coef[4]*sin(2*pi*time(tss)))
    dfRem[i,] <- c(df$lat[i], df$lon[i], rem)
    dfTr[i,] <- c(df$lat[i], df$lon[i],coef[1]+coef[2]* 1:length(tss))

  }
  lst <- list(dfSeas,dfRem,dfTr, dfSeasCoef)
  names(lst) <- c('Seasonality', 'Remainder', 'Trend', 'Seasonality_coefficients')
  lst
}

#'
#'
#' ## 1.2 Fit ARMA model
#' This function automatically fits an ARMA model without seasonal component.
#'
#' Function inputs:
#'
#' * tsx: a time series object for which an ARMA model needs to be fitted.
#'
#' Function outputs:
#'
#' * ARMA model coefficients. A list containing the ARMA coefficients and their order (number of AR coefficients, non-seasonal differences and MA coefficients
#'
## ------------------------------------------------------------------------
# ARMA coefficients
getARMAcoef <- function(tsx){
  library(forecast)
  arm <- auto.arima(tsx, seasonal=F) #arm$arma #A compact form of the specification, as a vector giving
  #the number of AR, MA, seasonal AR and seasonal MA coefficients,
  #plus the period, and the number of non-seasonal, and seasonal differences.
  teller <- 0
  coefmod <- list()
  if(arm$arma[1]>0){
    coefmod <- c(coefmod, ar = as.numeric(arm$coef[1:arm$arma[1]]))
    teller <- teller + arm$arma[1]
  }
  if(arm$arma[2]>0){
    coefmod <- c(coefmod, ma = as.numeric(arm$coef[(1+teller):(teller+arm$arma[2])]))
    #teller <- teller + arm$arma[2]
  }
  coefmod <- c(coefmod, order = as.numeric(c(arm$arma[1], arm$arma[6], arm$arma[2])))
  coefmod
}

#'
#' # 2. Time series simulation
#'
#' ## 2.1 Simulate one time series with disturbance
#' Function to simulate a time series with disturbance.
#'
#' Function inputs:
#'
#' * nyr = number of years that need to be simulated
#' * nobsyr = number of observations per year that will be simulated
#' * tMiss = timing of missing values [observation number]. If tMiss equals NA, no missing values are introduced.
#' * nDr = number of drought years that are introduced [i.e. setting seasonality of a year equal to its minimum value]. These drought years are randomly chosen.
#' * seasAv= average seasonality profile
#' * seasAmp = seasonality amplitude
#' * trAv = offset value of time series
#' * remSd = standard deviation of the remainder
#' * remMod = ARMA model of remainder
#' * distMag = magnitude of the disturbance
#' * distT = timing of the disturbance [observation number]
#' * distRec = duration of the recovery [number of observations]</font>
#'
#' Function outputs:
#'
#' * years for which a drought was introduced
#' * a time series object, containing the simulated seasonality, trend, remainder, disturbance, and the sum of these components.
#'
## ----simults-------------------------------------------------------------
simulTS <- function(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remMod, distMag, distT, distRec){
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
  simDist <- simulDist(distT, distRec, distMag, nobsyr*nyr)
  #if ((distT+distRec) <= (nobsyr*nyr)){
  # simDist <- c(rep(0,(distT-1)), seq(distMag, 0, length.out = distRec), rep(0, ((nobsyr*nyr)-distRec -distT +1)))
  #}
  #if ((distT+distRec) > (nobsyr*nyr)){
  # simDist <- c(rep(0,(distT-1)), seq(distMag, 0, length.out = distRec))
  #simDist <- simDist[1:(nobsyr*nyr)]
  #}
  #-------------------------------------------------
  # Sum components
  simTS <-simSeas+simTr+simRem+simDist
  # set missing values
  if (is.na(tMiss[1])==F){
    simTS[tMiss] <- NA
  }
  list(ydr, ts(t(rbind(simSeas, simTr, simRem, simDist, simTS)), frequency = nobsyr))
}

#' ## 2.2 Simulate disturbance component
#'
#' Simulation of a disturbance as a step function with linear recovery trend. The time series equals zero when it is not being disturbed. The disturbance introduces at observation number *distT* a negative value with magnitude *distMag*. The linear recovery takes *distRec* observations.
#'
#' Function inputs:
#'
#' * distT = timing of the disturbance [observation number]
#' * distRec = duration of the recovery period [number of observations]
#' * distMag = magnitude of the disturbance
#' * nobs = number of observations per year that will be simulated
#'
#' Function outputs:
#'
#' * a vector with the simulated disturbance component
#'
#'
## ----simulDist-----------------------------------------------------------
simulDist <- function(distT, distRec, distMag, nobs){
  if ((distT+distRec) <= (nobs)){
    simDist <- c(rep(0,(distT-1)), seq(distMag, 0, length.out = distRec), rep(0, ((nobs)-distRec -distT +1)))
  }
  if ((distT+distRec) > (nobs)){
    simDist <- c(rep(0,(distT-1)), seq(distMag, 0, length.out = distRec))
    simDist <- simDist[1:(nobs)]
  }
  simDist
}


#' ## 2.3 Simulate a set of time series
#'
#' Simulation of *nrep* disturbance time series.
#'
#' Function inputs:
#'
#' * nrep = number of time series to simulate
#' * nyr = number of years that need to be simulated
#' * nobsyr = number of observations per year that will be simulated
#' * nDr = number of drought years that are introduced [i.e. setting seasonality of a year equal to its minimum value]. These drought years are randomly chosen for each of the simulated time series.
#' * seasAv= average seasonality profile
#' * seasAmp = seasonality amplitude
#' * trAv = offset value of time series
#' * remSd = standard deviation of the remainder
#' * distMaglim = limits of the disturbance magnitude, should be a vector with the minimum and maximum value. If the minimum equals the maximum value, the disturbance magnitude is fixed for each simulated time series (and equal to the minimum value). When the minimum value does not equal the maximum value, a disturbance magnitude is randomly chosen in the given interval for each simulated time series.
#' * distTy = year of the disturbance. If distTy equals one, the disturbance will take place in the first year. The exact disturbance date (day or year) is randomly chosen per time series.
#' * distReclim = limits of the recovery duration of the recovery [number of observations], should be a vector with the minimum and maximum value. If the minimum equals the maximum value, the recovery period is fixed for each simulated time series (and equal to the minimum value). When the minimum value does not equal the maximum value, a recovery period is randomly chosen in the given interval for each simulated time series.</font>
#' * remcoef = list of ARMA models of remainder component. For each simulated time series, a remainder model is randomly chosen from the list.
#' * mval = number of missing values to be introduced. If mval equals NA, no missing values are introduced. For missing values with a random interval (see mvaldist), this should equal the fraction of missing values (mval equal to 0.1 will result in an NA value for 10% of the time series). For missing values having a regular interval, every mval observations one value is kept (eg for a daily time series, a mval equal to 5 will result in one observation every 5 days).
#' * mvaldist = the distribution of the missing values. Should equal 'random'or 'interval'.
#'
#' Function outputs:
#'
#' * a list with the simulated time series, offest, seasonality, remainder, disturbance component and parameters used for the simulation. The time series (components) are stored as matrix where each row is a time series and the columns are associated with the observation numbers.
#'
## ----simulCase-----------------------------------------------------------
simulCase <- function(nrep, nyr, nobsYr, nDr,seasAv,seasAmp,
                      trAv, remSd,distMaglim,distTy,distReclim, remcoef, mval, mvaldist){
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
    sts <- simulTS(nyr, nobsYr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remcoef[[modi]], distMag, distT, distRec)
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
                           m_distRec)
  names(TSsimParam) <- c("number_yrs", "obs_per_year", "number_droughts",
                         "seas_amp", "trend_av", 'rem_sd',
                         'dist_magn', 'dist_time', 'dist_rec')
  TSsimParam$year_drought <- m_year_dr
  TSsimParam$rem_coef <- m_remcoef
  list( TSsim, TSsimTr, TSsimSeas, TSsimRem, TSsimDist, TSsimParam)
}# function to simulate a case with n repetitions and store all settings

