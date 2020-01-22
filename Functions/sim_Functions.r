loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

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

# ARMA coefficients
getARMAcoef <- function(tsx){
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
    coefmod <- c(coefmod, ma = as.numeric(arm$coef[(1+teller):(teller+arm$arma[1])]))
    teller <- teller + arm$arma[2]
  }
  coefmod <- c(coefmod, order = as.numeric(c(arm$arma[1], arm$arma[6], arm$arma[2])))
  coefmod
}


simulTS <- function(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remMod, distMag, distT, distRec){
  #-------------------------------------------------
  # simulate seasonality
  simSeas <- rep(as.numeric(seasAv[1:(nobsyr*2)]), times=ceiling(nyr/2))
  if (ceiling(nyr/2) != (nyr/2)){
    simSeas <- simSeas[1:(length(simSeas)-nobsyr)]
  }
  simSeas <- simSeas/max(seasAv)*seasAmp
    
  # introduce drought years (seasonality equal to minimum value)
  if (nDr==0) ydr <- 0
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
    
  #-------------------------------------------------
  # simulate disturbance
  if ((distT+distRec) <= (nobsyr*nyr)){
    simDist <- c(rep(0,(distT-1)), seq(distMag, 0, length.out = distRec), rep(0, ((nobsyr*nyr)-distRec -distT +1)))
  }
  if ((distT+distRec) > (nobsyr*nyr)){
    simDist <- c(rep(0,(distT-1)), seq(distMag, 0, length.out = distRec))
    simDist <- simDist[1:(nobsyr*nyr)]
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



