#-------------------------------------------------
# Simulate time series

library(gdata)
#-------------------------------------------------
# inputs
ifolder <- '/home/wanda/Desktop/Data/RETURN/20190801_SimulationFixed/CharData/'
ofolder <- '/home/wanda/Desktop/Data/RETURN/20190816_SimulationVar/SimulatedTS/'

basename <- c('S1_Sample_EcoReg_','_noDist_Tree_50_lossYr_18_scl_30_npnt_1000_DESCENDING_10_H_IW_VH')
regions <- c(567)
pol <- 'VH'

#-------------------------------------------------
# functions 
simulTS <- function(nyr, nobsyr, tMiss, nDr, seasAv, seasAmp, trAv, remSd, remMod, distMag, distT, distRec){
  #-------------------------------------------------
  # simulate seasonality
  simSeas <- rep(as.numeric(seasAv[1:(nobsyr*2)]), times=ceiling(nyr/2))
  if (ceiling(nyr/2) != (nyr/2)){
    simSeas <- simSeas[1:(length(simSeas)-365)]
  }
  simSeas <- simSeas/max(seasAv)*seasAmp
  # introduce drought years (without seasonality)
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
  # set missing values
  if (is.na(tMiss[1])==F){
    simSeas[tMiss] <- NA
  }
  #-------------------------------------------------
  # simulate trend
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
  
  list(ydr, ts(t(rbind(simSeas, simTr, simRem, simDist, simSeas+simTr+simRem+simDist)), frequency = nobsyr))
}
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
simulCase <- function(nrep, nyr, nobsYr, nDr,seasAv,seasAmp,
                      trAv, remSd,distMaglim,distTy,distReclim, remcoef){
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
  # missing values to be introduced
  tMiss <- 1:(nyr*nobsYr) # missing values
  tMiss <- tMiss[-seq(1,(nyr*nobsYr),6)]
  
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
    distT <- ((distTy-1)*365) + sample((1:365), 1) 
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

#-------------------------------------------------
# import data
for (ei in length(regions)){
  ecoreg <- regions[ei]
  keep(ifolder, ofolder, basename, regions, pol, ecoreg,simulTS, loadRData, simulCase, ei, sure=T)
  
  seasmax <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_seas', pol, 'max.rda')))
  seasmean <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_seas', pol, 'mean.rda')))
  trmean <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_Tr', pol, 'mean.rda')))
  remsd <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_Rem_', pol, 'sd.rda')))
  remcoef <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_Rem_', pol, 'coef.rda')))
  
  #-------------------------------------------------
  # standard settings and subsequently let each variable vary:
  STnrep <- 1000 
  STnyr <- 6
  STnobsYr <- 365
  STnDr <- 0 
  STseasAv <- seasmean
  STseasAmp <- quantile(seasmax, .5)
  STtrAv <- mean(trmean)
  STremSd <- quantile(remsd, .5)
  STdistMaglim <- -c(1,9) # magnitude between 1 and 9 dB
  STdistT <- 3 # disturbance in year 3
  STdistReclim <- c(0.5,1.5)*365 # recovery time between 0.5 and 1.5 years
  
  VnDr <- c(0,1) # number of droughts
  Vnyr <- seq(3,15,3)# number of years
  VseasAmp <- quantile(seasmax, c(.05, .275, .5, .725, .95))# seasonal amplitude
  VremSd <- quantile(remsd, c(.05, .275, .5, .725, .95))# sd of remainder
  VdistMag <- -seq(1,9,2)# magnitude disturbance
  VdistT <- seq(1,5,1)# timing disturbance
  VdistRec <- seq(0.5,1.5,0.5)*365# recovery disturbance
  
  sttngs <- list(list(STnrep,
                      STnyr,
                      STnobsYr,
                      STnDr, 
                      STseasAv, 
                      STseasAmp,
                      STtrAv, 
                      STremSd, 
                      STdistMaglim,STdistT, STdistReclim), 
                 list(VnDr,  
                      Vnyr,  
                      VseasAmp, 
                      VremSd,  
                      VdistMag, 
                      VdistT,  
                      VdistRec))
  names(sttngs) <- c('standard', 'variable')
  names(sttngs[[1]]) <- c('STnrep',
  'STnyr',
  'STnobsYr',
  'STnDr', 
  'STseasAv', 
  'STseasAmp',
  'STtrAv', 
  'STremSd', 
  'STdistMaglim','STdistT', 'STdistReclim')
  names(sttngs[[2]]) <- c('VnDr',  
                          'Vnyr',  
                          'VseasAmp', 
                          'VremSd',  
                          'VdistMag', 
                          'VdistT',  
                          'VdistRec')
  save(sttngs, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_settings.rda')))
  
  #-------------------------------------------------
  # Simulate time series - drought
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(VnDr)){
    sc <- simulCase(STnrep, STnyr, STnobsYr, VnDr[i],STseasAv,STseasAmp,
              STtrAv, STremSd,STdistMaglim,STdistT,STdistReclim, remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
    rm(sc)
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_dr.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_dr.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_dr.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_dr.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_dr.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_dr.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
  #-------------------------------------------------
  # Simulate time series - length time series
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(Vnyr)){
    sc <- simulCase(STnrep, Vnyr[i], STnobsYr, STnDr,STseasAv,STseasAmp,
                    STtrAv, STremSd,STdistMaglim,STdistT,STdistReclim, remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_len.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_len.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_len.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_len.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_len.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_len.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
  #-------------------------------------------------
  # Simulate time series - seasonal amplitude
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(VseasAmp)){
    sc <- simulCase(STnrep, STnyr, STnobsYr, STnDr,STseasAv,VseasAmp[i],
                    STtrAv, STremSd,STdistMaglim,STdistT,STdistReclim, remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_seasAmp.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_seasAmp.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_seasAmp.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_seasAmp.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_seasAmp.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_seasAmp.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
  #-------------------------------------------------
  # Simulate time series - remainder SD
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(VremSd)){
    sc <- simulCase(STnrep, STnyr, STnobsYr, STnDr,STseasAv,STseasAmp,
                    STtrAv, VremSd[i],STdistMaglim,STdistT,STdistReclim, remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_remSd.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_remSd.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_remSd.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_remSd.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_remSd.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_remSd.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
  #-------------------------------------------------
  # Simulate time series - dist magnitude
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(VdistMag)){
    sc <- simulCase(STnrep, STnyr, STnobsYr, STnDr,STseasAv,STseasAmp,
                    STtrAv, STremSd,c(VdistMag[i],VdistMag[i]),STdistT,STdistReclim, remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_distMag.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_distMag.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_distMag.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_distMag.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_distMag.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_distMag.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
  #-------------------------------------------------
  # Simulate time series - dist timing
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(VdistT)){
    sc <- simulCase(STnrep, STnyr, STnobsYr, STnDr,STseasAv,STseasAmp,
                    STtrAv, STremSd,STdistMaglim,VdistT[i],STdistReclim, remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_distT.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_distT.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_distT.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_distT.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_distT.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_distT.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
  #-------------------------------------------------
  # Simulate time series - dist recovery
  simTS <- list()
  simTSTr <- list()
  simTSSeas <- list()
  simTSRem <- list()
  simTSDist <- list()
  simTSParam <- list()
  
  for (i in 1:length(VdistRec)){
    sc <- simulCase(STnrep, STnyr, STnobsYr, STnDr,STseasAv,STseasAmp,
                    STtrAv, STremSd,STdistMaglim,STdistT,c(VdistRec[i],VdistRec[i]), remcoef)
    simTS[[i]] <- sc[[1]]
    simTSTr[[i]] <- sc[[2]]
    simTSSeas[[i]] <- sc[[3]]
    simTSRem[[i]] <- sc[[4]]
    simTSDist[[i]] <- sc[[5]]
    simTSParam[[i]] <- sc[[6]]
  }
  save(simTS, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTS_distRec.rda')))
  save(simTSTr, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSTr_distRec.rda')))
  save(simTSSeas, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSSeas_distRec.rda')))
  save(simTSRem, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSRem_distRec.rda')))
  save(simTSDist, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSDist_distRec.rda')))
  save(simTSParam, file = file.path(ofolder, paste0(basename[1], ecoreg,basename[2], '_', pol, '_simTSParam_distRec.rda')))
  rm(simTS,simTSTr,simTSSeas, simTSRem, simTSDist, simTSParam)
  
}

# plot(ts(simTS[[1]][1,seq(1,365*6,6)], frequency=365/6))
# ti <- 1
# plot(ts(simTS[[ti]][1,seq(1,(Vnyr[ti]*365),6)], frequency = 365/6), ylab = 'Simulated time series')
# lines(ts(simTS[[ti]][10,seq(1,(Vnyr[ti]*365),6)], frequency = 365/6), col='red')
# lines(ts(simTS[[ti]][100,seq(1,(Vnyr[ti]*365),6)], frequency = 365/6), col='green')

# ti <- 1
# plot(ts(simTS[[ti]][1,seq(1,(6*365),6)], frequency = 365/6), ylab = 'Simulated time series')
# lines(ts(simTS[[ti]][10,seq(1,(6*365),6)], frequency = 365/6), col='red')
# lines(ts(simTS[[ti]][100,seq(1,(6*365),6)], frequency = 365/6), col='green')
