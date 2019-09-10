#-------------------------------------------------
# Test recovery indicators
library(devtools)
library(rgrowth)
library(bfast)
library(signal)
library(foreach)
library(doParallel)

#-------------------------------------------------
# inputs
ifolder <- '/home/wanda/Dropbox/20190816_SimulationVar/SimulatedTS/'
ofolder <- '/home/wanda/Dropbox/20190816_SimulationVar/SimulatedTSStab/'
ecoreg <- 567
basename <- c('S1_Sample_EcoReg_','_noDist_Tree_50_lossYr_18_scl_30_npnt_1000_DESCENDING_10_H_IW_VH')
pol <- 'VH'
simcase <- 'len' #'distMag' distRec distT remSd seasAmp len  dr

#-------------------------------------------------
# functions

# Load RData with user specified name
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Frazier
# tsio: vector of observations (time series)
# tdist: observation number of disturbance 
# obspyr: number of observations in 1 year
calcFrazier <- function(tsio, tdist, obspyr){
  if((tdist>(1+2*obspyr)) & (tdist < (length(tsio)-obspyr-1))){
    Vpre <- mean(tsio[(tdist-(2*obspyr)):(tdist-1)], na.rm=T)
    V0 <- tsio[tdist]
    Ddist <- Vpre-V0
    ARI <- max(tsio[(tdist):(tdist+(1*obspyr))], na.rm=T) - V0
    RRI <- ARI/Ddist
    R80P <- max(tsio[(tdist):(tdist+(1*obspyr))], na.rm=T)/(Vpre*0.8)
    YrYr <- (mean(tsio[(tdist):(tdist+(1*obspyr))], na.rm=T)-V0)/0.5
    lst <- list(RRI, R80P, YrYr)
    names(lst) <- c('RRI', 'R80P', 'YrYr')
  }else{
    lst <- list(NA, NA, NA)
    names(lst) <- c('RRI', 'R80P', 'YrYr')
  }
  lst
}

# BFAST
# tsio: vector of observations (time series)
# obspyr: number of observations in one year
calcBFrec <- function(tsio, obspyr){
  tsi <- ts(tsio, frequency = obspyr)
  bf <- bfast(tsi, h=0.15, max.iter = 1, season='harmonic')#, hpc = 'foreach'
  
  if (bf$nobp$Vt ==F){
    # Frazier
    tr <- as.numeric(bf$output[[1]]$Tt)
    tbp <- 1+as.numeric(bf$output[[1]]$Vt.bp[1])
    if((tbp>(1+2*obspyr)) & (tbp < (length(tr)-obspyr-1))){
      frz <- calcFrazier(as.numeric(tr), tbp, floor(obspyr))
    }else{
      frz <- list(NA, NA, NA)
      names(frz) <- c('RRI', 'R80P', 'YrYr')
      }
    
    # Trends
    # bf$output[[1]]$Vt.bp
    if (length(bf$output[[1]]$Vt.bp)==1){
      bp1 <- as.numeric(bf$output[[1]]$Vt.bp[1])
      bp2 <- length(tsi)
      sl <- (tr[bp2-1] - tr[bp1+1])/(bp2-bp1-2)
    }
    if(length(bf$output[[1]]$Vt.bp)>1){
      bp1 <- as.numeric(bf$output[[1]]$Vt.bp[1])
      bp2 <- as.numeric(bf$output[[1]]$Vt.bp[2])
      sl <- (tr[bp2-1] - tr[bp1+1])/(bp2-bp1-2)
    }
    
  }else{
    frz <- list(NA, NA, NA)
    names(frz) <- c('RRI', 'R80P', 'YrYr')
    sl <- NA}
  list(sl, frz)
}

# BFAST with varying settings
calcBFrecSet <- function(tsio, obspyr){
  tsi <- ts(tsio, frequency = obspyr)
  hrange <- c(0.01,0.05,0.1,0.15,0.2)
  srange <- c('none', 'harmonic')#'dummy'
  bflist <- list()
  # loop over all setting combinations
  slf <- -1
  magnf <- 1
  for (hi in 1:length(hrange)){
    for (si in 1:length(srange)){
      bf <- tryCatch( {bfast(tsi, h=hrange[hi], max.iter = 1, season=srange[si])}, error = function(e){bf<- NA})
      if(as.logical(is.na(bf)[1])){next}
      if (bf$nobp$Vt ==F){
        magnf <- bf$Mags[1,3]
        bpt <- bf$output[[1]]$Vt.bp[1]
        trt <- bf$output[[1]]$Tt
        slf <- as.numeric(trt[bpt+2] - trt[bpt+1])
        print(slf)
        print(magnf)
        if ((slf>0) & (magnf<0)) break
      }
      if ((slf>0) & (magnf<0)) break
    }
    if ((slf>0) & (magnf<0)) break
  }
  # compute metrics
  if (bf$nobp$Vt ==F){
    # Frazier
    tr <- as.numeric(bf$output[[1]]$Tt)
    tbp <- 1+as.numeric(bf$output[[1]]$Vt.bp[1])
    if((tbp>(1+2*obspyr)) & (tbp < (length(tr)-obspyr-1))){
      frz <- calcFrazier(as.numeric(tr), tbp, floor(obspyr))
    }else{
      frz <- list(NA, NA, NA)
      names(frz) <- c('RRI', 'R80P', 'YrYr')
    }
    bp1 <- as.numeric(bf$output[[1]]$Vt.bp[1])
    sl <- (tr[bp1+2] - tr[bp1+1])
    
  }else{
    frz <- list(NA, NA, NA)
    names(frz) <- c('RRI', 'R80P', 'YrYr')
    sl <- NA}
  list(sl, frz)
}

# Rgrowth with known disturbance time
# tsio: vector of observations (time series)
# obspyr: number of observations per year
# starty: year of first observation
# dt: number of days between observations
# s: Rgrowth setting - post-regrowth stability criterion (min number of years stable after recovery)
# w: Rgrowth setting - minimum years after disturbance before recovery
calcRgrowth <- function(tsio, obspyr, starty, dt, s, w){
  ts.Date <- as.Date(paste0(starty, '-01-01')) + seq(1, (dt*length(tsio)),dt) - 1
  ts.zoo <- zoo(tsio, ts.Date)
  ts.TS <- ts(tsio, start=c(starty,1), frequency=obspyr)
  bfm <- bfastmonitor(ts.TS, start = c((starty+1),1), formula = response~harmon, order = 1, plot = F)
  
  out <- tryCatch(
    {reg <- tsreg(ts.zoo, change = bfm$breakpoint, 
                   startOffset = "floor",
                   h = 0.5, 
                   history='all', 
                   s=s, #post-regrowth stability criterion (min number of years stable after recovery)
                   w=w, # minimum years after disturbance before recovery 
                   plot=F)
      # time disturbance
      m_tdst <-reg$disturbance
      # time recovery
      m_trec <- reg$regrowth_onset
      # duration recovery
      m_lrec <- m_trec - m_tdst
      
      list(m_tdst, m_trec, m_lrec)
    },
    error = function(e){
      m_tdst <-NA
      # time recovery
      m_trec <- NA
      # duration recovery
      m_lrec <- NA
      
      list(m_tdst, m_trec, m_lrec)
    }
  )
  out
}

# Rgrowth with unknown disturbance time
# tsio: vector of observations (time series)
# obspyr: number of observations per year
# starty: year of first observation
# tdist: disturbance timing (given in observation number)
# dt: number of days between observations
# s: Rgrowth setting - post-regrowth stability criterion (min number of years stable after recovery)
# w: Rgrowth setting - minimum years after disturbance before recovery
calcRgrowthRec <- function(tsio, obspyr, starty, tdist, dt, s, w){
  ts.Date <- as.Date(paste0(starty, '-01-01')) + seq(1, (dt*length(tsio)),dt) - 1
  ts.zoo <- zoo(tsio, ts.Date)
  ts.TS <- ts(tsio, start=c(starty,1), frequency=obspyr)
  # bfm <- bfastmonitor(ts.TS, start = c((starty+1),1), formula = response~harmon, order = 1, plot = F)
  # nobs <- obspyr*dt
  out <- tryCatch(
    {reg <- tsreg(ts.zoo, change = starty+(tdist/nobspyr), #(tdist/nobs)
                  startOffset = "floor",
                  h = 0.5, 
                  history='all', 
                  s=s, #post-regrowth stability criterion (min number of years stable after recovery)
                  w=w, # minimum years after disturbance before recovery 
                  plot=F)
    # time disturbance
    m_tdst <-reg$disturbance
    # time recovery
    m_trec <- reg$regrowth_onset
    # duration recovery
    m_lrec <- m_trec - m_tdst
    
    m_lrec
    },
    error = function(e){
      # duration recovery
      m_lrec <- NA
      
       m_lrec
    }
  )
  out
}

#-------------------------------------------------
# import time series
simTS <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTS_', simcase, '.rda')))
simTSTr <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTSTr_', simcase, '.rda')))
simTSSeas <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTSSeas_', simcase, '.rda')))
simTSDist <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTSDist_', simcase, '.rda')))
simTSRem <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTSRem_', simcase, '.rda')))
simTSParam <- loadRData(file = file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTSParam_', simcase, '.rda')))

#-------------------------------------------------
# measure recovery - Rgrowth
m_tdst <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_trec <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_lrec <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_lrecRec <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_tdst <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_trec <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_lrec <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_lrecRec <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))

for (i in 1:length(simTS)){# parameter
  print(i)
  for (ii in 1:dim(simTS[[i]])[1]){# ts
    dt <- 6
    tsi <- simTS[[i]][ii, seq(1,dim(simTS[[i]])[2],dt)]
    nobspyr <- 365/dt
    tdisti <- simTSParam[[i]][ii,]$dist_time
    tdistcorri <- tdisti/dt
    treci <- simTSParam[[i]][ii,]$dist_rec
    treccorri <- treci/dt
    
    rec_rgrowth <- calcRgrowth(tsi, nobspyr, 2014, dt, 0.5, 0.5)
    m_tdst[ii,i] <- rec_rgrowth[[1]] # timing disturbance
    m_trec[ii,i] <- rec_rgrowth[[2]] # timing recovery 
    m_lrec[ii,i] <- rec_rgrowth[[3]] # recovery period
    s_tdst[ii,i] <- 2014 + (tdisti/365) # true timing disturbance
    s_trec[ii,i] <- 2014 + ((tdisti + treci)/365) # true timing recovery 
    s_lrec[ii,i] <- s_trec[ii,i] - s_tdst[ii,i] # true recovery period
    
    m_lrecRec[ii,i] <- calcRgrowthRec(tsi, nobspyr, 2014,tdistcorri, dt, 0.5, 0.5)# recovery period
  }
}
meas_Rgrowth <- list(m_tdst, m_trec, m_lrec, s_tdst, s_trec, s_lrec, m_lrecRec)
names(meas_Rgrowth) <- c("m_tdst", "m_trec", "m_lrec", "s_tdst", "s_trec", "s_lrec", "m_lrecRec")

#-------------------------------------------------
# measure recovery  
m_RRI <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_R80p <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_YrYr <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_RRIsm <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_R80psm <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_YrYrsm <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_RRI <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_R80p <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_YrYr <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
for (i in 1:length(simTS)){# parameter
  print(i)
  for (ii in 1:dim(simTS[[i]])[1]){# ts
    dt <- 6
    tdist <- simTSParam[[i]]$dist_time[ii]
    tdistcorr <- 1+ceiling(tdist/dt)*dt
    tsref <- simTSTr[[i]][ii, ] + simTSDist[[i]][ii, ]
    tsi <- simTS[[i]][ii, ]
    
    # smooth time series
    tsf <- tsi
    tsf[is.na(tsf)==F] <- sgolayfilt(tsf[seq(1,dim(simTS[[i]])[2],dt)], p=2, n=11)
    rec_rgrowthsm <- calcFrazier(tsf, tdistcorr, 365)
    m_RRIsm[ii,i] <- rec_rgrowthsm[[1]]
    m_R80psm[ii,i] <- rec_rgrowthsm[[2]]
    m_YrYrsm[ii,i] <- rec_rgrowthsm[[3]]
    
    rec_rgrowth <- calcFrazier(tsi, tdistcorr, 365)
    m_RRI[ii,i] <- rec_rgrowth[[1]]
    m_R80p[ii,i] <- rec_rgrowth[[2]]
    m_YrYr[ii,i] <- rec_rgrowth[[3]]
    
    rec_rgrowthref <- calcFrazier(tsref, tdist, 365)
    s_RRI[ii,i] <- rec_rgrowthref[[1]]
    s_R80p[ii,i] <- rec_rgrowthref[[2]]
    s_YrYr[ii,i] <- rec_rgrowthref[[3]]
  }
}
meas_Frazier <- list(m_RRI, m_R80p, m_YrYr, m_RRIsm, m_R80psm, m_YrYrsm, s_RRI, s_R80p, s_YrYr)
names(meas_Frazier) <- c("m_RRI", "m_R80p", "m_YrYr", "m_RRIsm", "m_R80psm", "m_YrYrsm", "s_RRI", "s_R80p", "s_YrYr")

#-------------------------------------------------
# BFAST slope + Frazier
m_SLbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_RRIbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_R80pbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_YrYrbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_SLbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_RRIbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_R80pbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_YrYrbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))

for (i in 1:length(simTS)){# parameter
  print(i)
  cl<-makeCluster(11)
  registerDoParallel(cl)
  
  tmpout <- foreach(ii=1:dim(simTS[[i]])[1], .packages = "bfast") %dopar% {
    tdist <- simTSParam[[i]]$dist_time[ii]
    tdistcorr <- 1+ceiling(tdist/6)
    tsio <- simTS[[i]][ii, seq(1,dim(simTS[[i]])[2],6)]
    
    bfrec <- calcBFrec(tsio, 365/6)
    tsioref <- simTSDist[[i]][ii, seq(1,dim(simTS[[i]])[2],6)] + simTSTr[[i]][ii, seq(1,dim(simTS[[i]])[2],6)]
    bfrecref <- calcFrazier(tsioref, tdistcorr, 365/6)
  
    outp <- list(bfrec[[1]], bfrec[[2]]$RRI, bfrec[[2]]$R80P,bfrec[[2]]$YrYr,
                 bfrecref$RRI, bfrecref$R80P, bfrecref$YrYr, tsioref[tdistcorr+2] - tsioref[tdistcorr+1])
    names(outp) <- c('m_SLbf', 'm_RRIbf', 'm_R80pbf', 'm_YrYrbf',
                     's_RRIbf', 's_R80pbf', 's_YrYrbf', 's_SLbf')
    outp
  }
  stopCluster(cl)
  
  m_SLbf[,i] <- unlist(sapply(tmpout, '[', 'm_SLbf')) #bfrec[[1]]
  m_RRIbf[,i] <- unlist(sapply(tmpout, '[', 'm_RRIbf')) #bfrec[[2]]$RRI
  m_R80pbf[,i] <- unlist(sapply(tmpout, '[', 'm_R80pbf')) #bfrec[[2]]$R80P
  m_YrYrbf[,i] <- unlist(sapply(tmpout, '[', 'm_YrYrbf')) #bfrec[[2]]$YrYr
  s_RRIbf[,i] <- unlist(sapply(tmpout, '[', 's_RRIbf')) #bfrecref$RRI
  s_R80pbf[,i] <- unlist(sapply(tmpout, '[', 's_R80pbf')) #bfrecref$R80P
  s_YrYrbf[,i] <- unlist(sapply(tmpout, '[', 's_YrYrbf')) #bfrecref$YrYr
  s_SLbf[,i] <- unlist(sapply(tmpout, '[', 's_SLbf')) #tsioref[tdistcorr+2] - tsioref[tdistcorr+1]
}

meas_bfast <- list(m_SLbf, m_RRIbf, m_R80pbf, m_YrYrbf, s_RRIbf, s_R80pbf, s_YrYrbf, s_SLbf)
names(meas_bfast) <- c("m_SLbf", "m_RRIbf", "m_R80pbf", "m_YrYrbf", "s_RRIbf", "s_R80pbf", "s_YrYrbf", "s_SLbf")

#-------------------------------------------------
# BFAST slope + Frazier with optimised settings
m_SLbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_RRIbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_R80pbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
m_YrYrbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_SLbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_RRIbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_R80pbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))
s_YrYrbf <- matrix(NA, nrow = dim(simTS[[1]])[1], ncol = length(simTS))

for (i in 1:length(simTS)){# parameter
  print(i)
  cl<-makeCluster(11)
  registerDoParallel(cl)
  
  tmpout <- foreach(ii=1:dim(simTS[[i]])[1], .packages = "bfast") %dopar% {
    tdist <- simTSParam[[i]]$dist_time[ii]
    tdistcorr <- 1+ceiling(tdist/6)
    tsio <- simTS[[i]][ii, seq(1,dim(simTS[[i]])[2],6)]
    
    bfrec <- calcBFrecSet(tsio, 365/6)
    tsioref <- simTSDist[[i]][ii, seq(1,dim(simTS[[i]])[2],6)] + simTSTr[[i]][ii, seq(1,dim(simTS[[i]])[2],6)]
    bfrecref <- calcFrazier(tsioref, tdistcorr, 365/6)
    
    outp <- list(bfrec[[1]], bfrec[[2]]$RRI, bfrec[[2]]$R80P,bfrec[[2]]$YrYr,
                 bfrecref$RRI, bfrecref$R80P, bfrecref$YrYr, tsioref[tdistcorr+2] - tsioref[tdistcorr+1])
    names(outp) <- c('m_SLbf', 'm_RRIbf', 'm_R80pbf', 'm_YrYrbf',
                     's_RRIbf', 's_R80pbf', 's_YrYrbf', 's_SLbf')
    outp
  }
  stopCluster(cl)
  
  m_SLbf[,i] <- unlist(sapply(tmpout, '[', 'm_SLbf')) #bfrec[[1]]
  m_RRIbf[,i] <- unlist(sapply(tmpout, '[', 'm_RRIbf')) #bfrec[[2]]$RRI
  m_R80pbf[,i] <- unlist(sapply(tmpout, '[', 'm_R80pbf')) #bfrec[[2]]$R80P
  m_YrYrbf[,i] <- unlist(sapply(tmpout, '[', 'm_YrYrbf')) #bfrec[[2]]$YrYr
  s_RRIbf[,i] <- unlist(sapply(tmpout, '[', 's_RRIbf')) #bfrecref$RRI
  s_R80pbf[,i] <- unlist(sapply(tmpout, '[', 's_R80pbf')) #bfrecref$R80P
  s_YrYrbf[,i] <- unlist(sapply(tmpout, '[', 's_YrYrbf')) #bfrecref$YrYr
  s_SLbf[,i] <- unlist(sapply(tmpout, '[', 's_SLbf')) #tsioref[tdistcorr+2] - tsioref[tdistcorr+1]
}

meas_bfast_optim <- list(m_SLbf, m_RRIbf, m_R80pbf, m_YrYrbf, s_RRIbf, s_R80pbf, s_YrYrbf, s_SLbf)
names(meas_bfast_optim) <- c("m_SLbf", "m_RRIbf", "m_R80pbf", "m_YrYrbf", "s_RRIbf", "s_R80pbf", "s_YrYrbf", "s_SLbf")

#-------------------------------------------------
# export
save(meas_bfast_optim, file=file.path(ofolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASbfastOptim_', simcase, '.rda')))
save(meas_bfast, file=file.path(ofolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASbfast_', simcase, '.rda')))
save(meas_Frazier, file=file.path(ofolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASfrazier_', simcase, '.rda')))
save(meas_Rgrowth, file=file.path(ofolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASrgrowth_', simcase, '.rda')))


