#' ---
#' title: "Functions to evaluate the performance"
#' author: "Wanda De Keersmaecker"
#' date: "2/12/2020"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' ## R squared
## ------------------------------------------------------------------------
rsq <- function(x, y) {
  if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3)){
    rs <- summary(lm(y~x))$r.squared
  }else(rs <- NA)
  rs
}

#' 
#' ## RMSE
## ------------------------------------------------------------------------
# calculate RMSE
rmse <- function(val, meas){
  sqrt(mean((val - meas)^2,na.rm=TRUE))
} 

#' 
#' ## MAPE
## ------------------------------------------------------------------------
# calculate MAPE
mape <- function(val, meas){
  mean(abs(val - meas)/abs(val),na.rm=TRUE)
} 

#' 
#' ## Calculate performance of recovery indicators
## ------------------------------------------------------------------------
# Derive performance (R2 or RMSE) from recovery indicators                   
calcPerf <- function(val, meas, sttngs, recSttngs, metr, perf){
  lst <- list()
  simcases <- names(meas)
  #vsimcases <- paste0('V',simcases)
  
  for(sci in 1:length(meas)){
    vls <- list()
    for(rpi in 1:length(meas[[1]])){
      if((metr == 'SL') & ((recSttngs$input[rpi] == 'raw') | (recSttngs$input[rpi] == 'smooth'))){                
      }else{
        val[[sci]][[rpi]][is.infinite(val[[sci]][[rpi]])] <- NA
        meas[[sci]][[rpi]][is.infinite(meas[[sci]][[rpi]])]<-NA
        if(perf == 'R2'){
          vl <- sapply(1:dim(meas[[sci]][[rpi]])[2], function(ii) rsq(val[[sci]][[rpi]][,ii], meas[[sci]][[rpi]][,ii]))    
        }
        if(perf == 'RMSE'){
          vl <- sapply(1:dim(meas[[sci]][[rpi]])[2], function(ii) rmse(val[[sci]][[rpi]][,ii], meas[[sci]][[rpi]][,ii]))    
        }
        if(perf == 'MAPE'){
          vl <- sapply(1:dim(meas[[sci]][[rpi]])[2], function(ii) mape(val[[sci]][[rpi]][,ii], meas[[sci]][[rpi]][,ii]))    
        }
        tmp <- data.frame(t(vl))
        #print(typeof(tmp))
        names(tmp) <- sttngs[[simcases[sci]]]
        tmp2 <- melt(tmp)
        tmp2$Metric <- factor(metr)
        tmp2$Dense <- factor(recSttngs$freq[rpi])
        tmp2$Smooth <- revalue(factor(recSttngs$input[rpi]), c("BFAST"="segmented", 'smooth'='smoothed'))
        tmp2$Period <- revalue(factor(recSttngs$nDist[rpi]), c("1"="Short", "12"="Long"))
        tmp2$Period[tmp2$Dense == 'annual'] <- 'Long'
        #tmp2$Method <- factor(paste0(metr, ', ',recSttngs$freq[rpi], ', ', recSttngs$input[rpi],', ', recSttngs$nDist[rpi]))
        #print(typeof(tmp2))
        vls <- rbind(vls,tmp2)
      }
    }
    lst[[names(meas)[sci]]] <- vls
  }
  lst
}

