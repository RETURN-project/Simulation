#' R squared
#'
#' @param x vector of x values
#' @param y vector of y values
#'
#' @return the R squared between the x and y variables
#' @export
rsq <- function(x, y) {
  if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3)){
    rs <- summary(lm(y~x))$r.squared
  }else(rs <- NA)
  rs
}

#' RMSE
#'
#' @param val vector of x values
#' @param meas vector of y values
#'
#' @return the RMSE of the two vectors
#' @export
rmse <- function(val, meas){
  sqrt(mean((val - meas)^2,na.rm=TRUE))
}

#' MAPE
#'
#' @param val vector of x values
#' @param meas vector of y values
#'
#' @return the MAPE of the two vectors
#' @export
mape <- function(val, meas){
  mean(abs(val - meas)/abs(val),na.rm=TRUE)
}

#' Derive performance (R2, MAPE or RMSE) from recovery indicators
#'
#' @param val  ground truth
#' @param meas measured
#' @param sttngs simulation settings file
#' @param recSttngs recovery settings file
#' @param metr recovery metric being evaluated
#' @param perf performance indicator
#'
#' @return list of performance indicator
#' @export
#' @import reshape2
#' @import plyr
calcPerf <- function(val, meas, sttngs, recSttngs, metr, perf){
  lst <- list()
  simcases <- names(meas)

  for(sci in 1:length(meas)){# evaluated parameters
    vls <- list()
    for(rpi in 1:length(meas[[1]])){# recovery settings
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
        names(tmp) <- sttngs[[simcases[sci]]]
        tmp2 <- melt(tmp)
        tmp2$Metric <- factor(metr)
        tmp2$Dense <- factor(recSttngs$freq[rpi])
        tmp2$Smooth <- factor(recSttngs$input[rpi])#revalue(factor(recSttngs$input[rpi]), c("BFAST"="segmented", 'smooth'='smoothed'))
        tmp2$Period <- factor(recSttngs$nDist[rpi])#revalue(factor(recSttngs$nDist[rpi]), c("1"="Short", "12"="Long"))
        vls <- rbind(vls,tmp2)
      }
    }
    vls$Period <- revalue(vls$Period, c("1"="Short", "12"="Long"))

    levels(vls$Period) <- c('Short', 'Long')
    vls$Period[vls$Dense == 'annual'] <- "Long"
    lst[[names(meas)[sci]]] <- vls
  }
  lst
}

