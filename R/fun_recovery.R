#' Calculate recovery metrics from a time series with known disturbance date. The calcFrazier function derives the RRI, R80P and YrYr recovery indicators, defined by Frazier et al. (2018). The indicators are originally developped for annual long-term time series of optical vegetation indices (the indicators are shown in the figures below). Yet, in order to be able to derive the indicators as well for dense and/or short time series, a modified version is suggested. Here, the user can define the time period before, during and after the disturbance that is used to derive the indicators. To reduce the interference of the seasonal pattern of dense time series, the chosen time period should cover blocks of n years. Moreover, given the potentially high noise levels of dense time series, the mean value instead of the maximum value was used in the formulas. (Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018). Analyzing spatial and temporal variability in short-term rates of post-fire vegetation return from Landsat time series. Remote Sensing of Environment, 205, 32-45.)
#'
#' @param tsio vector of observations (time series with a fixed observation frequency)
#' @param tdist observation number of disturbance, indicating the timing of the disturbance
#' @param obspyr number of observations per year
#' @param shortDenseTS TRUE or FALSE. In case TRUE, the metrics are adjusted to be compatible with short, dense time series
#' @param nPre If shortDenseTS is TRUE, number of years prior to the disturbance used to calculate the pre-disturbance value
#' @param nDist If shortDenseTS is TRUE, number of months used to quantify the time series value during the disturbance
#' @param nPostMin If shortDenseTS is TRUE,  the post-disturbance condition is quantified starting from nPostMin years after the disturbance
#' @param nPostMax If shortDenseTS is TRUE, max number of years after the disturbance used to quantify the post-disturbance condition
#'
#' @return a list containing the RRI recovery indicator, R80p recovery indicator and YrYr recovery indicator
#' @export
calcFrazier <- function(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax){
  #library(strucchange)
  #library(bfast)
    if(shortDenseTS){# Metrics adjusted for short, dense time series
        # check if there are enough observations before and after the disturbance to calculate the metrics
        if((tdist>((nPre*obspyr))) & (tdist < (length(tsio)-(nPostMax*obspyr)+1))){
            # Vpre = pre-disturbance value, mean of observations within nPre year period prior to disturbance
            Vpre <- mean(tsio[(tdist-(nPre*obspyr)):(tdist-1)], na.rm=T)
            # V0 =  value during disturbance (over a period of one month)
            V0 <- mean(tsio[tdist:(tdist+ (nDist*round(obspyr/12))-1)], na.rm=T)
            # Ddist =  decrease due to disturbance (~impact)
            Ddist <- Vpre-V0
            # Post-disturbance value
            Vpost <- mean(tsio[(tdist+(nPostMin*obspyr)):(tdist+(nPostMax*obspyr))], na.rm=T)
            # ARI: difference between maximum value within nPost years after disturbance and the disturbance value
            ARI <- Vpost - V0
            # RRI: Relative Recovery Index (~recovery relative to impact)
            RRI <- ARI/Ddist
            # R80p recovery index (~ ability to reach 80% of pre-disturbance value)
            R80P <- Vpost/(Vpre*0.8)
            # YrYR recovery index (~ related to slope)
            YrYr <- (Vpost-V0)/((nPostMax+nPostMin)/2)
            # make list of recovery indicators as output of the function
            lst <- list(RRI, R80P, YrYr)
            names(lst) <- c('RRI', 'R80P', 'YrYr')
            # give NA as output if not able to calculate the recovery indicatores
        }else{
            lst <- list(NA, NA, NA)
            names(lst) <- c('RRI', 'R80P', 'YrYr')
        }
    }else{#original metrics, typically applied on long time series with annual observations
        # check if there are enough observations before and after the disturbance to calculate the metrics
        if((tdist>((2*obspyr))) & (tdist < (length(tsio)-(5*obspyr)+1))){
            # Vpre = pre-disturbance value, mean of observations within nPre year period prior to disturbance
            Vpre <- mean(tsio[(tdist-(2*obspyr)):(tdist-1)], na.rm=T)
            # V0 =  value during disturbance (over a period of one month)
            V0 <- tsio[tdist]
            # Ddist =  decrease due to disturbance (~impact)
            Ddist <- Vpre-V0
            # ARI: difference between maximum value within nPost years after disturbance and the disturbance value
            ARI <- max(tsio[(tdist +(4*obspyr)):(tdist+(5*obspyr))], na.rm=T) - V0
            # RRI: Relative Recovery Index (~recovery relative to impact)
            RRI <- ARI/Ddist
            if(is.infinite(RRI)){RRI <- NA}
            # R80p recovery index (~ ability to reach 80% of pre-disturbance value)
            R80P <- max(tsio[(tdist +(4*obspyr)):(tdist+(5*obspyr))], na.rm=T)/(Vpre*0.8)
            if(is.infinite(R80P)){R80P <- NA}
            # YrYR recovery index (~ related to slope)
            YrYr <- (tsio[tdist+(5*obspyr)]-V0)/5
            # make list of recovery indicators as output of the function
            lst <- list(RRI, R80P, YrYr)
            names(lst) <- c('RRI', 'R80P', 'YrYr')
        }else{
            lst <- list(NA, NA, NA)
            names(lst) <- c('RRI', 'R80P', 'YrYr')
        }
    }
    lst
}

#' Post-disturbance slope and recovery metrics derived from BFAST0n trend segments. The calcBFASTrec function derives a set of recovery indicators after fitting a segmented trend in the time series. Using the breakpoints function of the strucchange package, a segmented trend is fitted (hereafter called BFAST0n trend segments). The detected break showing the largest change (in absolute values) is assumed to represent the disturbance. Using the segmented trend and detected disturbance date, the RRI, R80p, YrYr and the slope of the post-disturbance trend segment are derived as recovery indicators.
#'
#' @param tsio vector of observations (time series)
#' @param obspyr number of observations in one year
#' @param h This parameter defines the minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment.
#' @param shortDenseTS TRUE or FALSE. In case TRUE, the metrics are adjusted to be compatible with short, dense time series. In case FALSE, the input time series is assumed to have annual observations and at least 2 and 5 pre- and post-disturbance years, respectively.
#' @param nPre If shortDenseTS is TRUE, number of years prior to the disturbance used to calculate the pre-disturbance value
#' @param nDist If shortDenseTS is TRUE, number of months used to quantify the time series value during the disturbance
#' @param nPostMin If shortDenseTS is TRUE, min number of years after the disturbance used to quantify the recovery
#' @param nPostMax If shortDenseTS is TRUE, max number of years after the disturbance used to quantify the recovery
#'
#' @return a list containing  the RRI, R80p, YrYr recovery indicator derived from the BFAST0n trend segments and slope of the trend segment after the disturbance (sl).
#' @export
#' @import strucchange
#' @import stats
calcBFASTrec <- function(tsio, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax){
  # Create time series object, needed as input for BFAST
  tsi <- ts(tsio, frequency = obspyr)
  # Convert the time series object into a dataframe, needed for the breakpoints function
    datapp <- bfastpp(tsi, order = 1, lag = NULL, slag = NULL,
                  na.action = na.omit, stl = 'none')
    # Apply BFAST0n on time series: find breaks in the regression
    bp <- breakpoints(response ~ trend, data = datapp, h = h)##, breaks = nbrks
    # Check if BFAST0n found breakpoints
    if(is.na(bp$breakpoints[1])){# no breakpoint found
        tr <- fitted(bp, 0)
        sl <- (tr[2] - tr[1])
        frz <- list(NA, NA, NA, sl)
        names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
    }else{# at least one breakpoint found
        # Extract BFAST trend component and breaks
        cf <- coef(bp)
        # Extract BFAST trend component and breaks
        tbp <- bp$breakpoints #observation number of break
        #tr <- rep(NA,length(tsi))
        indna <- which(is.na(tsi)==F)
        tbp <- indna[tbp]   # correct observation number for missing values
        #tr[is.na(tsi)==F] <- fitted(bptst, length(tbptst))
        #Derive trend component without missing values
        bpf <- c(0, tbp, length(tsi))
        trf <- rep(NA,length(tsi))
        for(ti in 1:(length(bpf)-1)){
            trf[(bpf[ti]+1):bpf[ti+1]] <- cf[ti,1] + ((cf[ti,2]*((bpf[ti]+1):bpf[ti+1])))
        }
        # Find the major break
        dbr <- trf[tbp+1]-trf[tbp]
        tbp <- tbp[which(abs(dbr) == max(abs(dbr)))]
        # Calculate Frazier recovery metrics on BFAST trend component
        frz <- calcFrazier(as.numeric(trf), (tbp+1), floor(obspyr), shortDenseTS, nPre, nDist, nPostMin, nPostMax)
        # Calculate the post-disturbance slope of the BFAST trend component (first segment after break)
        sl <- (trf[tbp+3] - trf[tbp+2])
        frz <- c(frz, sl)
        names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
    }
    frz
}

#'
