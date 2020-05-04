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
#' @param seas TRUE or FALSE, include seasonal term when detecting breaks?
#' @param breaks 'BIC' or 'LWZ': criteria used to define the optimal number of breaks
#'
#' @return a list containing  the RRI, R80p, YrYr recovery indicator derived from the BFAST0n trend segments and slope of the trend segment after the disturbance (sl).
#' @export
#' @import strucchange
#' @import stats
calcBFASTrec <- function(tsio, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, breaks = 'BIC', seas = F){
  # Create time series object, needed as input for BFAST
  tsi <- ts(tsio, frequency = obspyr)
  # Convert the time series object into a dataframe, needed for the breakpoints function
    datapp <- bfastpp(tsi, order = 1, lag = NULL, slag = NULL,
                  na.action = na.omit, stl = 'none')
    nreg <- switch(seas, 5, 2)
    # Test if enough observations are available to use time series segmentation
    if(round(length(tsio[is.na(tsio)==F]) * h) > nreg){
      # set_fast_options()
      # Apply BFAST0n on time series: find breaks in the regression
      if (seas){
        bp <- breakpoints(response ~ trend + harmon, data = datapp, h = h, breaks = breaks)#, breaks = breaks
      } else{
        bp <- breakpoints(response ~ trend, data = datapp, h = h, breaks = breaks)##, breaks = breaks
      }
      # Check if BFAST0n found breakpoints
      if(is.na(bp$breakpoints[1])){# no breakpoint found
        frz <- list(NA, NA, NA, NA)
        names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
      }else{# at least one breakpoint found
        # Extract trend component and breaks
        cf <- coef(bp, breaks = breaks)#, breaks = breaks
        # Extract trend component and breaks
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
        # Calculate Frazier recovery metrics on trend component
        frz <- calcFrazier(as.numeric(trf), (tbp+1), floor(obspyr), shortDenseTS, nPre, nDist, nPostMin, nPostMax)
        # Calculate the post-disturbance slope of the trend component (first segment after break)
        sl <- (trf[tbp+3] - trf[tbp+2])
        frz <- c(frz, sl)
        names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
      }
    }else{
      frz <- list(NA, NA, NA, NA)
      names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
    }

    frz
}

#
#' convert to matrix of performance indicators to a dataframe
#'
#' @param mat matrix of recovery indicators (rows represent the values of the tested time series characteristic, colums represent the different recovery indicators)
#' @param setvr time series characteristic being evaluated
#' @param metric recovery metric being evaluated
#' @param freq vector of temporal frequencies per recovery indicator
#' @param input vector of preprocessing techniques per recovery indicator
#' @param nDist vector of the time span during the disturbance used to measure recovery
#'
#' @return dataframe
#' @export
toDF <- function(mat, setvr, metric, freq, input, nDist, breaks, seas){
  tst <- as.data.frame(t(mat))
  names(tst) <- setvr
  tst$Metric <- factor(metric)
  tst$Dense <- factor(freq)
  tst$Smooth <- factor(input)
  tst$Period <- revalue(factor(nDist), c("1"="Short", "12"="Long"))#factor(recSttngs$nDist)#
  tst$Breaks <- factor(breaks)
  tst$Seas <- factor(seas)
  tst
}

#' Run sensitivity analysis for particular parameter of interest
#'
#' @param vr integer: defines the target parameter for the sensitivity study. The vr'th 'eval' parameter in the settings list is selected as target
#' @param sttngs list of settings
#' @param pars list of simulation parameters
#' @param funSet list of recovery indicator settings
#' @param ofolder folder where the output files should be stored
#' @param basename basename for the output files
#'
#' @return saves performance indicators to the output folder (R2, RMSE, MAPE) and indicators of the relation between the measured and simulated recovery indicators (slope, intercept, p value of normality test of residuals). For each performace indicator a matrix is saved where the rows refer to the evaluated values of the paramter of interest and the columns to the various recovery indicator settings.
#' @export
#'
evalParam <- function(vr, sttngs, pars, funSet, ofolder, basename){
  evr <- sttngs$general$eval[vr]# name of parameter that will be evaluated in the simulation
  parvr <- pars[[evr]]

  RRI_rmse <- matrix(NA,length(parvr), length(funSet[[1]]))
  RRI_mape <- matrix(NA,length(parvr), length(funSet[[1]]))
  RRI_rsq <- matrix(NA,length(parvr), length(funSet[[1]]))
  RRI_nTS <- matrix(NA,length(parvr), length(funSet[[1]]))
  RRI_int <- matrix(NA,length(parvr), length(funSet[[1]]))
  RRI_slope <- matrix(NA,length(parvr), length(funSet[[1]]))
  RRI_norm <- matrix(NA,length(parvr), length(funSet[[1]]))

  R80p_rmse <- matrix(NA,length(parvr), length(funSet[[1]]))
  R80p_mape <- matrix(NA,length(parvr), length(funSet[[1]]))
  R80p_rsq <- matrix(NA,length(parvr), length(funSet[[1]]))
  R80p_nTS <- matrix(NA,length(parvr), length(funSet[[1]]))
  R80p_int <- matrix(NA,length(parvr), length(funSet[[1]]))
  R80p_slope <- matrix(NA,length(parvr), length(funSet[[1]]))
  R80p_norm <- matrix(NA,length(parvr), length(funSet[[1]]))

  YrYr_rmse <- matrix(NA,length(parvr), length(funSet[[1]]))
  YrYr_mape <- matrix(NA,length(parvr), length(funSet[[1]]))
  YrYr_rsq <- matrix(NA,length(parvr), length(funSet[[1]]))
  YrYr_nTS <- matrix(NA,length(parvr), length(funSet[[1]]))
  YrYr_int <- matrix(NA,length(parvr), length(funSet[[1]]))
  YrYr_slope <- matrix(NA,length(parvr), length(funSet[[1]]))
  YrYr_norm <- matrix(NA,length(parvr), length(funSet[[1]]))

  SL_rmse <- matrix(NA,length(parvr), length(funSet[[1]]))
  SL_mape <- matrix(NA,length(parvr), length(funSet[[1]]))
  SL_rsq <- matrix(NA,length(parvr), length(funSet[[1]]))
  SL_nTS <- matrix(NA,length(parvr), length(funSet[[1]]))
  SL_int <- matrix(NA,length(parvr), length(funSet[[1]]))
  SL_slope <- matrix(NA,length(parvr), length(funSet[[1]]))
  SL_norm <- matrix(NA,length(parvr), length(funSet[[1]]))


  # iterate over values of evaluated parameter and simulate nrep time series per combination of all other variables
  for (i in 1:length(parvr)){#length(setvr)

    m_RRIi <- matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))#
    m_R80pi <-  matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
    m_YrYri <-  matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
    m_SLi <-  matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
    s_RRIi <- matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
    s_R80pi <-  matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
    s_YrYri <-  matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
    s_SLi <-  matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))

    # iterate over the parameter settings and simulate each time a time series
    for (pari in 1: length(parvr[[1]][[1]])){
      # simulate time series for a parameter combination
      sc <- simulCase(parvr[[i]]$nrep[pari], parvr[[i]]$nyr[pari], parvr[[i]]$nobsYr[pari], parvr[[i]]$nDr[pari], parvr[[i]]$seasAv[[pari]], parvr[[i]]$seasAmp[pari],parvr[[i]]$trAv[pari], parvr[[i]]$remSd[pari], c(parvr[[i]]$distMag[pari],parvr[[i]]$distMag[pari]), parvr[[i]]$distT[pari], c(parvr[[i]]$distRec[pari],parvr[[i]]$distRec[pari]), parvr[[i]]$remcoef[pari], parvr[[i]]$missVal[pari], parvr[[i]]$DistMissVal[pari], parvr[[i]]$distType[pari])

      # iterate over the recovery indicator settings
      for (rset in 1:length(funSet[[1]])){#1:length(funSet[[1]])
        tsi <- sc[[1]][1,]
        tsseas <-sc[[3]][1,]
        obspyr <- sc[[6]][1,]$obs_per_year
        tdist <- sc[[6]][1,]$dist_time
        tsref <- sc[[2]][1,] + sc[[5]][1,]
        nobs <- (sc[[6]][1,]$number_yrs)* obspyr
        tm <- 1:nobs

        inp <- funSet$input[rset]
        frq <- funSet[['freq']][rset]
        shortDenseTS <- funSet[['shortDenseTS']][rset]
        nPre <- funSet[['nPre']][rset]
        nDist <- funSet[['nDist']][rset]
        nPostMin <- funSet[['nPostMin']][rset]
        nPostMax <- funSet[['nPostMax']][rset]
        h <- funSet[['h']][rset]
        seas <- funSet[['seas']][rset]
        breaks <- funSet[['breaks']][rset]

        if (frq == 'annual'){
          #convert time series to annual values by selecting date closest to seasonal max
          tsi <- toAnnualTS(tsseas, tsi, obspyr)
          tsref <- toAnnualTS(tsseas, tsref, obspyr)
          tdist <- ceiling(tdist/obspyr)
          obspyr <- 1
        }
        if (frq == 'quarterly'){
          #convert time series to quarterly resolution
          dts <- seq(as.Date('2000-01-01'), by = '1 days', length = nobs)
          tsi <- toRegularTS(tsi, dts, 'mean', 'quart')
          tsref <- toRegularTS(tsref, dts, 'mean', 'quart')
          tdist <- ceiling(tdist/obspyr)
          obspyr <- 4
        }

        if (inp == 'smoothed'){
          temp.zoo<-zoo(tsi,(1:length(tsi)))
          m.av<-rollapply(temp.zoo, 150, mean, na.rm = T, fill = NA)
          tsi <- as.numeric(m.av)

        }

        if((inp == 'smoothed') | (inp == 'raw')){
          outp <- calcFrazier(tsi, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
          m_RRIi[rset,pari] <- outp$RRI# measured RRI
          m_R80pi[rset,pari] <- outp$R80P# measured R80p
          m_YrYri[rset,pari] <- outp$YrYr# measured YrYR
          rm(outp)
        }
        if(funSet$input[rset] == 'segmented'){
          outp <- calcBFASTrec(tsi, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, breaks, seas)
          m_RRIi[rset,pari] <- outp$RRI# measured RRI
          m_R80pi[rset,pari] <- outp$R80P# measured R80p
          m_YrYri[rset,pari] <- outp$YrYr# measured YrYR
          m_SLi[rset,pari] <- outp$Sl# measured YrYR
          #print(outp)
          rm(outp)
        }

        # reference indicators
        outp <- calcFrazier(tsref, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
        s_RRIi[rset,pari] <- outp$RRI#simulated (true) RRI
        s_R80pi[rset,pari] <- outp$R80P#simulated (true) R80p
        s_YrYri[rset,pari] <- outp$YrYr
        if(inp == 'segmented'){
          s_SLi[rset,pari] <- tsref[tdist+3] - tsref[tdist+2]
        }
      }
      rm(sc)
    }
    # evaluate recovery indicators
    # R2
    RRI_rsq[i,] <- sapply(1:dim(s_RRIi)[1], function(it) rsq(s_RRIi[it,], m_RRIi[it,]))
    R80p_rsq[i,] <- sapply(1:dim(s_R80pi)[1], function(it) rsq(s_R80pi[it,], m_R80pi[it,]))
    YrYr_rsq[i,] <- sapply(1:dim(s_YrYri)[1], function(it) rsq(s_YrYri[it,], m_YrYri[it,]))
    SL_rsq[i,] <- sapply(1:dim(s_SLi)[1], function(it) rsq(s_SLi[it,], m_SLi[it,]))

    # MAPE
    RRI_mape[i,] <- sapply(1:dim(s_RRIi)[1], function(it) mape(s_RRIi[it,], m_RRIi[it,]))
    R80p_mape[i,] <- sapply(1:dim(s_R80pi)[1], function(it) mape(s_R80pi[it,], m_R80pi[it,]))
    YrYr_mape[i,] <- sapply(1:dim(s_YrYri)[1], function(it) mape(s_YrYri[it,], m_YrYri[it,]))
    SL_mape[i,] <- sapply(1:dim(s_SLi)[1], function(it) mape(s_SLi[it,], m_SLi[it,]))

    # RMSE
    RRI_rmse[i,] <- sapply(1:dim(s_RRIi)[1], function(it) rmse(s_RRIi[it,], m_RRIi[it,]))
    R80p_rmse[i,] <- sapply(1:dim(s_R80pi)[1], function(it) rmse(s_R80pi[it,], m_R80pi[it,]))
    YrYr_rmse[i,] <- sapply(1:dim(s_YrYri)[1], function(it) rmse(s_YrYri[it,], m_YrYri[it,]))
    SL_rmse[i,] <- sapply(1:dim(s_SLi)[1], function(it) rmse(s_SLi[it,], m_SLi[it,]))

    # nTS
    RRI_nTS[i,] <- apply(m_RRIi, 1, function(x){sum(is.na(x)==F)/length(x)})
    R80p_nTS[i,] <- apply(m_R80pi, 1, function(x){sum(is.na(x)==F)/length(x)})
    YrYr_nTS[i,] <- apply(m_YrYri, 1, function(x){sum(is.na(x)==F)/length(x)})
    SL_nTS[i,] <- apply(m_SLi, 1, function(x){sum(is.na(x)==F)/length(x)})

    # intercept
    RRI_int[i,] <- sapply(1:dim(s_RRIi)[1], function(it) linFit(s_RRIi[it,], m_RRIi[it,]))[1,]
    R80p_int[i,] <- sapply(1:dim(s_R80pi)[1], function(it) linFit(s_R80pi[it,], m_R80pi[it,]))[1,]
    YrYr_int[i,] <- sapply(1:dim(s_YrYri)[1], function(it) linFit(s_YrYri[it,], m_YrYri[it,]))[1,]
    SL_int[i,] <- sapply(1:dim(s_SLi)[1], function(it) linFit(s_SLi[it,], m_SLi[it,]))[1,]

    # slope
    RRI_slope[i,] <- sapply(1:dim(s_RRIi)[1], function(it) linFit(s_RRIi[it,], m_RRIi[it,]))[2,]
    R80p_slope[i,] <- sapply(1:dim(s_R80pi)[1], function(it) linFit(s_R80pi[it,], m_R80pi[it,]))[2,]
    YrYr_slope[i,] <- sapply(1:dim(s_YrYri)[1], function(it) linFit(s_YrYri[it,], m_YrYri[it,]))[2,]
    SL_slope[i,] <- sapply(1:dim(s_SLi)[1], function(it) linFit(s_SLi[it,], m_SLi[it,]))[2,]

    # Shapiro-Wilk normality test
    RRI_norm[i,] <- sapply(1:dim(s_RRIi)[1], function(it) linFit(s_RRIi[it,], m_RRIi[it,]))[3,]
    R80p_norm[i,] <- sapply(1:dim(s_R80pi)[1], function(it) linFit(s_R80pi[it,], m_R80pi[it,]))[3,]
    YrYr_norm[i,] <- sapply(1:dim(s_YrYri)[1], function(it) linFit(s_YrYri[it,], m_YrYri[it,]))[3,]
    SL_norm[i,] <- sapply(1:dim(s_SLi)[1], function(it) linFit(s_SLi[it,], m_SLi[it,]))[3,]
  }

  RRI_rsqDF <- toDF(RRI_rsq, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_rsqDF <- toDF(R80p_rsq, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_rsqDF <- toDF(YrYr_rsq, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_rsqDF <- toDF(SL_rsq, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  RRI_mapeDF <- toDF(RRI_mape, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_mapeDF <- toDF(R80p_mape, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_mapeDF <- toDF(YrYr_mape, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_mapeDF <- toDF(SL_mape, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  RRI_rmseDF <- toDF(RRI_rmse, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_rmseDF <- toDF(R80p_rmse, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_rmseDF <- toDF(YrYr_rmse, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_rmseDF <- toDF(SL_rmse, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  RRI_nTSDF <- toDF(RRI_nTS, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_nTSDF <- toDF(R80p_nTS, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_nTSDF <- toDF(YrYr_nTS, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_nTSDF <- toDF(SL_nTS, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  RRI_intDF <- toDF(RRI_int, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_intDF <- toDF(R80p_int, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_intDF <- toDF(YrYr_int, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_intDF <- toDF(SL_int, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  RRI_slopeDF <- toDF(RRI_slope, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_slopeDF <- toDF(R80p_slope, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_slopeDF <- toDF(YrYr_slope, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_slopeDF <- toDF(SL_slope, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  RRI_normDF <- toDF(RRI_norm, names(parvr), 'RRI', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  R80p_normDF <- toDF(R80p_norm, names(parvr), 'R80p', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  YrYr_normDF <- toDF(YrYr_norm, names(parvr), 'YrYr', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)
  SL_normDF <- toDF(SL_norm, names(parvr), 'SL', funSet$freq, funSet$input, funSet$nDist, funSet$breaks, funSet$seas)

  # export the performance indicators
  save(RRI_rsqDF, file = file.path(ofolder, paste0(basename, '_RRI_R2_' , evr, '.rda')))
  save(R80p_rsqDF, file = file.path(ofolder, paste0(basename, '_R80p_R2_' , evr, '.rda')))
  save(YrYr_rsqDF, file = file.path(ofolder, paste0(basename, '_YrYr_R2_' , evr, '.rda')))
  save(SL_rsqDF, file = file.path(ofolder, paste0(basename, '_SL_R2_' , evr, '.rda')))

  save(RRI_rmseDF, file = file.path(ofolder, paste0(basename, '_RRI_RMSE_' , evr, '.rda')))
  save(R80p_rmseDF, file = file.path(ofolder, paste0(basename, '_R80p_RMSE_' , evr, '.rda')))
  save(YrYr_rmseDF, file = file.path(ofolder, paste0(basename, '_YrYr_RMSE_' , evr, '.rda')))
  save(SL_rmseDF, file = file.path(ofolder, paste0(basename, '_SL_RMSE_' , evr, '.rda')))

  save(RRI_mapeDF, file = file.path(ofolder, paste0(basename, '_RRI_MAPE_' , evr, '.rda')))
  save(R80p_mapeDF, file = file.path(ofolder, paste0(basename, '_R80p_MAPE_' , evr, '.rda')))
  save(YrYr_mapeDF, file = file.path(ofolder, paste0(basename, '_YrYr_MAPE_' , evr, '.rda')))
  save(SL_mapeDF, file = file.path(ofolder, paste0(basename, '_SL_MAPE_' , evr, '.rda')))

  save(RRI_nTSDF, file = file.path(ofolder, paste0(basename, '_RRI_nTS_' , evr, '.rda')))
  save(R80p_nTSDF, file = file.path(ofolder, paste0(basename, '_R80p_nTS_' , evr, '.rda')))
  save(YrYr_nTSDF, file = file.path(ofolder, paste0(basename, '_YrYr_nTS_' , evr, '.rda')))
  save(SL_nTSDF, file = file.path(ofolder, paste0(basename, '_SL_nTS_' , evr, '.rda')))

  save(RRI_intDF, file = file.path(ofolder, paste0(basename, '_RRI_int_' , evr, '.rda')))
  save(R80p_intDF, file = file.path(ofolder, paste0(basename, '_R80p_int_' , evr, '.rda')))
  save(YrYr_intDF, file = file.path(ofolder, paste0(basename, '_YrYr_int_' , evr, '.rda')))
  save(SL_intDF, file = file.path(ofolder, paste0(basename, '_SL_int_' , evr, '.rda')))

  save(RRI_slopeDF, file = file.path(ofolder, paste0(basename, '_RRI_slope_' , evr, '.rda')))
  save(R80p_slopeDF, file = file.path(ofolder, paste0(basename, '_R80p_slope_' , evr, '.rda')))
  save(YrYr_slopeDF, file = file.path(ofolder, paste0(basename, '_YrYr_slope_' , evr, '.rda')))
  save(SL_slopeDF, file = file.path(ofolder, paste0(basename, '_SL_slope_' , evr, '.rda')))

  save(RRI_normDF, file = file.path(ofolder, paste0(basename, '_RRI_norm_' , evr, '.rda')))
  save(R80p_normDF, file = file.path(ofolder, paste0(basename, '_R80p_norm_' , evr, '.rda')))
  save(YrYr_normDF, file = file.path(ofolder, paste0(basename, '_YrYr_norm_' , evr, '.rda')))
  save(SL_normDF, file = file.path(ofolder, paste0(basename, '_SL_norm_' , evr, '.rda')))

}

