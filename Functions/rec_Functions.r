library(strucchange)
library(bfast)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

toAnnualTS <- function(tsseas, tsi, obspyr){
    seasi <- rowMeans(matrix(tsseas, nrow = obspyr),na.rm=T)# average seasonality
    smax <- which(seasi == max(seasi, na.rm=T))# yearly observation number with max seas
    tsmi <- matrix(tsi, nrow = obspyr)
    dst <- abs(matrix(rep(1:obspyr,times = (length(tsi)/obspyr)), nrow = obspyr)-smax)# distance of observations to seasonal max
    dst[is.na(tsmi)] <- NA #set distance of NA observations equal to NA, these can not be selected
    rsel <- as.matrix(apply(dst, 2, which.min))# row numbers of observations to be selected, i.e. those closest to seasonal max
    toNA <- unlist(lapply(rsel, identical, integer(0)))# years without observation: assign temporary the first observation of the year
    rsel[toNA] <- 1
    rsel <- unlist(rsel) # rows to be selected
    csel <- 1:dim(tsmi)[2] # columns to be selected
    tsyr <- tsmi[rsel + nrow(tsmi) * (csel - 1)]# get values of yearly time series
    tsyr[toNA] <- NA# years without observation: set value to NA
    tsyr
}



calcFrazier <- function(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax){
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
            ARI <- max(tsio[(tdist +(4+obspyr)):(tdist+(5*obspyr))], na.rm=T) - V0
            # RRI: Relative Recovery Index (~recovery relative to impact)
            RRI <- ARI/Ddist
            if(is.infinite(RRI)){RRI <- NA}
            # R80p recovery index (~ ability to reach 80% of pre-disturbance value)
            R80P <- max(tsio[(tdist +(4+obspyr)):(tdist+(5*obspyr))], na.rm=T)/(Vpre*0.8)
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

