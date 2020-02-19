#' ---
#' title: "General functions"
#' author: "Wanda De Keersmaecker"
#' date: "2/6/2020"
#' output: html_document
#' ---
#'
#'
#' # 1. Load RData file and returns it
#' This function is a substitute for the load function, allowing to assign a user defined variable name when loading a RData file.
#'
#' Function inputs:
#'
#' * fileName: the path to the Rdata file that needs to be loaded</font>
#'
## ------------------------------------------------------------------------
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' # 2. Time series compression
#' The TScompress function compresses a time series vector. The purpose is to avoid storing many NA values._
#'
#' Function inputs:
#'
#' * ts: a vector that will be compressed
#'
#' Function outputs:
#'
#' * a vector containing, in order, the length of the input vector, the number of observations without NA, the observations without NA, and the values of the observations that are no NA
#'
## ------------------------------------------------------------------------
TScompress <- function(ts){
    c(length(ts), length(which(!is.na(ts))), which(!is.na(ts)), ts[!is.na(ts)])
}

#' # 3. Time series decompression
#' The TSdecompress recovers a vector that has been compressed by the TScompress function in its original format.
#'
#' Function inputs:
#'
#' * ts: a vector compressed be the TScompress function
#'
#' Function outputs:
#'
#' * the vector restored in its original format
#'
## ------------------------------------------------------------------------
TSdecompress <- function(ts){
    vec <- rep(NA,ts[1])
    if(ts[2] > 0){
        vec[ts[3:(ts[2]+2)]] <- ts[(ts[2]+3):((2*ts[2])+2)]
    }
    vec
}

#' # 3. Convert time series to annual frequency
#' The toAnnualTS function converts a time series with n observations per year to an annual time series (one observation per year). The main concept is to select observations per year closest to a given day of year that have no missing value (NA). Here, the day of year for which the seasonality is maximum is being used.*
#'
#' Function inputs:
#'  - tsseas: vector of observations (time series) representing the seasonal component of the time series to be converted
#'  - tsi: vector of observations (time series) that needs to be converted to an annual time series
#'  - obspyr: number of observations per year of the time series to be converted
#'
#'
#'  Function outputs:
#'  - tsyr: vector of observations (time series) with annual observation frequency
#'
#'
## ------------------------------------------------------------------------
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


