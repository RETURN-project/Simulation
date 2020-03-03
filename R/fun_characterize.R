#' Decompose time series into trend, seasonality and remainder: This function decomposes time series into three components using BFAST01 functionality: trend, seasonality and remainder. Trends are fitted using linear regression without breaks, seasonality is fitted using a first order harmonic function and the remainder equals the anomalies (i.e. time series - trend - seasonality).
#'
#' @param df a dataframe with time series that need to be decomposed. The dataframe needs to be structured as follows: each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the time series values for each observation date.
#' @param nyr number of years of the input time series
#' @param nobsYr number of observations per year of the input time series
#'
#' @return a list containing the estimated seasonality, remainder, trend and seasonality coefficients. The seasonality is a dataframe with the seasonality of each pixel. Each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the seasonality values for each observation date. The trend and remainder are dataframes with the trend and remainder of each pixel (dataframe is structured in the same way as the seasonality). Seasonality_coefficients is a dataframe with the coeficients of the fitted harmonic function. Each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the coefficients of the fitted harmonic function.
#' @export
#' @import bfast
decompTSbfast <- function(df, nyr, nobsYr){
  #library(bfast)
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


#'  Fit ARMA model: This function automatically fits an ARMA model without seasonal component.
#'
#' @param tsx a time series object for which an ARMA model needs to be fitted.
#'
#' @return ARMA model coefficients. A list containing the ARMA coefficients and their order (number of AR coefficients, non-seasonal differences and MA coefficients
#' @export
#' @import forecast
getARMAcoef <- function(tsx){
  #library(forecast)
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

