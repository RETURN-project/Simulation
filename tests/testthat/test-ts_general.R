context("General functions")

test_that("dense to annual", {
  tsseas <- rep(c(1,4,5,2),2)
  tsi <- c(2,5,7,4,5,6,7,3,4,5,6,7,2,1,3,4,5,3,6,9)
  obspyr <- 4

  tsa <- toAnnualTS(tsseas, tsi, obspyr)
  annual <- tsi[seq(3,20,by=4)]

  diff <- sum(tsa - annual)

  expect_equal(diff, 0, tolerance = 1e-4)
})

test_that("dense to annual with missing values", {
  tsseas <- rep(c(1,4,5,2),2)
  tsi <- c(2,5,NA,4,5,6,7,3,4,5,6,7,2,1,3,4,5,3,6,9)
  obspyr <- 4

  tsa <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/12)
  annual <- tsi[seq(3,20,by=4)]

  expect_equal(tsa, annual, tolerance = 1e-4)
})

test_that("dense to annual with varying intra-annual selection period", {
  tsseas <- rep(c(1,4,5,2),2)
  tsi <- c(2,5,NA,4,5,6,7,3,4,5,6,7,2,1,3,4,5,NA,NA,9)
  obspyr <- 4

  tsa1 <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/12)
  annual1 <- tsi[seq(3,20,by=4)]
  tsa2 <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/2)
  annual2 <- c(5, 7, 6, 3, 9)

  expect_equal(tsa1, annual1, tolerance = 1e-4)
  expect_equal(tsa2, annual2, tolerance = 1e-4)
})

test_that("ts compression", {
  tsi <- rep(NA,20)
  tsi[c(1,5,19)] <- c(7,3,9)

  tsc <- TScompress(tsi)
  tsd <- TSdecompress(tsc)
  diff <- sum(tsi - tsd, na.rm = T)

  expect_equal(diff, 0, tolerance = 1e-4)
})


test_that("Temporal aggregation", {
  # time series
  tsi <- c(1,2,1,20,1,30,-12,-2,-11,-21,-10,-30,-9,-39,-8,-48,-7,-57,-6,-66)

  # dates associated with observations
  dts <- as.Date(c('2001-01-02','2001-01-03','2001-02-02','2001-02-04','2001-03-02','2001-03-04','2001-04-02','2001-04-05','2001-05-02','2001-05-12',
                   '2001-06-02','2001-06-03','2001-07-02','2001-07-22','2001-08-02','2001-08-12','2001-09-02','2001-09-22','2001-12-02','2001-10-02'))

  brmomax <- toRegularTS(tsi,dts, fun = 'max', resol = 'monthly')
  brmomean <- toRegularTS(tsi,dts, fun = 'mean', resol = 'monthly')
  brday <- toRegularTS(tsi,dts, fun = 'max', resol = 'daily')
  brquart <- toRegularTS(tsi,dts, fun = 'max', resol = 'quart')

  # case 1 - monthly max
  expect_equal(as.numeric(brmomax), c(2,20,30,-2,-11,-10,-9,-8,-7,-66,NA,-6), tolerance = 1e-4)
  # case 2 - monthly mean
  expect_equal(as.numeric(brmomean), c(1.5,10.5,15.5,-7.0,-16.0,-20.0,-24.0,-28.0,-32.0,-66.0,NA,-6.0), tolerance = 1e-4)
  # case 3 - daily
  expect_equal(as.numeric(brday)[c(1,2,32,34,60,62,91,94,121,131,152,153,182,202,213,223,244,264,274,335)],
               c(1,2,1,20,1,30,-12,-2,-11,-21,-10,-30,-9,-39,-8,-48,-7,-57,-66,-6), tolerance = 1e-4)
  # case 4 - quarterly
  expect_equal(as.numeric(brquart),c(30,-2, -7,-6))
})

