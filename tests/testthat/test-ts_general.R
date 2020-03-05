context("General functions")

test_that("dense to annual", {
  tsseas <- c(1,4,5,2)
  tsi <- c(2,5,7,4,5,6,7,3,4,5,6,7,2,1,3,4,5,3,6,9)
  obspyr <- 4

  tsa <- toAnnualTS(tsseas, tsi, obspyr)
  annual <- tsi[seq(3,20,by=4)]

  diff <- sum(tsa - annual)

  expect_equal(diff, 0, tolerance = 1e-4)
})

test_that("ts compression", {
  tsi <- rep(NA,20)
  tsi[c(1,5,19)] <- c(7,3,9)

  tsc <- TScompress(tsi)
  tsd <- TSdecompress(tsc)
  diff <- sum(tsi - tsd, na.rm = T)

  expect_equal(diff, 0, tolerance = 1e-4)
})


