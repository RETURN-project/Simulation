context("Performance evaluation")

test_that("perfect fit", {

  meas <- c(3,4,2,7,5,9)
  val <-meas

  MAPE <- mape(val, meas)
  R2 <- rsq(val, meas)
  RMSE <- rmse(val, meas)

  expect_equal(MAPE, 0, tolerance = 1e-4)
  expect_equal(RMSE, 0, tolerance = 1e-4)
  expect_equal(R2, 1, tolerance = 1e-4)
})

test_that("systematic error of 1", {

  meas <- c(3,4,2,7,5,9)
  val <-meas + 1

  MAPE <- mape(val, meas)
  R2 <- rsq(val, meas)
  RMSE <- rmse(val, meas)

  expect_equal(MAPE, mean(1/val), tolerance = 1e-4)
  expect_equal(RMSE, 1, tolerance = 1e-4)
  expect_equal(R2, 1, tolerance = 1e-4)
})

test_that("random error",{
  set.seed(197)
  meas <- runif(10000, min=0, max=100)
  set.seed(198)
  val <- runif(10000, min=0, max=100)

  R2 <- rsq(val, meas)
  expect_equal(R2, 0, tolerance = 1e-3)
})
