context("Simulation of benchmarking time series")

test_that("Exponential decay", {

  ts <- seq(0, 25, by = 0.1)
  pert <- 1.5
  tpert <- 0
  thalf <- 1

  ys <- exponential(ts, pert = pert, tpert = tpert, thalf = thalf)

  expect_equal(abs(min(ys)), 0, tolerance = 1e-4)
  expect_equal(abs(max(ys)), abs(pert), tolerance = 1e-4)

})

test_that("Linear decay", {

  ts <- seq(0, 5, by = 0.01)
  pert <- 1.5
  ys <- piecewise(ts, pert = pert)

  expect_equal(abs(min(ys)), 0, tolerance = 1e-4)
  expect_equal(abs(max(ys)), abs(pert), tolerance = 1e-4)
})
