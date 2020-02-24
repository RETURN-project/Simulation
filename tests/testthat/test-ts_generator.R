context("Simulation of benchmarking time series")

test_that("Exponential decay", {

  ts <- seq(0, 10, by=0.1)
  pert <- 1.5
  ys <- exponential(ts, pert=pert)

  expect_equal(min(ys), 0, tolerance = 1e-4)
  expect_equal(max(ys), pert, tolerance = 1e-4)
})
