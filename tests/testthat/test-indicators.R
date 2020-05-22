context("Recovery indicator: YrYr")

test_that("Straight line time series", {
  # Generate a time series (straight line of slope m and intercept c)
  ts <- seq(-2, 10, by = 0.1)
  c <- 10
  m <- 3.14
  ys <- c + m * ts

  expect_equal(yryr(ts, ys), m, tolerance = 1e-6)
  expect_equal(yryr(ts, ys, tpert = -1, deltat = 3), m, tolerance = 1e-6)
})

test_that("Sinusoidal time series", {
  # Generate a time series (straight line of slope m and intercept c)
  ts <- seq(-2*pi, 10*pi, by = 0.01*pi)
  ys <- sin(ts)

  # Sampling after whole periods should yield a 0 mean slope
  expect_equal(yryr(ts, ys, tpert = 0, deltat = 1*pi), 0, tolerance = 1e-6)
  expect_equal(yryr(ts, ys, tpert = 0, deltat = 2*pi), 0, tolerance = 1e-6)
  expect_equal(yryr(ts, ys, tpert = pi, deltat = 2*pi), 0, tolerance = 1e-6)
})

test_that("NAs are tolerated", {
  # Generate a time series (straight line of slope m and intercept c)
  ts <- seq(-2*pi, 10*pi, by = 0.01*pi)
  ts[30] <- NA
  ys <- sin(ts)

  # Sampling after whole periods should yield a 0 mean slope
  expect_equal(yryr(ts, ys, tpert = 0, deltat = 1*pi), 0, tolerance = 1e-6)
  expect_equal(yryr(ts, ys, tpert = 0, deltat = 2*pi), 0, tolerance = 1e-6)
  expect_equal(yryr(ts, ys, tpert = pi, deltat = 2*pi), 0, tolerance = 1e-6)
})

test_that("Non matching lengths", {
  # Generate a bad time series (different vector sizes)
  ts <- seq(-2, 10, by = 0.1)
  ys <- seq(-2, 10, by = 0.2)

  expect_error(yryr(ts, ys),
               "'x' and 'y' lengths differ")
})

test_that("Out of bounds", {
  # Generate a bad time series (different vector sizes)
  ts <- seq(-2, 10, by = 0.1)
  ys <- seq(-2, 10, by = 0.1)

  expect_error(yryr(ts, ys, tpert = 9, deltat = 2),
               "*bounds*")
})

context("Recovery indicator: R80p")

test_that("Fully recovered time series", {
  # Generate a time series
  ts <- seq(-2, 10, by = 0.1)
  ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.5)

  expect_equal(r80p(ts, ys), 1 / 0.8, tolerance = 1e-2)
})

test_that("NAs are tolerated", {
  # Generate a time series
  ts <- seq(-2, 10, by = 0.1)
  ts[5] <- NA # Introduce NA
  ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.5)

  expect_equal(r80p(ts, ys), 1 / 0.8, tolerance = 1e-2)
})

context("Recovery indicator: RRI")

test_that("Fully recovered time series", {
  # Generate a time series
  ts <- seq(-2, 10, by = 0.1) # as a vector of times
  ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.25) # plus a vector of values

  expect_equal(rri(ts, ys), 1, tolerance = 1e-4)
})

test_that("NAs are tolerated", {
  # Generate a time series
  ts <- seq(-2, 10, by = 0.1)
  ts[5] <- NA # Introduce NA
  ys <- exponential(ts, pert = -2, offset = 1, thalf = 0.25) # plus a vector of values

  expect_equal(rri(ts, ys), 1, tolerance = 1e-2)
})
