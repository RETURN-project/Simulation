context("Recovery indicators")

test_that("Frazier - annual - too short time series", {

  tsio <- c(rep(0,1), seq(-1, 0), rep(0,1))
  tdist <- 2
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)

  expect_equal(metrics$RRI, NA)
  expect_equal(metrics$R80P, NA)
  expect_equal(metrics$YrYr, NA)
})

test_that("Frazier - annual", {

  tsio <- c(rep(1,2), seq(-5, -1), rep(-2,1))
  tdist <- 3
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dnbr <- 6
  ari <- 4

  rrim <- ari/dnbr
  r80pm <- -1/(0.8*pre)
  yryrm <- (-2 + 5)/5

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - dense", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 1
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- max(tsio[73:84])
  dnbr <- pre-dist
  ari <- post-dist

  rrim <- ari/dnbr
  r80pm <- post/(0.8*pre)
  yryrm <- (mean(tsio[73:84]) - dist)/(4*12)

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - segmented", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 1
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.1

  metrics <- calcBFASTrec(tsio, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- max(tsio[73:84])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (mean(tsio[73:84]) - dist)/(4*12)

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})


test_that("Frazier - segmented annual - long", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0,8))
  tdist <- 9
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.2
  seas <- F

  metrics <- calcBFASTrec(tsio, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, seas = F)
  pre <- 1
  dnbr <- 6
  ari <- 2.5

  rrim <- ari/dnbr
  r80pm <- tsio[14]/(0.8*pre)
  yryrm <- (tsio[14] - tsio[9])/5

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - segmented annual - short", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0,8))
  tdist <- 9
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 1
  nPostMax <- 1
  h <- 0.2
  seas <- F

  metrics <- calcBFASTrec(tsio, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, seas = F)
  pre <- 1
  dnbr <- 6
  ari <- 0.5

  rrim <- ari/dnbr
  r80pm <- tsio[10]/(0.8*pre)
  yryrm <- (tsio[10] - tsio[9])/1

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})
