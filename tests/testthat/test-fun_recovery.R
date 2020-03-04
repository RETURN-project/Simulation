context("Recovery indicators")

test_that("Frazier - annual - too short time series", {

  tsio <- c(rep(0,1), seq(-1, 0), rep(0,1))
  tdist <- 2
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 1
  nDist <- 1
  nPostMin <- 1
  nPostMax <- 1

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
  nPre <- 1
  nDist <- 1
  nPostMin <- 1
  nPostMax <- 1

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dnbr <- 6
  ari <- 4

  rri <- ari/dnbr
  r80p <- -1/(0.8*pre)
  yryr <- (-2 + 5)/5

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})

test_that("Frazier - dense", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 12
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- mean(tsio[73:85])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (post - dist)/4.5


  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})

test_that("Frazier - segmented", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 12
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.1

  metrics <- calcBFASTrec(tsio, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- mean(tsio[73:85])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (post - dist)/4.5
  sl <- 4/60

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
  expect_equal(metrics$Sl, sl, tolerance = 1e-2)
})
