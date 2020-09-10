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

test_that("plot sensitivity - basic test",{
  # fls <- list.files(path = "/home/wanda/Documents/Rscripts/BenchmarkRecovery/data/testdata/", pattern = "toyset_", full.names = TRUE)
  # for (i in 1:length(fls)){
  #   load(fls[i])}
  dat<- as.data.frame(list(
    Metric = as.factor(c('RRI', 'RRI',  'R80p','R80p','YrYr','YrYr','RRI','RRI','R80p','R80p','YrYr','YrYr','RRI','RRI','R80p','R80p','YrYr','YrYr','RRI','RRI','R80p', 'R80p','YrYr','YrYr')),
    Dense = as.factor(rep('anual',24)),
    Smooth = as.factor(rep('raw',24)),
    Period = as.factor(rep(c('Short','Long'),12)),
    Seas = as.factor(rep(F,24)),
    variable = as.factor(rep(c(0,0.08,0.02,0.03), each = 6)),
    value = c(0.0652595473, 0.0228917011, 0.1122423965, 0.1037905358, 0.0438699878, 0.0071708630, 0.0004348493, 0.0011731170, 0.1652815743, 0.1187443186,
              0.1445397388, 0.0341160032, 0.0447431895, 0.0063563672, 0.2544452752, 0.1292675690, 0.1124770530, 0.0019788734, 0.1212960214, 0.0437967248,
              0.1994489563, 0.2174260878, 0.0622801004, 0.0041313903)
  ))
  lbls <- c("raw, BAP")
  xlbl <- "Disturbance magnitude"
  ylbl <- "R²"
  scales <- 'free_y'
  pl <- plotSens(dat, lbls, xlbl, ylbl, scales)
  expect_identical(pl$labels$x, xlbl)
  expect_identical(pl$labels$y, ylbl)
  expect_identical(pl$labels$colour, "interaction(Smooth, Dense, Seas)")
  expect_identical(pl$labels$group, "interaction(Smooth, Dense, Seas)")
})
