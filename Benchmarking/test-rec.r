library(testthat)

test_that("Test RRI dense raw",{
nfolder <- 'C:\\Users\\keers001\\Dropbox\\output\\Jupyter_notebook\\'
system(paste0('jupyter-nbconvert.exe ', nfolder, 'rec_Functions.ipynb  --to script'))#, intern=FALSE
source(file.path(nfolder, 'rec_Functions.r'))
tsi <- c(rep(5,24), rep(1,12),rep(3,24), rep(4,26))
tdist <- 25
obspyr <- 12

shortDenseTS <- T
nPre <- 2
nDist <- 12
nPostMin <- 3
nPostMax <-5

rec <- calcFrazier(tsi, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
    expect_equal(rec[['RRI']], (4-1)/(5-1))
    expect_equal(rec[['R80P']], (4)/(5*0.8))
    expect_equal(rec[['YrYr']], (4-1)/4)
})

test_that("Test RRI annual raw",{
nfolder <- 'C:\\Users\\keers001\\Dropbox\\output\\Jupyter_notebook\\'
system(paste0('jupyter-nbconvert.exe ', nfolder, 'rec_Functions.ipynb  --to script'))#, intern=FALSE
source(file.path(nfolder, 'rec_Functions.r'))
tsi <- c(rep(5,24), rep(1,12),rep(3,24), rep(4,36))
tsi <- toAnnualTS(rep(c(1,2,3,4,5,6,7,6,5,4,3,2),8),tsi,12)
tdist <- ceiling(25/12)
obspyr <- 1

shortDenseTS <- F
nPre <- 2
nDist <- 12
nPostMin <- 3
nPostMax <-5

#(4-1)/(5-1)#RRI

#(4)/(5*0.8)#R80p

#(4-1)/5#YrYr
rec <- calcFrazier(tsi, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
    expect_equal(rec[['RRI']], (4-1)/(5-1))
    expect_equal(rec[['R80P']], (4)/(5*0.8))
    expect_equal(rec[['YrYr']], (4-1)/5)
})





