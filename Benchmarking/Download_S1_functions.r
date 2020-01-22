RunS1app <- function(mail, pword, lat, lon, buff, startdate, enddate, polVV, polVH, polDIFF, 
            orbitASC, orbitDES, prepSlope, prepSpeckle, SpeckleMulti, BSGamma,
            ScaleDB, DLperDate, maxWait, maxWaitShort, SCshot, SCpath){
  
  #-------------------
  # function to swith active window
    shortSleep <- 5
  myswitch <- function (remDr, windowId) 
  {
    qpath <- sprintf("%s/session/%s/window", remDr$serverURL, 
                     remDr$sessionInfo[["id"]])
    remDr$queryRD(qpath, "POST", qdata = list(handle = windowId))
  }
  
  #-------------------
  # login with google
    library('RSelenium')
    library('pingr')
    #porti <- 4444
    #portAv <- sum(is.na(pingr::ping_port("localhost", porti)))
    #while(portAv != 3){
    #    porti <- porti+1
    #    portAv <- sum(is.na(pingr::ping_port("localhost", porti)))
    #}
    rD <- rsDriver(verbose = FALSE, browser=c("chrome"), chromever="79.0.3945.36")#, extraCapabilities = eCaps)#port=as.integer(porti),
    remDr <- rD$client#78.0.3904.105,79.0.3945.16,79.0.3945.36,80.0.3987.16
    
    remDr$setTimeout(type = "implicit", milliseconds = maxWait*1000)
    remDr$setTimeout(type = "page load", milliseconds = maxWait*1000)
  
  
  remDr$navigate("https://wurnrt-s1ard.appspot.com/")
  remDr$getTitle()
  
  curWait <- 0
  while((! startsWith(remDr$getTitle()[[1]], "GEE Sentinel-1 Analysis Ready Data")) && (curWait < maxWaitShort)){
    Sys.sleep(1)
    curWait <- curWait + 1
  }
  Sys.sleep(shortSleep)
  
  webElems <- remDr$findElements(using = "class", value = "button")
  webElems[[1]]$clickElement()
  rm(webElems)
  
  #-------------------
  # enter emailadres
  currWin <- remDr$getCurrentWindowHandle()
  allWins <- unlist(remDr$getWindowHandles())
  curWait <- 0
  while((length(allWins)==1) && (curWait < maxWaitShort)){
    Sys.sleep(1)
    curWait <- curWait + 1
    allWins <- unlist(remDr$getWindowHandles())
  }
  otherWindow <- allWins[!allWins %in% currWin[[1]]]
  myswitch(remDr, otherWindow)
  
  curWait <- 0
  while((! startsWith(remDr$getTitle()[[1]], "Inloggen - Google Accounts")) && (curWait < maxWaitShort)){
    Sys.sleep(1)
    curWait <- curWait + 1
  }
  Sys.sleep(shortSleep)
  webElems <- remDr$findElements(using = "name", value = "identifier")
  # webElems[[1]]$highlightElement()
  webElems[[1]]$sendKeysToElement(list(mail))
  rm(webElems)
  webElems <- remDr$findElements(using = "id", value = "identifierNext")
  webElems[[1]]$clickElement()
  rm(webElems)
  
  #-------------------
  # enter password
  webElems <- remDr$findElements(using = "name", value = "password")
  curWait <- 0
  while((length(webElems)==0) && (curWait < maxWaitShort)){
    Sys.sleep(1)
    curWait <- curWait + 1
    webElems <- remDr$findElements(using = "name", value = "password")
  }
  Sys.sleep(shortSleep)
  webElems[[1]]$sendKeysToElement(list(pword))
  rm(webElems)
  webElems <- remDr$findElements(using = "id", value = "passwordNext")
  webElems[[1]]$clickElement()
  rm(webElems)
  
  #-------------------
  # enter data
  allWins <- unlist(remDr$getWindowHandles())
  while((length(allWins)>1) && (curWait < maxWaitShort)){
    Sys.sleep(1)
    curWait <- curWait + 1
    allWins <- unlist(remDr$getWindowHandles())
  }
  Sys.sleep(shortSleep)
  
  myswitch(remDr, allWins)
  # remDr$setWindowSize(1920, 1080)
  
  curWait <- 0
  while((! startsWith(remDr$getTitle()[[1]], "GEE Sentinel-1 Analysis Ready Data")) && (curWait < maxWaitShort)){
    Sys.sleep(1)
    curWait <- curWait + 1
  }
  Sys.sleep(shortSleep)
  
  #-------------------
  #lat
  webElems <- remDr$findElements(using = "name", value = "lat")
  webElems[[1]]$clearElement()
  webElems[[1]]$sendKeysToElement(list(lat))
  rm(webElems)
  
  #-------------------
  #lon
  webElems <- remDr$findElements(using = "name", value = "lng")
  webElems[[1]]$clearElement()
  webElems[[1]]$sendKeysToElement(list(lon))
  rm(webElems)
  
  #-------------------
  #buffer
  webElems <- remDr$findElements(using = "name", value = "buf")
  webElems[[1]]$clearElement()
  webElems[[1]]$sendKeysToElement(list(buff))
  rm(webElems)
  
  #-------------------
  #start date
  webElems <- remDr$findElements(using = "id", value = "startdate")
  webElems[[1]]$clearElement()
  webElems[[1]]$sendKeysToElement(list(startdate, "\uE007"))
  rm(webElems)
  
  #-------------------
  #end date
  webElems <- remDr$findElements(using = "id", value = "enddate")
  webElems[[1]]$clearElement()
  webElems[[1]]$sendKeysToElement(list(enddate, "\uE007"))
  rm(webElems)
  
  #-------------------
  #polarization
  #VV
  remDr$executeScript("window.scrollTo(0,document.body.scrollHeight);")
  webElems <- remDr$findElements(using = "name", value = "pol")
  
  if(! webElems[[1]]$isElementSelected()[[1]] && polVV){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! polVV){#VV
    webElems[[1]]$clickElement()
  }
  
  #VH
  if(! webElems[[2]]$isElementSelected()[[1]] && polVH){#VV
    webElems[[2]]$clickElement()
  }
  if(webElems[[2]]$isElementSelected()[[1]] && ! polVH){#VV
    webElems[[2]]$clickElement()
  }
  
  #VV-VH
  if(! webElems[[3]]$isElementSelected()[[1]] && polDIFF){#VV
    webElems[[3]]$clickElement()
  }
  if(webElems[[3]]$isElementSelected()[[1]] && ! polDIFF){#VV
    webElems[[3]]$clickElement()
  }
  rm(webElems)
  # if(SCshot){remDr$screenshot(file = paste0(SCpath,'_A.png'))}
  #-------------------
  #orbit
  
  webElems <- remDr$findElements(using = "name", value = "orbit")
  #Ascending
  if(! webElems[[1]]$isElementSelected()[[1]] && orbitASC){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! orbitASC){#VV
    webElems[[1]]$clickElement()
  }
  
  #Descending
  if(! webElems[[2]]$isElementSelected()[[1]] && orbitDES){#VV
    webElems[[2]]$clickElement()
  }
  if(webElems[[2]]$isElementSelected()[[1]] && ! orbitDES){#VV
    webElems[[2]]$clickElement()
  }
  rm(webElems)
  
  webElem <- remDr$findElement("css", "body")
  webElem$sendKeysToElement(list(key = "end"))
  webElems <- remDr$findElements(using = "name", value = "orbit")

  #Ascending
  if(! webElems[[1]]$isElementSelected()[[1]] && orbitASC){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! orbitASC){#VV
    webElems[[1]]$clickElement()
  }

  #Descending
  if(! webElems[[2]]$isElementSelected()[[1]] && orbitDES){#VV
    webElems[[2]]$clickElement()
  }
  if(webElems[[2]]$isElementSelected()[[1]] && ! orbitDES){#VV
    webElems[[2]]$clickElement()
  }
  rm(webElems)
  
  #-------------------
  #preprocessing
  # Slope correction
  webElems <- remDr$findElements(using = "name", value = "preproc")
  if(! webElems[[1]]$isElementSelected()[[1]] && prepSlope){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! prepSlope){#VV
    webElems[[1]]$clickElement()
  }
  # speckle filter
  if(! webElems[[2]]$isElementSelected()[[1]] && prepSpeckle){#VV
    webElems[[2]]$clickElement()
  }
  if(webElems[[2]]$isElementSelected()[[1]] && ! prepSpeckle){#VV
    webElems[[2]]$clickElement()
  }
  rm(webElems)
  
  # speckle filter - Refined Lee vs Multi-temporal Lee
  if(prepSpeckle){
    webElems <- remDr$findElements(using = "name", value = "leefilter")
    if(! webElems[[2]]$isElementSelected()[[1]] && SpeckleMulti){#VV
      webElems[[2]]$clickElement()
    }
    if(webElems[[2]]$isElementSelected()[[1]] && ! SpeckleMulti){#VV
      webElems[[1]]$clickElement()
    }
  }
  rm(webElems)
  
  
  #-------------------
  #backscatter
  webElems <- remDr$findElements(using = "name", value = "backscatter")
  if(! webElems[[1]]$isElementSelected()[[1]] && BSGamma){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! BSGamma){#VV
    webElems[[2]]$clickElement()
  }
  rm(webElems)
  
  #-------------------
  #scale
  webElems <- remDr$findElements(using = "name", value = "imgscale")
  if(! webElems[[1]]$isElementSelected()[[1]] && ScaleDB){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! ScaleDB){#VV
    webElems[[2]]$clickElement()
  }
  rm(webElems)
  
  #-------------------
  #download type
  webElems <- remDr$findElements(using = "name", value = "download")
  if(! webElems[[1]]$isElementSelected()[[1]] && DLperDate){#VV
    webElems[[1]]$clickElement()
  }
  if(webElems[[1]]$isElementSelected()[[1]] && ! DLperDate){#VV
    webElems[[2]]$clickElement()
  }
  rm(webElems)
  
 #    remDr$executeScript("document.body.style.zoom = '50%'")
  # if(SCshot){remDr$screenshot(file = paste0(SCpath,'.png'))}
  
  #-------------------
  #submit
  Sys.sleep(2)
  webElems <- remDr$findElements(using = "id", value = "submitbut")
  webElems[[1]]$clickElement()
  rm(webElems)
  
  #-------------------
  # Wait until page is done to close
  # webElems <- remDr$findElements(using = "class", value = "panel2")
  # curWait <- 0
  
  # while((! webElems[[1]]$isElementDisplayed()[[1]]) && (! startsWith(webElems[[1]]$getElementText()[[1]], "Results\nCreating")) && (curWait < maxWait)){
  #   Sys.sleep(10)
  #   curWait <- curWait + 10
    
  #   webElems <- remDr$findElements(using = "class", value = "panel2")
  # }
  Sys.sleep(15*60)# shortSleep
    
  save(1, file = paste0(SCpath,'.png'))
  #-------------------
  #close
  remDr$close()
  rm(rD)
  gc()

  
}



