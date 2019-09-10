library(ggplot2)
library(reshape)
#-------------------------------------------------
# inputs
ifolder <- '/home/wanda/Dropbox/20190816_SimulationVar/SimulatedTSStab/'
setfolder <- '/home/wanda/Dropbox/20190816_SimulationVar/SimulatedTS/'
figfolder <- '/home/wanda/Dropbox/20190816_SimulationVar/Figures'
ecoreg <- 567
basename <- c('S1_Sample_EcoReg_','_noDist_Tree_50_lossYr_18_scl_30_npnt_1000_DESCENDING_10_H_IW_VH')
pol <- 'VH'
simcase <- 'len'# dr'distRec distMag distT len seasAmp remSd

#-------------------------------------------------
# functions

# Load RData with user specified name
loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
}

#-------------------------------------------------
# import data
# rmse <- loadRData(file=file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_RMSE_', simcase, '.rda')))
# mape <- loadRData(file=file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MAPE_', simcase, '.rda')))
meas_bfast <- loadRData(file=file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASbfastOptim_', simcase, '.rda')))
meas_frazier <- loadRData(file=file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASfrazier_', simcase, '.rda')))
meas_rgrowth <- loadRData(file=file.path(ifolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_MEASrgrowth_', simcase, '.rda')))
sttngs <- loadRData(file=file.path(setfolder, paste0(basename[1], ecoreg, basename[2], '_', pol, '_simTS_settings.rda')))

simnames <- list('Disturbance magnitude [dB]',
                 'Number of droughts [-]',
                 'Time series length [year]',
                 'Seasonal amplitude [dB]',
                 'SD remainder [dB]',
                 'Disturbance timing [year]',
                 'Recovery period [year]')
names(simnames) <- c('distMag',
                     'dr', 
                     'len',
                     'seasAmp',
                     'remSd',
                     'distT',
                     'distRec')

lvls <- sttngs[[2]]
    
names(lvls) <- c('dr', 
                 'len',
                 'seasAmp',
                 'remSd',
                 'distMag',
                 'distT',
                 'distRec')
lvls$seasAmp <- lvls$seasAmp*2
lvls$distRec <- lvls$distRec/365
lvls$seasAmp <- round(lvls$seasAmp,2)
lvls$remSd <- round(lvls$remSd,2)

#-------------------------------------------------
# missing values
 
countMissVal <- function(meas, levels, method){
        vl <- colSums(is.na(meas)==T)/dim(meas)[1]*100
        tmp <- data.frame(t(vl))
        names(tmp) <- levels
        tmp2 <- melt(tmp)
        tmp2$Method <- factor(method) 
        tmp2
}

missVal <- countMissVal(meas_frazier$m_RRI, lvls[[simcase]], 'RRI')
missVal <- rbind(missVal, countMissVal(meas_frazier$m_RRIsm, lvls[[simcase]], 'RRIsmooth'))
missVal <- rbind(missVal, countMissVal(meas_frazier$m_R80p, lvls[[simcase]], 'R80p'))
missVal <- rbind(missVal, countMissVal(meas_frazier$m_R80psm, lvls[[simcase]], 'R80psmooth'))
missVal <- rbind(missVal, countMissVal(meas_frazier$m_YrYr, lvls[[simcase]], 'YrYr'))
missVal <- rbind(missVal, countMissVal(meas_frazier$m_YrYrsm, lvls[[simcase]], 'YrYrsmooth'))
missVal <- rbind(missVal, countMissVal(meas_rgrowth$m_lrec, lvls[[simcase]], 'Recov_period_BFmon'))
missVal <- rbind(missVal, countMissVal(meas_rgrowth$m_lrecRec, lvls[[simcase]], 'Recov_period'))
missVal <- rbind(missVal, countMissVal(meas_bfast$m_SLbf, lvls[[simcase]], 'BF_slope'))
missVal <- rbind(missVal, countMissVal(meas_bfast$m_RRIbf, lvls[[simcase]], 'BF_RRI'))
missVal <- rbind(missVal, countMissVal(meas_bfast$m_R80pbf, lvls[[simcase]], 'BF_R80p'))
missVal <- rbind(missVal, countMissVal(meas_bfast$m_YrYrbf, lvls[[simcase]], 'BF_YrYr'))


# vl <- colSums(is.na(meas_frazier$m_RRI)==T)/dim(meas_frazier$m_RRI)[1]*100
# tmp <- data.frame(t(vl))
# names(tmp) <- lvls[[simcase]]
# tmp2 <- melt(tmp)
# tmp2$Method <- factor('RRI')
# 
# df <- rbind(tmp2,tmp3)

png(file.path(figfolder, paste0('PercMV_', simcase)),
    width = 837, height = 345, units = "px",)
ggplot(missVal,aes(variable,value,fill=Method))+
        geom_bar(stat="identity",position='dodge') +  coord_cartesian(ylim = c(-0.25, 100))+
        scale_fill_manual(values=c("#F32F1B", "#904325", "#FCBD1C", "#B88C11",
                                   "#83F665","#60B649", "#79ACFC","#476A9C",
                                   "#FB80CA", "#7B3C63","#F6685C","#86352E"))+ 
        xlab(simnames[[simcase]]) +
        ylab('Percentage missing values [%]')
dev.off()

#-------------------------------------------------
# errors
library(reshape)
# RRI

prepDFAbsError <- function(meas, val, levels, method){
        tmp <- data.frame(abs(meas - val))
        names(tmp) <- levels
        tmp2 <- melt(tmp)
        tmp2$Method <- factor(method)
        tmp2
}
prepDFAbsPercError <- function(meas, val, levels, method){
        tmp <- data.frame(abs(meas - val)/abs(val)*100)
        names(tmp) <- levels
        tmp2 <- melt(tmp)
        tmp2$Method <- factor(method)
        tmp2
}

absPercErrors <- prepDFAbsPercError(meas_frazier$m_RRI, meas_frazier$s_RRI, lvls[[simcase]], 'RRI')
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_frazier$m_RRIsm, meas_frazier$s_RRI, lvls[[simcase]], 'RRIsmooth'))
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_frazier$m_R80p, meas_frazier$s_R80p, lvls[[simcase]], 'R80p'))
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_frazier$m_R80psm, meas_frazier$s_R80p, lvls[[simcase]], 'R80psmooth'))
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_frazier$m_YrYr, meas_frazier$s_YrYr, lvls[[simcase]], 'YrYr'))
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_frazier$m_YrYrsm, meas_frazier$s_YrYr, lvls[[simcase]], 'YrYrsmooth'))
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_rgrowth$m_lrec, meas_rgrowth$s_lrec, lvls[[simcase]], 'Recov_period_BFmon'))
absPercErrors <- rbind(absPercErrors, prepDFAbsPercError(meas_rgrowth$m_lrecRec, meas_rgrowth$s_lrec, lvls[[simcase]], 'Recov_period'))

png(file.path(figfolder, paste0('APE_', simcase)),
    width = 837, height = 345, units = "px",)
ggplot(absPercErrors, aes(variable, value, fill=Method)) +
        geom_boxplot() + coord_cartesian(ylim = c(-0.25, 100)) + 
        scale_fill_manual(values=c("#F32F1B", "#904325", "#FCBD1C", "#B88C11",
                                   "#83F665","#60B649", "#79ACFC","#476A9C",
                                   "#FB80CA", "#7B3C63","#F6685C","#86352E"))+ 
        xlab(simnames[[simcase]]) +
        ylab('Absolute Percentage Error [%]')
dev.off()

png(file.path(figfolder, paste0('APEmean_', simcase)),
    width = 837, height = 345, units = "px",)
ggplot(absPercErrors, aes(variable, value)) +
        geom_boxplot()  + coord_cartesian(ylim = c(-0.25, 100))+
        xlab(simnames[[simcase]]) +
        ylab('Absolute Percentage Error [%]')        
dev.off()


#-------------------------------------------------
# R2
# rsq <- function(x, y) summary(lm(y~x))$r.squared
# rsq(m_tdst[,1], s_tdst[,1])
# plot(meas_frazier$m_R80psm[,1],meas_frazier$s_R80p[,1])
rsq <- function(x, y) {
    if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3)){
        rs <- summary(lm(y~x))$r.squared
    }else(rs <- NA)
    rs
}

calcRsq <- function(val, meas, levels, method){
    vl <- sapply(1:dim(meas)[2], function(i) rsq(val[,i], meas[,i]))
    tmp <- data.frame(t(vl))
    names(tmp) <- levels
    tmp2 <- melt(tmp)
    tmp2$Method <- factor(method) 
    tmp2
}

dfRsq <- calcRsq(meas_frazier$s_RRI, meas_frazier$m_RRI, lvls[[simcase]], 'RRI')
dfRsq <- rbind(dfRsq, calcRsq(meas_frazier$s_RRI, meas_frazier$m_RRIsm, lvls[[simcase]], 'RRIsmooth'))
dfRsq <- rbind(dfRsq, calcRsq(meas_bfast$m_RRIbf, meas_bfast$s_RRIbf,lvls[[simcase]], 'RRI_BF'))
dfRsq <- rbind(dfRsq, calcRsq(meas_frazier$s_R80p,meas_frazier$m_R80p,lvls[[simcase]], 'R80p'))
dfRsq <- rbind(dfRsq, calcRsq(meas_frazier$s_R80p, meas_frazier$m_R80psm, lvls[[simcase]], 'R80psmooth'))
dfRsq <- rbind(dfRsq, calcRsq(meas_bfast$m_R80pbf, meas_bfast$s_R80pbf,lvls[[simcase]], 'R80p_BF'))
dfRsq <- rbind(dfRsq, calcRsq(meas_frazier$s_YrYr, meas_frazier$m_YrYr, lvls[[simcase]], 'YrYr'))
dfRsq <- rbind(dfRsq, calcRsq(meas_frazier$s_YrYr,meas_frazier$m_YrYrsm, lvls[[simcase]], 'YrYrsmooth'))
dfRsq <- rbind(dfRsq, calcRsq(meas_bfast$m_YrYrbf, meas_bfast$s_YrYrbf,lvls[[simcase]], 'YrYr_BF'))
dfRsq <- rbind(dfRsq, calcRsq(meas_bfast$m_SLbf, meas_bfast$s_SLbf,lvls[[simcase]], 'slope_BF'))
dfRsq <- rbind(dfRsq, calcRsq(meas_rgrowth$s_lrec, meas_rgrowth$m_lrec, lvls[[simcase]], 'Rec_period_BF'))
dfRsq <- rbind(dfRsq, calcRsq(meas_rgrowth$s_lrec, meas_rgrowth$m_lrecRec,lvls[[simcase]], 'Rec_period'))

png(file.path(figfolder, paste0('Rsq_', simcase)),
    width = 837, height = 345, units = "px",)
ggplot(dfRsq, aes(variable,value,fill=Method))+
    geom_bar(stat="identity",position='dodge') +  coord_cartesian(ylim = c(0, 1))+
    scale_fill_manual(values=c("#F32F1B", "#904325", "#FCBD1C", "#B88C11",
                               "#83F665","#60B649", "#79ACFC","#476A9C",
                               "#FB80CA", "#7B3C63","#F6685C","#86352E"))+ 
    xlab(simnames[[simcase]]) +
    ylab('R²')
dev.off()

cls <- c("#381814","#5A2C3A","#694B6B","#597499","#219EB0","#1BC5A7","#81E686","#EAFC64")
png(file.path(figfolder, paste0('Rsq_meth_', simcase)),
    width = 1265, height = 345, units = "px",)
ggplot(dfRsq, aes(Method,value,fill=variable))+
    geom_bar(stat="identity",position='dodge') +  coord_cartesian(ylim = c(0, 1))+
    scale_fill_manual(values=cls[round(seq(1,length(cls),length.out=length(lvls[[simcase]])))])+ 
    xlab('Method') +
    ylab('R²')+
    guides(fill=guide_legend(title=simnames[[simcase]]))
dev.off()