################################
################################
# chow
# Chow test for structural break
# and test unbroken slope change
################################
################################
chow = function(t,x,tbeg=NULL,tend=NULL,plot=TRUE,buffer=4){
    buffer = floor(buffer)
    if (buffer < 4){
        buffer = 4
    }
    ########################
    # remove missing entries
    ########################
    ztemp = data.frame(t,x)
    ztemp = na.omit(ztemp)
    if (length(tbeg) == 1){
        ztemp = ztemp[ztemp$t >= tbeg,]
    }
    if (length(tend) == 1){
        ztemp = ztemp[ztemp$t <= tend,]
    }
    t = ztemp$t; x = ztemp$x
    #######
    # do it
    #######
    ntot = length(t)
    if (ntot < 8) stop("Not enough data!")
    ndx = (buffer+1):(ntot-buffer)
    numndx = length(ndx)
    Ftest = rep(0,numndx)    # Chow test
    pvals = numeric(numndx)  # Chow test
    F1 = rep(0,numndx)       # unbroken test
    p1 = rep(0,numndx)       # unbroken test
    xfit = lm(x~t)
    Stot = sum(xfit$res^2)
    this.ndx = 0
    for (jbreak in ndx){
        ###########
        # Chow test
        ###########
        t1 = t[1:(jbreak-1)]
        x1 = x[1:(jbreak-1)]
        t2 = t[jbreak:ntot]
        x2 = x[jbreak:ntot]
        xfit1 = lm(x1~t1)
        xfit2 = lm(x2~t2)
        S1 = sum(xfit1$res^2)
        S2 = sum(xfit2$res^2)
        this.ndx = this.ndx+1
        ff = (Stot-S1-S2)*(ntot-4)/(S1+S2)/2
        pp = pf(ff,2,ntot-4,lower.tail=FALSE)
        Ftest[this.ndx] = ff
        pvals[this.ndx] = pp
        ###############
        # unbroken test
        ###############
        phi = t-t[jbreak]
        phi[phi < 0] = 0
        xfit = lm(x~t+phi)
        Sub = sum(xfit$res^2)
        ff1 = (Stot-Sub)*(ntot-3)/Sub
        pp1 = pf(ff1,1,ntot-3,lower.tail=FALSE)
        F1[this.ndx] = ff1
        p1[this.ndx] = pp1
    }
    Ftest[!is.finite(Ftest)] = 0
    F1[!is.finite(F1)] = 0
    pmin = min(pvals)
    ndx.chow = which(pvals == pmin)
    pmin1 = min(p1)
    ndx.unbroken = which(p1 == pmin1)
    tmax.chow = t[ndx[ndx.chow]]
    tmax.unbroken = t[ndx[ndx.unbroken]]
    ###########
    # best fits
    ###########
        ########
        # broken
        ########
        t1 = t[t < tmax.chow]; x1 = x[t < tmax.chow]
        t2 = t[t >= tmax.chow]; x2 = x[t >= tmax.chow]
        x1lin = lm(x1~t1); x2lin = lm(x2~t2)
        ##########
        # unbroken
        ##########
        phi = t-tmax.unbroken
        phi[phi < 0] = 0
        xulin = lm(x~t+phi)
    #############
    # output data
    #############
    zout = list(tmax.chow=tmax.chow,
        tmax.unbroken=tmax.unbroken,
        pmin.chow=pmin,pmin.unbroken=pmin1,
        ndx=ndx,t=t[ndx],
        Ftest.chow=Ftest,pval.chow=pvals,
        Ftest.unbroken=F1,pval.unbroken=p1,
        fit.chow=c(x1lin$fit,x2lin$fit),fit.unbroken=xulin$fit)
    if (plot){
        ######
        # plot
        ######
        plot(t,x,type="o",pch=20,cex=1.2,col="cornsilk4",
            xlab="",ylab="")
        lines(t1,x1lin$fit,lwd=3,col="blue")
        lines(t2,x2lin$fit,lwd=3,col="blue")
        lines(t,xulin$fit,lwd=3,col="red")
    }
    cat("minimum p-value\n")
    cat(pmin,"   Chow test        at",tmax.chow,"\n")
    cat(pmin1,"   Unbroken trend   at",tmax.unbroken,"\n\n")
    invisible(zout)
}

montechow = function(nsize,ntimes=500){
    t = 0:(nsize-1)
    pchow = numeric(ntimes)
    punbr = numeric(ntimes)
    count = 0
    for (j in 1:ntimes){
        x = rnorm(nsize)
        q = chow(t,x,plot=F)
        punbr[j] = q$pmin.unbroken
        pchow[j] = q$pmin.chow
        count = count+1
        if (count > 9){
            count = 0
            plot(0,0,main=paste(j,"/",ntimes))
        }
    }
    pchow = sort(pchow)
    punbr = sort(punbr)
    tp = (1:ntimes)/(ntimes+1)
    zout = data.frame(true.pval=tp,chow.pval=pchow,unbroken.pval=punbr)
    zout = round(zout,6)
    invisible(zout)
}
