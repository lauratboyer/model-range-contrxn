## Theortcl-model_Range-contrxn_Figures-base.r
## Defines plot functions to show basic outputs of spatial
## population model with density-dependent dispersal
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: June 16, 2014
## Time-stamp: <2014-07-14 10:56:32 Laura>

########################################################
# format information on run parameters to include unobstrusively
# in plot.emat() graph
format.run.info <- function(ri=envpop$run.info) {

    ri.tog <- paste(names(ri), ri, sep="=")
    pu <- par("usr")

    gx <- grconvertX(c(0.02,0.99), from="ndc")
    gy <- grconvertY(0.01, from="ndc")
    gwidth <- gx[2]-gx[1]
    srf <- max(strwidth(ri.tog, vfont=c("sans serif","plain")))
    sh <- max(strheight(ri.tog, vfont=c("sans serif","plain")))
    num.str <- floor(gwidth/srf)
    num.lines <- ceiling(length(ri.tog)/num.str)

    text(gx[1],gy+3*sh, "*", offset=0, xpd=NA, pos=4, col="grey50")
    text(gx[1],gy+3*sh, paste(ri.tog[1:num.str], collapse=", "),
         xpd=NA, pos=4, col="grey50")
    text(gx[1],gy+1.5*sh, paste(ri.tog[(num.str+1):(2*num.str)], collapse=", "),
         xpd=NA, pos=4, col="grey50")
    text(gx[1],gy, paste(ri.tog[(2*num.str+1):length(ri.tog)], collapse=", "),
         xpd=NA, pos=4, col="grey50")

}

## 1. Panel plot overview of population abundance by cell over time
## with zoom on last 50 timesteps, map of cell layout + K,
## and option to plot ratio of N/K
plot.emat <- function(emat=envpop$mat, show.nk=FALSE) {

    plines <- function(xl) {

        ymax <- ifelse(show.nk, 1.25, max(K))
        plot(1:10, type="n", xlim=xl, ylim=c(0,ymax), las=1)
        abline(h=K, col=c(col.mat.transp), lwd=0.5)

        cellv <- 1:ncell
        mi <- arrayInd(cellv, .dim=c(grid.width, grid.width))
        dmm <- sapply(cellv, function(rr) lines(emat[mi[rr,1],
                                                     mi[rr,2],]/ifelse(show.nk,K[rr],1),
                                                lwd=2, col=col.mat.transp[rr]))
    }
    check.dev.size(6.4, 7.75)
    par(mai=rep(0.3,4), omi=c(0.65,0.5,0.35,0.5), family="HersheySans")
    layout(rbind(3,c(1,2)),height=c(2,1),width=c(2,1))
    im.mat <- 1: (nrow(emat)*ncol(emat))
    attr(im.mat,"dim") <- dim(emat)[1:2]
    ymax <- max(emat, na.rm=TRUE)
    plines(c(450, 490))

    pmai <- par()$mai
    par(mai=c(0.3,0,0.45,0.3))
    image(1:nrow(emat),1:ncol(emat),
          im.mat, col=c(col.mat),asp=1,axes=FALSE)
    box()
    text(cell.index$xx, cell.index$yy, K, cex=0.5,
         col="white", vfont=c("sans serif","bold"))
    mtext("Cell layout and K", adj=0)
    abline(h=(1:nrow(emat))-0.5)
    abline(v=(1:nrow(emat))-0.5)
    par(mai=pmai)
    plines(c(0, ts.max))
    lab <- sprintf("Baseline emigration rate: %s;
 maximum emigration rate %s", emig.base, emig.max)
    lab2 <- ifelse(sd(K)==0, "No environment heterogeneity",
                   sprintf("min K: %s, max K: %s, sd K: %s",
                           min(K), max(K), round(sd(K),1)))
    mtext(lab,adj=0,line=2)
    mtext(sprintf("Habitat heterogeneity type: %s", habtype), adj=0, line=1)
    mtext(lab2,adj=0)

    fname <- sprintf("Theo-mod_range-contrxn_emigbase-%s_emigmax-%s_habtype-%s",
                     emig.base, emig.max, habtype)
    fname <- gsub("\\.","",fname)

    format.run.info()
    dev.copy2pdf(file=paste(fname,".pdf",sep=""))
}

########################################################################
########################################################################
## Panel plot. Immigrants by cell vs emigrants
## Would be nice to add lower two panels with emigration rate
plot.immig <- function(imat1=envpop$immig.store,
                       emat1=envpop$emig.store, zoom.x) {

    check.dev.size(10,4.5); par(mfrow=c(1,2)); par(family="HersheySans",mai=rep(0.5,4))
    ymax <- max(c(imat1,emat1), na.rm=TRUE)
    xl <- c(0,ts.max)
    if(!missing(zoom.x)) xl <- zoom.x
    plot(0, type="n", xlim=xl, ylim=c(0,ymax), las=1, ann=FALSE)
    dmm <- sapply(1:ncell, function(x) lines(1:ts.max, imat1[,x], col=col.mat[x]))
    mtext("# of migrants to cell")
    plot(0, type="n", xlim=xl, ylim=c(0,ymax), las=1, ann=FALSE)
    dmm <- sapply(1:ncell, function(x) lines(1:ts.max, emat1[,x], col=col.mat[x]))
    mtext(" of emigrants from cell")
}

########################################################################
########################################################################
## Single plot. Ratio of immigrants sent out vs emigrants received by cell
## When >1 cell is source, when <1 cell is sink
plot.ratio.ei <- function(imat1=envpop$immig.store,
                       emat1=envpop$emig.store) {

    check.dev.size(8,5); par(mfrow=c(1,1), mai=c(0.8,1,0.5,1.5)); par(family="HersheySans")
    #imat1[imat1==0] <- 0.1
    #emat1[emat1==0] <- 0.1
    ratio.mat <- emat1/imat1 # if e > i, source, else sink

    ymax <- range(ratio.mat, na.rm=TRUE)
    ymax <- c(0,2)
    plot(0, type="n", xlim=c(0,ts.max), ylim=ymax, las=1, ann=FALSE)
    abline(h=1)
    dmm <- sapply(1:ncell, function(x) lines(1:ts.max, ratio.mat[,x],
                                             lwd=3, col=col.mat[x]))
    mtext("Emigrants/immigrants")
    gy <- grconvertY(c(0.25, 0.75), from="npc")
    pu <- par("usr")
    text(pu[2], gy, c("Sink", "Source"), xpd=NA, pos=4)

    fname <- sprintf("Theo-mod_range-contrxn_EIratio_emigbase-%s_emigmax-%s_habtype-%s",
                     emig.base, emig.max, habtype)
    fname <- gsub("\\.","",fname)
    dev.copy2pdf(file=paste(fname,".pdf",sep=""))
}

########################################################################
########################################################################
## Map of cell abundance over time for defined time-steps
## mark as exctinct if decline over 95%
map.abund <- function(popmat=envpop$mat, K.thresh=0.05, global.scale=TRUE) {

    check.dev.size(6.35,7.75)
    par(family="HersheySans", mfrow=c(7,5),
        omi=c(0,0.25,0.25,0.25),mai=c(0.15,0.2,0.15,0.2))

    dn <- dimnames(popmat)
    tsv <- c(1:9, seq(10, 100, by=10), seq(200, 950, by=50))
    tsv <- tsv[tsv <= dim(envpop$mat)[3]]
    popmat <- popmat/K

    if(global.scale) {
        abund.breaks <- seq(min(popmat), max(popmat), length=21)
        if(abund.breaks[1]==0) {
            abund.breaks <- c(-0.1,abund.breaks)
        }else{
            abund.breaks <- c(-0.1, 0, abund.breaks) }
        colv <- c("dark grey",rev(heat_hcl(length(abund.breaks)-2)))
    }

    map.ts <- function(ts) {

        imat <- popmat[,,ts]
        imat[imat<K.thresh] <- 0

        if(!global.scale) {
        abund.breaks <- seq(min(imat), max(imat), length=21)
        if(abund.breaks[1]==0) {
            abund.breaks <- c(-0.1,abund.breaks)
        }else{
            abund.breaks <- c(-0.1, 0, abund.breaks) }
        colv <- c("dark grey",rev(heat_hcl(length(abund.breaks)-2)))
    }
        image(1:nrow(imat),1:ncol(imat), imat,
          asp=1,axes=FALSE,breaks=abund.breaks,col=colv)
        box()
        mtext(paste("t=",ts,sep=""),adj=0)

        abline(h=(1:nrow(imat))-0.5)
        abline(v=(1:nrow(imat))-0.5)
    }

    dmm <- sapply(tsv, map.ts)
}

########################################################################
########################################################################
plot.cell.layout <- function(emat=envpop$mat) {

    check.dev.size(5,5)
    par(mai=rep(0.3,4), family="HersheySans")
    im.mat <- 1: (nrow(emat)*ncol(emat))
    attr(im.mat,"dim") <- dim(emat)[1:2]
    ymax <- max(emat, na.rm=TRUE)
    image(1:nrow(emat),1:ncol(emat),
          im.mat, col=c(col.mat),asp=1,axes=FALSE,
          xlab="", ylab="")
    box()
    text(cell.index$xx, cell.index$yy,
         paste(paste("#", cell.index$index,sep=""),
               K, r.growth, sep="//"), cex=0.5,
         col="black", vfont=c("sans serif","bold"))
    mtext("Cell layout and K", adj=0)
    abline(h=(1:nrow(emat))-0.5)
    abline(v=(1:nrow(emat))-0.5)
}


########################################################################
########################################################################
plot.tm.aor <- function(popmat=envpop$mat, pres.thresh=0) {

    abund.ts <- apply(popmat, 3, sum)
    area.ts <- apply((floor(popmat))>pres.thresh, 3, sum)

    plot(abund.ts, area.ts)


    invisible(list(abund=abund.ts, area=area.ts))
}

###################################################
###################################################
plot.fish.biom <- function(pmat=envpop$mat) {

    plot(c(F.array[,,ts.max]), c(pmat[,,ts.max]/K), ylim=c(0,1))





}
