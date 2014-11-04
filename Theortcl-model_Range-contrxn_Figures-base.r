## Theortcl-model_Range-contrxn_Figures-base.r
## Defines plot functions to show basic outputs of spatial
## population model with density-dependent dispersal
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: June 16, 2014
## Time-stamp: <2014-10-29 11:02:45 Laura>
########################################################
# Define general labels and such

emig.lab <- c("None", "Uniform", "Preferential")
names(emig.lab) <- c("base","emig","wpref")
ftype.lab <- c("Even F throughout the range","Higher F in core", "Higher F in edges")
names(ftype.lab) <- c("even","core","edge")

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

  if(class(emat)=="environment") emat <- emat$mat

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
plot.cell.layout <- function(add.info=TRUE, use.transp=FALSE) {

  #emat=envpop$mat
    check.dev.size(5,5)
    par(mai=rep(0.3,4), family="HersheySans")
    im.mat <- 1:25#(nrow(emat)*ncol(emat))
    attr(im.mat,"dim") <- c(5,5) #dim(emat)[1:2]
#    ymax <- max(emat, na.rm=TRUE)
#    image(1:nrow(emat),1:ncol(emat),
#          im.mat, col=c(col.mat),asp=1,axes=FALSE,
#          xlab="", ylab="")
    col.now <- col.mat
    if(use.transp) col.now <- col.mat.transp
 image(1:5,1:5,
          im.mat, col=col.now,asp=1,axes=FALSE,
          xlab="", ylab="")
    box()
    if(add.info) {
    text(cell.index$xx, cell.index$yy,
         paste(paste("#", cell.index$index,sep=""),
               K, r.growth, sep="//"), cex=0.5,
         col="black", vfont=c("sans serif","bold"))
    mtext("Cell layout and K", adj=0)
  }
#    abline(h=(1:nrow(emat))-0.5)
#    abline(v=(1:nrow(emat))-0.5)
    abline(h=(1:5)-0.5)
    abline(v=(1:5)-0.5)

    dev.copy2pdf(file="tm-range-dyn_cell-layout-map.pdf")

}


########################################################################
########################################################################
plot.tm.aor <- function(popmat=envpop$mat, pres.thresh=0, add=FALSE,
                        colv="black", pchv=19, ...) {

    abund.ts <- apply(popmat, 3, sum)
    area.ts <- apply((floor(popmat))>(pres.thresh*K), 3, sum)


    if(!add) {
        plot(abund.ts, area.ts, las=1, xlab="", ylab="", xaxt="n",
             col=colv, pch=pchv, ...)

        axis(1, labels=NA)
    }else{        points(abund.ts, area.ts, col=colv, pch=pchv, ...) }

    ## add 1:1 line
    #sl <- max(area.ts)/max(abund.ts)
    #abline(0, sl)
    #points(abund.ts, sl*abund.ts)

    invisible(list(abund=abund.ts, area=area.ts))
}

########################################################################
########################################################################
calc.aor.metric <- function(popmat=envpop$mat, pres.thresh=0.05, add=FALSE, ...) {

    abund.ts <- apply(popmat, 3, sum)
    area.ts <- apply((floor(popmat))>(pres.thresh*K), 3, sum)

    ## get 1:1 line
    sl <- max(area.ts)/max(abund.ts)
    diff <- sum(area.ts - sl*abund.ts)
    sign <- ifelse(diff < 0, -1, 1)

    s1 <- sign*1/log(diff) # metric -- retain sign but log to diminish differences
}

bp.aor.metric <- function() {


  metr.mat <- cbind(sapply(aor.emig, calc.aor.metric),
                    sapply(aor.wpref, calc.aor.metric),
                    sapply(aor.whighpref, calc.aor.metric))
  metr.mat <- metr.mat/metr.mat[1] # standardize to no AOR

  check.dev.size(5.6,6.25)
  par(mfrow=c(1,1), mai=c(0.5,1,0.25,1.2), family="HersheySans")
  layout(rbind(1,2),height=c(1.5,1))

  # top panel raw data for AOR
  pthresh <- 0.05
  plot.tm.aor(popmat=aor.emig$noss, pres.thresh=pthresh, col=NA)
  mtext("Range size in # of squares occupied",side=2,line=3)
  mtext("Total population abundance",side=1, line=1)

  abline(0, 100/1000000, col="grey")
  bpal <- c("grey",col2transp(brewer.pal(3, "YlGnBu")))

  lapply(1:4, function(x) plot.tm.aor(aor.emig[[x]], pres.thresh=pthresh, add=TRUE,
                                           colv=bpal[x]))
  lapply(1:4, function(x) plot.tm.aor(aor.wpref[[x]], pres.thresh=pthresh,
                                            add=TRUE, pchv=3, colv=bpal[x]))
  lapply(1:4, function(x) plot.tm.aor(aor.whighpref[[x]], pres.thresh=pthresh,
                                                add=TRUE, pchv=4, colv=bpal[x]))

  pu <- par("usr")
  legend.ltb.2(1.025*pu[2], pu[4],
               c("None","Slight","Medium","Strong", NA,
                 "Uniform", "Weak", "Strong"),
         col=c(bpal,NA,rep(bpal[3],3)), bty="n", xpd=NA,
               pch=c(rep(15,5),20,3,4), pt.cex=2)

  # bottom panel AOR metric
  par(mai=c(0.5,1,0.1,0.2))
  bpal <- col2transp(brewer.pal(3, "YlGnBu"),1)
  bp <- barplot(metr.mat, beside=TRUE, col=c("grey",bpal),
          border=NA, las=1, xlim=c(0,9),ylab="AOR strength metric",
          space=c(0.5,rep(0.2,3), 1.5, rep(0.2, 3), 1.5, rep(0.2,3)), width=0.5)
  abline(h=0)
  xpos <- apply(bp, 2, mean)
  text(xpos, -0.1, c("Uniform",
                     "Weak preferential",
                     "Strong preferential"), xpd=NA, cex=0.8)


  dev.copy2pdf(file="tm-range-contrxn_aor-scenario_metric-bp.pdf")

}
###################################################
###################################################
plot.fish.biom <- function(pmat=envpop$mat) {

    plot(c(F.array[,,ts.max]), c(pmat[,,ts.max]/K), ylim=c(0,1))
}

###################################################
###################################################
## Plot scenario time-series (3) (i.e. top panel of emat)
plot.scen.ts <- function(scen.list=ee.lowF, show.nk=FALSE) {

    Know <- attr(scen.list, "K")
    rgval <- attr(scen.list, "rgval")
    grid.width <- attr(scen.list, "grid.width")
    hab.struct <- ifelse(rgval[1]!=rgval[2], "Core-edge", "Even")

    plines <- function(wtype, xl, F.type) {

      # scenarios to show in plot:
      scen.sub <- scen.list[[F.type]]

      emat <- scen.sub[[wtype]]$mat
      ymax <- ifelse(show.nk, 1.25, max(Know))
      plot(1:10, type="n", xlim=xl, ylim=c(0,ymax), las=1,
           xaxt=ifelse(F.type=="edge","t","n"),
             yaxt=ifelse(wtype=="base", "t", "n"))
      abline(h=K*r.growth/r.mrt, col=c(col.mat.transp), lwd=0.5)

        cellv <- 1:ncell
        mi <- arrayInd(cellv, .dim=c(grid.width, grid.width))
        dmm <- sapply(cellv, function(rr) lines(emat[mi[rr,1],
                                                     mi[rr,2],]/ifelse(show.nk,K[rr],1),
                                                lwd=2, col=col.mat.transp[rr]))
    pu <- par("usr")
    lh <- strheight("I", vfont=c("sans serif","plain")) # line height for labels

      if(F.type=="even") text(mean(pu[1:2]), pu[4] + 5*lh,
           toupper(emig.lab[wtype]), col="royalblue4",
           vfont=c("sans serif","bold"), xpd=NA, cex=1.5)

      if(wtype=="base") text(pu[1],pu[4] + 1*lh,
           ftype.lab[F.type],xpd=NA, pos=4, cex=1.5,
           vfont=c("sans serif", "bold"), offset=0, col="royalblue3")
    }

    check.dev.size(9.3, 8)
    par(mfrow=c(3,3), mai=c(0.1,0.1,0.5,0.1), omi=c(0.65,0.75,0.65,0.1), family="HersheySans")

    sapply(names(scen.list), function(ft) {
      plines("base", c(0, ts.max), ft);
      plines("emig", c(0, ts.max), ft);
      plines("wpref", c(0, ts.max), ft)})

    mtext("Population size (# of individuals)", side=2, outer=TRUE, line=4)
    mtext("Timesteps", side=1, outer=TRUE, line=3)


    # General plot label
    px <- grconvertX(0, from="nic")
    py <- grconvertY(0.95, from="ndc")

    text(px, py, sprintf("%s habitat structure", hab.struct),
         pos=4, cex=2.25, xpd=NA)

    fname <- sprintf("tm-range-contrxn_ts-3-scen_HS-%s_FMx3.pdf",
                     hab.struct)
    dev.copy2pdf(file=fname)
}

# plot time-series of population biomass comparing two set-ups of environment differentiation
plot.scen.ts.comp2 <- function(scen.list=ee.lowF,
                               scen.list.2=ce.lowF, show.nk=FALSE) {

    ts.fstart <- attr(scen.list, "fishing.start")
    grid.width <- attr(scen.list, "grid.width")

    plines <- function(scen.now, wtype, xl) {

      F.type <- "even"
      scen.sub <- scen.now[[F.type]]
      Know <- attr(scen.now, "K")
      rgval <- attr(scen.now, "rgval")
      emat <- scen.sub[[wtype]]$mat

      ymax <- ifelse(show.nk, 1.25, 1.25*max(Know))
      plot(1:10, type="n", xlim=xl, ylim=c(0,ymax), las=1,
           xaxt=ifelse(counter==2,"t","n"),
             yaxt=ifelse(wtype=="base", "t", "n")
           )



      pu <- par("usr")
      polygon(c(pu[1],ts.fstart, ts.fstart,pu[1]), pu[c(3,3,4,4)],
              border=NA, col=col2transp("cornsilk3"))
        abline(h=Know[1]*rgval/r.mrt[1], lwd=1,
               col=c("tomato","dodgerblue3"))
      abline(h=0,col="grey") # horizontal line at n=0
        cellv <- 1:ncell
        mi <- arrayInd(cellv, .dim=c(grid.width, grid.width))

        dmm <- sapply(cellv, function(rr) lines(emat[mi[rr,1],
                                                     mi[rr,2],]/ifelse(show.nk,K[rr],1),
                                                lwd=2, col=col.mat.transp[rr]))
    pu <- par("usr")
    lh <- strheight("I", vfont=c("sans serif","plain")) # line height for labels
        mtext(toupper(emig.lab[wtype]),adj=1,col="royalblue")
      if(wtype=="base") mtext(toupper("Dispersal type:"),adj=0,col="grey30")

      if(wtype=="base" & counter ==1) mtext("Low difference in habitat quality",adj=0,line=2,cex=1.5)
      if(wtype=="base" & counter ==2) mtext("High difference in habitat quality",adj=0,line=2,cex=1.5)
      text(ts.fstart,pu[4]-1.75*lh,"F=0",col="white",cex=1.5,pos=2,vfont=c("sans serif", "bold"))

#      if(wtype=="base") text(pu[1],pu[4] + 1*lh,
#           ftype.lab[F.type],xpd=NA, pos=4, cex=1.5,
#           vfont=c("sans serif", "bold"), offset=0, col="royalblue3")
    }

    check.dev.size(9.3, 7.5)
    par(mfrow=c(2,3), mai=c(0.1,0.1,0.75,0.1), omi=c(0.65,0.75,0.15,0.1), family="HersheySans")

    counter <- 1
      sapply(c("base","emig","wpref"),
           function(ft) plines(scen.list, ft, c(0, ts.max)))

    counter <- 2
    sapply(c("base","emig","wpref"),
           function(ft) plines(scen.list.2, ft, c(0, ts.max)))

    mtext("Population size (# of individuals)", side=2, outer=TRUE, line=4)
    mtext("Timesteps", side=1, outer=TRUE, line=3)


    # General plot label
    px <- grconvertX(0, from="nic")
    py <- grconvertY(0.95, from="ndc")

#    text(px, py, sprintf("%s habitat structure", "lala"),
#         pos=4, cex=2.25, xpd=NA)

    fname <- "tm-range-contrxn_ts-3-scen_FMx3_two-envir-scenarios.pdf"
    dev.copy2pdf(file=fname)
}

########################################
## Proportion of population in core cells
prop.in.core <- function() {

  b4fish <- 250 # time-steps at equilibrium
  aftfish <- 900
  cc <- calc.core.size(100)$core.cells

  calc.dprop <- function(arr.now) {

  coreb4 <- sum(arr.now[,,b4fish][cc])/length(b4fish)
  allb4 <-  sum(arr.now[,,b4fish])/length(b4fish)

  coreaf <- sum(arr.now[,,aftfish][cc])/length(aftfish)
  allaf <-  sum(arr.now[,,aftfish])/length(aftfish)

  c(coreb4/allb4, coreaf/allaf) # difference in proportion of biomass in core cells
}

  r.emig <- sapply(aor.emig, calc.dprop)
  r.wpref <- sapply(aor.wpref, calc.dprop)
  r.whpref <- sapply(aor.whighpref, calc.dprop)
  mat <- cbind(r.emig, r.wpref, r.whpref)

  par(mfrow=c(1,1),mai=c(0.8,0.8,0.8,2),family="HersheySans")
  plot(1:10,type="n",ylim=c(0,11.5),xlim=c(0,1), axes=FALSE, ann=FALSE)
  abline(v=seq(0,0.8,by=0.2),col="cornsilk3")
  bwidth <- 0.3
  bp <- barplot(mat, beside=TRUE, border=NA,
          col=c("limegreen","tomato"), las=1,
          width=bwidth, space=c(rep(c(1, 0),4), 2,0,rep(c(1, 0),3),
                        2,0,rep(c(1,0),3)),
                horiz=TRUE, add=TRUE)
  text(1, bp[8], "Uniform dispersal",xpd=NA,pos=4)
  text(1, bp[16], "Weak preferential",xpd=NA,
       pos=4)
  text(1, bp[15], "dispersal",xpd=NA,
       pos=4)
  text(1, bp[24], "Strong preferential",xpd=NA,
       pos=4)
  text(1, bp[23], "dispersal",xpd=NA,
       pos=4)
  segments(1.01,bp[24]+bwidth/2,1.01,bp[17],xpd=NA,lwd=0.5)
  segments(1.01,bp[16]+bwidth/2,1.01,bp[9],xpd=NA,lwd=0.5)
  segments(1.01,bp[8]+bwidth/2,1.01,bp[1],xpd=NA,lwd=0.5)
}
