## Theortcl-model_Range-contrxn_Scenarios.r
## Setting up scenarios to compare model runs
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 14, 2014
## Time-stamp: <2015-09-17 08:33:34 lauratb>

########################################################
run.df.scen <- FALSE # set to TRUE to run basic scenarios
aor.scenario.rerun <- FALSE # set to TRUE to run AOR scenarios
########################################################
args.eelowF <- list(evenfish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.5, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.5, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.5, emig.base=0.1, pref.disp=1, add.r=TRUE)),
                    corefish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.5, fish.fact.edge=0.1, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.5, fish.fact.edge=0.1, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.5, fish.fact.edge=0.1, emig.base=0.1, pref.disp=1, add.r=TRUE)),
                    edgefish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.1, fish.fact.edge=0.5, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.1, fish.fact.edge=0.5, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.1, fish.fact=0.1, fish.fact.edge=0.5, emig.base=0.1, pref.disp=1, add.r=TRUE)))

args.celowF <- list(evenfish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.5, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.5, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.5, emig.base=0.1, pref.disp=1, add.r=TRUE)),
                    corefish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.5, fish.fact.edge=0.1, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.5, fish.fact.edge=0.1, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.5, fish.fact.edge=0.1, emig.base=0.1, pref.disp=1, add.r=TRUE)),
                    edgefish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.1, fish.fact.edge=0.5, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.1, fish.fact.edge=0.5, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.05, fish.fact=0.1, fish.fact.edge=0.5, emig.base=0.1, pref.disp=1, add.r=TRUE)))


args.ce.low.lowF <- list(evenfish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.085, fish.fact=0.5, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.085, fish.fact=0.5, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.085, fish.fact=0.5, emig.base=0.1, pref.disp=1, add.r=TRUE)))
args.ce.high.lowF <- list(evenfish=list(nodisp=list(r.growth.core=0.1, r.growth.edge=0.02, fish.fact=0.5, emig.base=0),
                        evdisp=list(r.growth.core=0.1, r.growth.edge=0.02, fish.fact=0.5, emig.base=0.1, pref.disp=0),
                        prefdisp=list(r.growth.core=0.1, r.growth.edge=0.02, fish.fact=0.5, emig.base=0.1, pref.disp=1, add.r=TRUE)))


scenario.run <- function(rgcore = 0.1, rgedge=rgcore, Fval=0.5) {


## Constant environment
## fishing everywhere
  F.remove <- 0.1
    evenfish <- list()
    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, emig.base=0)
    evenfish$base <- envpop

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, emig.base=0.1)
    evenfish$emig <- envpop

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, emig.base=0.1, pref.disp=1, add.r.pref=TRUE)
    evenfish$wpref <- envpop

## ... with fishing only in core
    corefish <- list()
    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, fish.fact.edge=F.remove, emig.base=0)
    corefish$base <- envpop

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, fish.fact.edge=F.remove, emig.base=0.1)
    corefish$emig <- envpop

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, fish.fact.edge=F.remove, emig.base=0.1,
               pref.disp=1, add.r.pref=TRUE)
    corefish$wpref <- envpop

    ## ... with fishing only in edge
    edgefish <- list()
    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=F.remove, fish.fact.edge=Fval, emig.base=0)
    edgefish$base <- envpop

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=F.remove, fish.fact.edge=Fval, emig.base=0.1)
    edgefish$emig <- envpop

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=F.remove, fish.fact.edge=Fval, emig.base=0.1,
               pref.disp=1, add.r.pref=TRUE)
    edgefish$wpref <- envpop

    scenl <- list(even=evenfish, core=corefish, edge=edgefish)
    attr(scenl, "ts.max") <- ts.max
    attr(scenl, "Fval") <- Fval
    attr(scenl, "K") <- K
    attr(scenl, "rgval") <- c(core=rgcore, edge=rgedge)
    attr(scenl, "grid.width") <- grid.width
    attr(scenl, "core.layout") <- core.layout
    attr(scenl, "fishing.start") <- envpop$run.info$fishing.start
    return(scenl)
}

## For each scenario calculate the different in biomass
## expected compared to no dispersal version
if(run.df.scen) {
if(!exists("ee.lowF")) ee.lowF <- scenario.run()
#if(!exists("ee.highF")) ee.highF <- scenario.run(Fval=0.8)
if(!exists("ce.lowF")) ce.lowF <- scenario.run(rgedge=0.05)
if(!exists("ce.low.lowF")) ce.low.lowF <- scenario.run(rgedge=0.085)
if(!exists("ce.high.lowF")) ce.high.lowF <- scenario.run(rgedge=0.02)
#if(!exists("ce.highF")) ce.highF <- scenario.run(rgedge=0.1, Fval=0.8)
}

comp.scen.ec <- function(sim.list=ee.lowF, scen="emig",
                         calc.type="absolute") { #scen is emig or wpref

    tsm <- attr(sim.list, "ts.max")
    var.in <- expand.grid(F.type=c("even","core","edge"),
                          region.focus=c("core","edge"),
                          stringsAsFactors=FALSE)

    getstat <- function(ri) { # wregion is core or edge

        F.type <- var.in[ri,1]
        wregion <- var.in[ri,2]
        sl.scen <- sim.list[[F.type]][[scen]]$mat[,,tsm]
        sl.base <- sim.list[[F.type]]$base$mat[,,tsm]

        varbl <- paste(wregion, ".cells", sep="")
        if(calc.type=="absolute") {
            valnow <- ((sl.scen-sl.base)/K)[core.layout[[varbl]]]
        }else{
            valnow <- ((sl.scen/sl.base))[core.layout[[varbl]]]-1
        }

          c(scen.type=paste(paste(F.type,"F",sep=""), varbl, sep="."),
            mean=mean(valnow), sd=sd(valnow),
            min=min(valnow), max=max(valnow))
    }

    out <- as.data.frame(t(sapply(1:nrow(var.in), getstat)))

    class(out$mean) <- "numeric"
    class(out$sd) <- "numeric"
    class(out$min) <- "numeric"
    class(out$max) <- "numeric"
    out

}

######################################################
######################################################
scen.sourcesink <- function(sim.list=ee.lowF,
                            output.type="stats", scen="emig") { #scen is emig or wpref

    tsm <- attr(sim.list, "ts.max")
    var.in <- expand.grid(F.type=c("even","core","edge"),
                          region.focus=c("core","edge"),
                          stringsAsFactors=FALSE)

    getstat <- function(ri) { # wregion is core or edge

        F.type <- var.in[ri,1]
        wregion <- var.in[ri,2]
        immig.eq <- sim.list[[F.type]][[scen]]$immig.store[tsm,]
        emig.eq <- sim.list[[F.type]][[scen]]$emig.store[tsm,]

        varbl <- paste(wregion, ".cells", sep="")
        valnow <- (emig.eq/immig.eq)[core.layout[[varbl]]]-1

        if(output.type == "stats") {
          c(scen.type=paste(paste(F.type,"F",sep=""), varbl, sep="."),
          c(mean=mean(valnow), sd=sd(valnow)))
        }else{
          valnow }
    }

    if(output.type == "stats") {
      out <- as.data.frame(t(sapply(1:nrow(var.in), getstat)))
      class(out$mean) <- "numeric"
      class(out$sd) <- "numeric"
    }else{
      out <- lapply(1:nrow(var.in), getstat)
      lout <- sapply(out, length)
      out <- data.frame(F.type=rep(var.in[,1],lout),
                 region=rep(var.in[,2], lout),
                 ssval=unlist(out))
      out$treatment <- paste(out$F.type,"F.",out$region,sep="")
    }
    out

}

## make custom boxplot showing the mean +/- SD with min max as lines
custom.bp <- function(wdat, rvar, wfact) {

  dsplit <- split(wdat[,rvar], wdat[,wfact]) # split by treatement
  ordr <- 1:3; names(ordr) <- c("evenF","coreF","edgeF")
  trt.val <- ordr[gsub("(.*)\\..*", "\\1",names(dsplit))]
  dsplit <- dsplit[names(dsplit)[order(trt.val)]]


  bxp.stats <- function(x) { # calculate stats for each box

    mn <- mean(x)
    sdx <- sd(x)
    vals <- c(min(x), mn - sdx, mn, mn + sdx, max(x))
    # adjust in case min/max is greater/smaller than mean +/- SD
    if(vals[1]>vals[2]) vals[2] <- vals[1] #
    if(vals[5]<vals[4]) vals[4] <- vals[5] #
    vals
  }

  n <- sapply(dsplit, length)
  list(stats=sapply(dsplit, bxp.stats), n=n)
}



######################################################
######################################################

bp.scenario <- function(wmetric="Bdecl", scen.list=list(ee.lowF, ce.lowF)) {

  if(!(wmetric %in% c("Bdecl", "SrcSink"))) stop()

    check.dev.size(7, 6)
    par(family="HersheySans", mfrow=c(2,1),
        omi=c(0.1,0.75,0.4,1.2), mai=c(0.1,0.1,0.1,0.2))
    Fval <- attr(scen.list[[1]], "Fval")

   draw.panel <- function(wscen) {

     counter <<- counter + 1
     rgval <- attr(wscen, "rgval")
     hab.struct <- ifelse(rgval[1]!=rgval[2], "Core-edge", "Even")

     # select function to compute metric to show
     funk2use <- ifelse(wmetric=="SrcSink", "scen.sourcesink",
                                        "comp.scen.ec")# option: comp.scen.ec, scen.sourcesink

     scen.emig <- do.call(funk2use, list(wscen, scen="emig", output.type="all"))
    scen.wpref <- do.call(funk2use, list(wscen, scen="wpref", output.type="all"))

     bp.emig <- custom.bp(scen.emig, "ssval", "treatment")
     bp.wpref <- custom.bp(scen.wpref, "ssval", "treatment")

     houp <- list(stats=cbind(bp.emig$stats, bp.wpref$stats),
          n=c(bp.emig$n, bp.wpref$n))
    spacev <- c(0.1, rep(c(0.1,1),2),0.1,3,rep(c(0.1,1),2),0.1)

     yl <- c(-35,35)
     if(funk2use == "scen.sourcesink") yl <- c(-0.42,0.42)

     bpval <- bxp(houp, ylim=yl, las=1, boxfill=col2transp(c("tomato1","royalblue1")),
         show.names=FALSE, whisklty=1, staplecol=NA,
                  at=c(1:6, 10:15))

#     bxp(bp.wpref, ylim=yl, las=1, col=c("tomato1","royalblue1"))

#    bp <- barplot(comp.mat,ylim=yl, beside=TRUE,
#                  las=1, col=NA,
#                  space=spacev, xlim=c(0,20), border=NA)
   abline(h=0, col="grey")
#   bp <- barplot(comp.mat,ylim=yl, beside=TRUE, las=1, add=TRUE,,
#              col=c(rep(c("tomato1","royalblue1"),3),
#              rep(c("tomato2","royalblue2"),3)),
#              space=spacev, xlim=c(0,20), density=c(rep(50,6), rep(NA,6)),
#              border=NA)
     box()
#abline(v=mean(bp[2,])-0.5)


   ## Add labels
     add.lab <- FALSE
     if(add.lab) {
   sh <- strheight("l", vfont=c("sans serif", "plain"))
   xpos <- apply(bp,2,mean)
     ypos.lab <- grconvertY(0.065, from="ndc")
     ypos.ttl <- grconvertY(1, from="nic")
   text(xpos, rep(ypos.lab, 6), c("even.F", "core.F", "edge.F"))
     if(counter ==1) {
       text(mean(xpos[1:3]), ypos.ttl, "Even dispersal",xpd=NA, cex=1.15)
       text(mean(xpos[4:6]), ypos.ttl, "Preferential dispersal",xpd=NA, cex=1.15)
       }

     pu <- par("usr")
     #text(mean(xpos[5:6]), 40+2*sh, "F: edges",xpd=NA, cex=1.5)
     vf <- c("sans serif","bold")
     text(pu[2],pu[4]-sh,hab.struct,  cex=1.25, pos=4, xpd=NA, vfont=vf)
     text(pu[2],pu[4]-3*sh,"habitat",  cex=1.25, pos=4, xpd=NA)
     text(pu[2],pu[4]-5*sh,"structure",  cex=1.25, pos=4, xpd=NA)

     if(counter ==2) {
       legend.ltb.2(pu[2],pu[3] + 5*sh,c("Core","Edge"),inset=c(0.2,0),
                 col=c("tomato","royalblue"),xpd=NA,pt.cex=2,pch=15,bty="n")
     }
 }
}

    counter <<- 0
   dmm <- lapply(scen.list, draw.panel)
    mtext("Relative % difference with no dispersal scenario, scaled against K",
          side=2,
          outer=TRUE, line=2.5,cex=1.05)
   fname <- sprintf("tm-range-dyn_histo-comp-scen_F=%s_%s.pdf", Fval*10, wmetric)
    pu <- par("usr")

   dev.copy2pdf(file=fname)

}

###############################################################
###############################################################
# same as bp.scenario but histogram format -- original function
histo.scenario <- function(wmetric="Bdecl", scen.list=list(ee.lowF, ce.lowF)) {

  if(!(wmetric %in% c("Bdecl", "SrcSink"))) stop()

    check.dev.size(7, 6)
    par(family="HersheySans", mfrow=c(2,1),
        omi=c(0.1,0.75,0.4,1.2), mai=c(0.1,0.1,0.1,0.2))
    Fval <- attr(scen.list[[1]], "Fval")

   draw.panel <- function(wscen) {

     counter <<- counter + 1
     rgval <- attr(wscen, "rgval")
     hab.struct <- ifelse(rgval[1]!=rgval[2], "Core-edge", "Even")

     # select function to compute metric to show
     funk2use <- ifelse(wmetric=="SrcSink", "scen.sourcesink",
                                        "comp.scen.ec")# option: comp.scen.ec, scen.sourcesink

     scen.emig <- do.call(funk2use, list(wscen, scen="emig"))
     scen.wpref <- do.call(funk2use, list(wscen, scen="wpref"))
     comp.mat <- rbind(scen.emig$mean, scen.wpref$mean)
     min.mat <- rbind(scen.emig$min, scen.wpref$min)
     max.mat <- rbind(scen.emig$max, scen.wpref$max)

     # reorder to have core edge cells with scenarios in even/core/edge F order
     comp.mat <- rbind(comp.mat[c(1,3,5,2,4,6)],
                       comp.mat[c(7,9,11,8,10,12)])
     min.mat <- rbind(min.mat[c(1,3,5,2,4,6)],
                       min.mat[c(7,9,11,8,10,12)])
     max.mat <- rbind(max.mat[c(1,3,5,2,4,6)],
                       max.mat[c(7,9,11,8,10,12)])

    spacev <- c(0.1, rep(c(0.1,1),2),0.1,3,rep(c(0.1,1),2),0.1)

     yl <- c(-0.35,0.35)
     if(funk2use == "scen.sourcesink") yl <- c(-1,1)

     bp <- barplot(comp.mat,ylim=yl, beside=TRUE,
                  las=1, col=NA,
                  space=spacev, xlim=c(0,20), border=NA)
     abline(h=0, col="grey")
     abline(h=c(-0.2,0.2),col="cornsilk3",lwd=0.5)
     bp <- barplot(comp.mat,ylim=yl, beside=TRUE, las=1, add=TRUE,
              col=c(rep(c("tomato1","royalblue1"),3),
              rep(c("tomato2","royalblue2"),3)),
              space=spacev, xlim=c(0,20), density=c(rep(50,6), rep(NA,6)),
              border=NA)
     box()

     segments(bp, c(min.mat), y1=c(max.mat), col="grey50")

     abline(v=mean(bp[2,])-0.5)


   ## Add labels
   sh <- strheight("l", vfont=c("sans serif", "plain"))
   xpos <- apply(bp,2,mean)
     ypos.lab <- grconvertY(0.065, from="ndc")
     ypos.ttl <- grconvertY(1, from="nic")
   text(xpos, rep(ypos.lab, 6), c("even.F", "core.F", "edge.F"))
     if(counter ==1) {
       text(mean(xpos[1:3]), ypos.ttl, "Even dispersal",xpd=NA, cex=1.15)
       text(mean(xpos[4:6]), ypos.ttl, "Preferential dispersal",xpd=NA, cex=1.15)
       }

     pu <- par("usr")
     #text(mean(xpos[5:6]), 40+2*sh, "F: edges",xpd=NA, cex=1.5)
     vf <- c("sans serif","bold")
     text(pu[2],pu[4]-sh,hab.struct,  cex=1.25, pos=4, xpd=NA, vfont=vf)
     text(pu[2],pu[4]-3*sh,"habitat",  cex=1.25, pos=4, xpd=NA)
     text(pu[2],pu[4]-5*sh,"structure",  cex=1.25, pos=4, xpd=NA)

     if(counter ==2) {
       legend.ltb.2(pu[2],pu[3] + 5*sh,c("Core","Edge"),inset=c(0.2,0),
                 col=c("tomato","royalblue"),xpd=NA,pt.cex=2,pch=15,bty="n")
     }
}

    counter <<- 0
   dmm <- lapply(scen.list, draw.panel)
    mtext("Relative % difference with no dispersal scenario, scaled against K",
          side=2,
          outer=TRUE, line=2.5,cex=1.05)
   fname <- sprintf("tm-range-dyn_histo-comp-scen_F=%s_%s.pdf", Fval*10, wmetric)
    pu <- par("usr")

   dev.copy2pdf(file=fname)

}
######################################################
######################################################
## Looks under what conditions core cells are source cells
## (and vice and versa)
## e.g. if we increase dispersal preference towards good environment
## we would expect core cells to receive more individuals,
## but this should be less strong if there is high fishing in the core
## effort

srcsinkX2 <- function(rg.core=0.1, rg.edge=0.05) {

  # habitat set-up defined in function arguments

  # define scenarios
  # scenario 1 is dispersal preference
  scen.prefdisp <- seq(0, 1, by=0.1)
  # scenario 2 is increasing fishing effort in the core
  # no fishing in the edges
  scen.fe <- seq(0, 1, by=0.1)

  # define sub-function that launches simulation and extracts
  # source-sink values (aggregated and by cell) for core/edge
  get.ss <- function(scen1.val, scen2.val) {

    # launch simulation
    tm.spatial.dyn(r.growth.core=rg.core, r.growth.edge=rg.edge,
                   fish.fact=scen2.val, fish.fact.edge=0,
                   pref.disp=scen1.val, add.r.pref=TRUE)

    # get ratio of emigration to immigration by cell
    em.mat <- envpop$emig.store[ts.max,]
    im.mat <- envpop$immig.store[ts.max,]
    ei.by.cells <- em.mat/im.mat-1
    ei.core <- sum(em.mat[core.layout$core.cells])/sum(im.mat[core.layout$core.cells])-1
    ei.edge <- sum(em.mat[core.layout$edge.cells])/sum(im.mat[core.layout$edge.cells])-1

    # get ratio of e-to-i by region (i.e. # of edge individuals received by core cells)
    pos <- calc.core.size(envpop$run.info$grid.width^2)
    xmat <- em.mat*envpop$immig.rate[,,ts.max] # matrix with individual exchange at equil
    # (senders in rows, receivers in columns)
    # take ratio of individuals that go from core to edge cells over
    # individuals that go from edge to core cells
    c2e <- sum(xmat[pos$core.cells,pos$edge.cells])
    e2c <- sum(xmat[pos$edge.cells,pos$core.cells])
    ei.reg <- c2e/e2c-1 # center on zero


    list(core=ei.core, edge=ei.edge, core.region=ei.reg, by.cells=ei.by.cells)
  }

  scenrun <- lapply(scen.prefdisp, function(pf)
                    lapply(scen.fe, function(sf) get.ss(pf, sf)))
  attr(scenrun,"scen.prefdisp") <- scen.prefdisp
  attr(scenrun,"scen.fe") <- scen.fe

  return(scenrun)
}

srcsinkX2.plot <- function(scen.list=srcsink.list, wvar="core") {

  nc <- length(scen.list[[1]]) # get number of scenarios
  nr <- length(scen.list[[1]][[1]])
  prefdisp <- attr(scen.list[[1]], "scen.prefdisp")
  feval <- attr(scen.list[[1]], "scen.fe")

  breakv <- c(-2,seq(-0.955,0.955,by=0.01),2)
  neutral.col <- "ivory"
  col.sink <- colorRampPalette(c(rev(brewer.pal(9,"Blues")[-(1:3)]),
                                 neutral.col, brewer.pal(9,"Greens")[-(1:3)]))
  colv <- col.sink(length(breakv)-1)

  scen.labs <- c("None", "Medium", "High")
  names(scen.labs) <- names(scen.list)

  draw.panel <- function(scen.now) {

    wscen <- scen.list[[scen.now]]

  # pref disp values in columns, fishing effort in rows
  mat <- sapply(1:nc,function(x)
                sapply(1:nr,function(y) wscen[[x]][[y]][wvar]))

  ## Make plot
  image(prefdisp, feval, mat, breaks=breakv, col=colv, las=1,
        yaxt=ifelse(scen.now=="low","t","n"), ann=FALSE,
        asp=1)
    box()
    mtext(scen.labs[scen.now])
  #contour(mat, breaks=breakv, add=TRUE)
  }

  check.dev.size(9.41, 3.79)
  par(family="HersheySans", mfrow=c(1,3),
      mai=c(0.2,0.1,0.2,0.1), omi=c(0.6,0.75,0.6,1.5))
  dmm <- lapply(c("low","med","high"), draw.panel)
  mtext("Preferential dispersal",side=1,outer=TRUE,line=2.5)
  mtext("Fishing effort in the core",side=2,outer=TRUE,line=3.5)
  mtext("Net individual movement in core cells under core-edge habitat difference",
        side=3,line=1.5, outer=TRUE, cex=1.25)

  # add legend
  pu <- par("usr")
  xpos <- grconvertX(0.85,from="ndc")

  legv <- rev(1+round(seq(-1, 1,by=0.2),1))
  leg.col <- rev(col.sink(length(legv)))
  lg <- legend.ltb.2(xpos,pu[4],xpd=NA,legend=legv,
                     col=leg.col,pch=15,bty="n",pt.cex=3.2,cex=1.2)

  tpos.x <- lg$text$x[1] + lg$rect$w
  text(tpos.x, lg$text$y[1], "Source", xpd=NA, cex=1.5)
  text(tpos.x, lg$text$y[length(lg$text$y)],
       "Sink", xpd=NA, cex=1.5)
  arrows(tpos.x, lg$text$y[2], y1=lg$text$y[length(lg$text$y)-1],
         length=0.1,xpd=NA, code=3)

  dev.copy2pdf(file="tm-range-dyn_srcsink-surface.pdf")
}

if(run.df.scen) {
if(!exists("srcsink.list")) {
  srcsink.scen.low <- srcsinkX2(0.1, 0.1)
  srcsink.scen.med <- srcsinkX2(0.1, 0.05)
  srcsink.scen.high <- srcsinkX2(0.1, 0.01)
  srcsink.list <- list(low=srcsink.scen.low, med=srcsink.scen.med,
                     high=srcsink.scen.high)
}
}

######################################################
######################################################


if(aor.scenario.rerun) {

  require(parallel)
  start.timer()
  rge.vect <- c(0.5, 0.3, 0.2, 0.1)
  # fishing mortality used in scenarios
  ff.all <- 0.5
  ff.edge <- ff.all
  aor.emig <- mclapply(rge.vect, function(rge) {
    tm.spatial.dyn(emig.base=0.1, r.growth.core=0.5, r.growth.edge=rge,
                   fish.fact=ff.all, fish.fact.edge=ff.edge, pref.disp=0, grid.width=10);
    envpop$mat})
  names(aor.emig) <- c("noss","lss","mss","hss")

  aor.wpref <- mclapply(rge.vect, function(rge) {
    tm.spatial.dyn(emig.base=0.1, r.growth.core=0.5, r.growth.edge=rge,
                   fish.fact=ff.all, fish.fact.edge=ff.edge, pref.disp=1, add.r.pref=TRUE, grid.width=10);
    envpop$mat})
  names(aor.wpref) <- c("noss","lss","mss","hss")


  aor.whighpref <- mclapply(rge.vect, function(rge) {
    tm.spatial.dyn(emig.base=0.1, r.growth.core=0.5, r.growth.edge=rge,
                   fish.fact=ff.all, fish.fact.edge=ff.edge, pref.disp=2, add.r.pref=TRUE, grid.width=10);
    envpop$mat})
  names(aor.whighpref) <- c("noss","lss","mss","hss")
  stop.timer()

 }


aor.scenarios <- function(pthresh=0.05) {

    emig.col <- col2transp("dark grey",0.8)
    wpref.col <- col2transp("royalblue3",0.8)
    wprefhigh.col <- col2transp("tomato",0.8)
    check.dev.size(8,8)
    par(mfrow=c(2,2), mai=rep(0.25, 4), omi=c(0.5, 0.5, 0.2, 0.1), family="HersheySans")

    plot.tm.aor(popmat=aor.emig$noss, pres.thresh=pthresh, col=emig.col)
    plot.tm.aor(popmat=aor.wpref$noss, pres.thresh=pthresh, add=TRUE, col=wpref.col)
    plot.tm.aor(popmat=aor.whighpref$noss, pres.thresh=pthresh, add=TRUE, col=wprefhigh.col)
    mtext("No core-edge spatial structure")

    plot.tm.aor(popmat=aor.emig$lss, pres.thresh=pthresh, col=emig.col)
    plot.tm.aor(popmat=aor.wpref$lss, pres.thresh=pthresh, add=TRUE, col=wpref.col)
    plot.tm.aor(popmat=aor.whighpref$lss, pres.thresh=pthresh, add=TRUE, col=wprefhigh.col)
    mtext("Low core-edge spatial structure")

    plot.tm.aor(popmat=aor.emig$mss, pres.thresh=pthresh, col=emig.col)
    plot.tm.aor(popmat=aor.wpref$mss, pres.thresh=pthresh, add=TRUE, col=wpref.col)
    plot.tm.aor(popmat=aor.whighpref$mss, pres.thresh=pthresh, add=TRUE, col=wprefhigh.col)
    mtext("Medium core-edge spatial structure")

    plot.tm.aor(popmat=aor.emig$hss, pres.thresh=pthresh, col=emig.col)
    plot.tm.aor(popmat=aor.wpref$hss, pres.thresh=pthresh, add=TRUE, col=wpref.col)
    plot.tm.aor(popmat=aor.whighpref$hss, pres.thresh=pthresh, add=TRUE, col=wprefhigh.col)
    mtext("High core-edge spatial structure")

    mtext("Range size (# of occupied squares)", outer=TRUE, side=2, line=2)
    mtext("Total population abundance", outer=TRUE, side=1, line=1)

    dev.copy2pdf(file="tm-range-contrxn_aor-scenario_emig-wpref-comp.pdf")
}

#######################################################
# Calculate net exchange of individuals between regions
core.edge.exch <- function(run.obj=envpop) {

  # get matrix positions for core/edge cells
  pos <- calc.core.size(run.obj$run.info$grid.width^2)

  byts <- function(ts) {

  # exchange matrix: senders in row, receivers in columns
  xmat <- run.obj$emig.store[ts,] * run.obj$immig.rate[,,ts]

  # to know how many individuals went from core cells to
  # edge cells, select core *rows* and edge *columns*
  # and take the sum
  c2e <- sum(xmat[pos$core.cells,pos$edge.cells])
  c2c <- sum(xmat[pos$core.cells,pos$core.cells]) # core to core
  e2c <- sum(xmat[pos$edge.cells,pos$core.cells]) # edge to core
  e2e <- sum(xmat[pos$edge.cells,pos$edge.cells]) # edge to edge

  # exchange ratio between regions is the number of core individuals
  # sent to edge cells over the # of indivs received from edge cells
  core.ei <- c2e/e2c # take the inverse for edge
}

  sapply(c(10,100,500,1000), byts)



}
