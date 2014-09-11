## Theortcl-model_Range-contrxn_BiomassIndic.r
## Comparing biomass in individual cells compared
## to overall biomass under different scenarios
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 18, 2014
## Time-stamp: <2014-07-31 10:08:10 Laura>

# core.edge.fact controls the different in growth between
# core and edge regions
bioF <- function(core.edge.fact=1, fmort=0.5,
                 emig.base=0.1, pref.disp=0,
                 fish.fact.edge=fmort/core.edge.fact,
                 add.r.pref=FALSE) {

  rg.core <- 0.1
  rg.edge <- rg.core/core.edge.fact

  # no fishing
  tm.spatial.dyn(r.growth.core=rg.core,
                 r.growth.edge=rg.edge,
                 emig.base=emig.base,
                 fish.fact=0, pref.disp=pref.disp,
                 add.r.pref=add.r.pref)
  pop.nofish <- envpop$mat[,,ts.max]

  # with fishing
  tm.spatial.dyn(r.growth.core=rg.core,
                 r.growth.edge=rg.edge,
                 emig.base=emig.base,
                 fish.fact=fmort, fish.fact.edge=fish.fact.edge,
                 pref.disp=pref.disp,
                 add.r.pref=add.r.pref)
  pop.wfish <- envpop$mat[,,ts.max]

  # get population B_cur/B_0
  propBcur <- sum(pop.wfish)/sum(pop.nofish)

  # by cell B_cur/B_0
  bycell.bb <- pop.wfish/pop.nofish

  list(tot=propBcur, spatial=bycell.bb)
}

################################################
################################################
## Compute population abundance at equilibrium
## under various core/edge growth structures
abund.vs.ce <- function(emig.base=0.1, pref.disp=0, add.r.pref=FALSE) {

    rg.core.vect <- seq(0.05, 0.25, by=0.1)
    rg.edge.factor <- seq(1,5,by=0.5)

  get.bio <- function(rg.core=0.1, edge.fact=0.1) {

    tm.spatial.dyn(r.growth.core=rg.core,
                   r.growth.edge=rg.core/edge.fact,
                   emig.base=emig.base,
                   pref.disp=pref.disp,
                   add.r.pref=add.r.pref)
    sum(envpop$mat[,,ts.max])
  }

  sum.all <- sapply(rg.core.vect, function(rgc)
                    sapply(rg.edge.factor, function(rec)
                           get.bio(rgc, rec)))
}

### Run biomass sims
if(!exists("run.bio.sims")) {

  bio.prm.df <- expand.grid(emig.base=c(0.05,0.1,0.15,0.2,0.25),
                            pref.disp=c(0,1))
  run.bio.sims <- lapply(1:nrow(bio.prm.df),
                         function(x) abund.vs.ce(emig.base=bio.prm.df[x,1],
                                                 pref.disp=bio.prm.df[x,2]))
}
if(!exists("run.bio.sims.wpref")) {

  bio.prm.df <- expand.grid(emig.base=c(0.05,0.1,0.15,0.2,0.25),
                            pref.disp=c(0,0.1))
  run.bio.sims.wpref <- lapply(1:nrow(bio.prm.df),
                         function(x) abund.vs.ce(emig.base=bio.prm.df[x,1],
                                                 pref.disp=bio.prm.df[x,2],
                                                 add.r.pref=TRUE))
}

### Compare core-edge biomass
comp.ce.biom <- function(bio.list=run.bio.sims) {

  t1 <- bio.list[[1]]
  par(mai=c(0.8,0.8,0.5,1.5), family="HersheySans")
  pu <- par("usr")
  plot(1:nrow(t1), type="n", las=1, ylim=c(0.5,1),
       xlab="Habitat quality difference between the core and the edge",
       ylab="Proportion of total population at equilibrium compared to no-dispersal scenario")
  abline(h=1)
  colv <- tim.colors(length(bio.list))

  add.lines <- function(tnum) {

    dnow <- bio.list[[tnum]]
    colpal <- colorRampPalette(c("grey",colv[tnum]))
    colnow <- (colpal(4))
    dmm <- sapply(2,
                  function(x) lines(dnow[,x]/sum(K),
                                    col=colnow[x+1]))
  }

  dmm <- lapply(1:10, add.lines)


  legend.ltb.2(pu[2],pu[4],legend=apply(bio.prm.df,1,paste, collapse=" -- "),
         xpd=NA, lty=1, col=colv, bty="n")

}

############################################################
## Compare regional vs local remaining biomass
## Scenarios: no environmental structure, even, core, edge fishing

reg.vs.local.bio <- function(pref.disp=0.1,
                             core.fact=1, edge.fact=1) {

  ffv <- seq(0.1,1,by=0.1)
  rg.core <- 0.1
  rg.edge <- 0.05



  # first run model with this configuration and no fishing:
  tm.spatial.dyn(emig.base=rg.core, r.growth.edge=rg.edge,
                 pref.disp=pref.disp, add.r.pref=TRUE)
  baseline <- envpop$mat

  fish.sim <- function(fish.fact) {

    # now add fishing:
    tm.spatial.dyn(emig.base=rg.core, r.growth.edge=rg.edge, pref.disp=pref.disp, add.r.pref=TRUE,
                   fish.fact=fish.fact*core.fact, fish.fact.edge=fish.fact*edge.fact)
    wfish <- envpop$mat
    reg.bprop <- sum(wfish[,,ts.max])/sum(baseline[,,ts.max])
    loc.bprop <- wfish[,,ts.max]/baseline[,,ts.max]
    loc.df <- data.frame(cell.type=rep(c("core","edge"),sapply(core.layout[2:3],length)),
                         val=loc.bprop[unlist(core.layout[2:3])])
    loc.df$abs.diff <- (reg.bprop - loc.df$val)/reg.bprop*100
    ldspl <- split(loc.df, loc.df$cell.type)
    bstatfunk <- function(x) c(median=median(x),mean=mean(x),sd=sd(x),min=min(x),max=max(x))
    sapply(ldspl, function(xx) bstatfunk(xx$abs.diff))
  }

  simlist <- lapply(ffv, fish.sim); names(simlist) <- ffv
  simlist
}

if(!exists("reg.loc.bleft.allscen")) {
  reg.loc.bleft.even <- list(reg.vs.local.bio(0),
                             reg.vs.local.bio(0.5), reg.vs.local.bio(1))
  reg.loc.bleft.coreF <- list(reg.vs.local.bio(0, edge.fact=0),
                             reg.vs.local.bio(0.5, edge.fact=0),
                             reg.vs.local.bio(1, edge.fact=0))
  reg.loc.bleft.edgeF <- list(reg.vs.local.bio(0, core.fact=0),
                             reg.vs.local.bio(0.5, core.fact=0),
                             reg.vs.local.bio(1, core.fact=0))

reg.loc.bleft.allscen <- list(reg.loc.bleft.even,
                              reg.loc.bleft.coreF,
                              reg.loc.bleft.edgeF)
}
### Plot difference between regional b_msy vs local b_msy for core and edge cells
### under different fishing/habitat/dispersal scenarios
plot.reg.vs.local.bio <- function(wl.all=reg.loc.bleft.allscen) {

  yl <- c(-60,60)
  num.f <- length(wl.all[[1]][[1]])
  val.f <- as.numeric(names(wl.all[[1]][[1]]))
  colreg <- col2transp(c("tomato","dodgerblue"))
  colmed <- c("tomato3","dodgerblue3")
  names(colreg) <- c("core","edge")
  names(colmed) <- c("core","edge")

  check.dev.size(9.95, 7.75)
  par(family="HersheySans", mfrow=c(3,3), mai=c(0.1,0.1,0.1,0.1),
      omi=c(0.6,0.6,0.6,0.4))
  ylabel <- "Relative difference, regional to local Bcur/B0: (reg-loc)/reg"

  add.panel <- function(wl) {

    counter <<- counter + 1
    plot(1:10,type="n",xlim=range(val.f),ylim=yl,las=1,
         ylab="",xlab="",
         xaxt=ifelse(counter>=7,"t","n"),
         yaxt=ifelse(counter%in%c(1,4,7) ,"t","n"))
    abline(h=seq(-40,40,by=20),
         col=c("cornsilk2","cornsilk3","cornsilk4",
           "cornsilk3","cornsilk2"))
    pu <- par("usr")
    if(counter == 1) mtext("No preferential dispersal")
    if(counter == 2) mtext("Low preferential dispersal")
    if(counter == 3) mtext("High preferential dispersal")
    if(counter == 1) mtext("Even F",side=2,line=3)
    if(counter == 4) mtext("F in core only",side=2,line=3)
    if(counter == 7) mtext("F in edges only",side=2,line=3)


  add.band <- function(wreg) {

    core.poly.y <- c(sapply(wl, function(x) x["max",wreg]),
                   rev(sapply(wl, function(x) x["min",wreg])))
    core.poly.x <- c(val.f,rev(val.f))
    polygon(core.poly.x, core.poly.y, col=colreg[wreg], density=NA, border=NA)
    lines(val.f, sapply(wl, function(x) x["median",wreg]),lwd=2,col=colmed[wreg])
  }

  add.band("core")
  add.band("edge")
  }

  counter <<- 0
  dmm <- lapply(wl.all[[1]], add.panel)
  dmm <- lapply(wl.all[[2]], add.panel)
  dmm <- lapply(wl.all[[3]], add.panel)
  mtext(ylabel, side=3, outer=TRUE, line=2.5, adj=0, cex=1.2)
  mtext("Increasing fishing mortality as factor of r_g in core",
        side=1, outer=TRUE, line=2.5)

  dev.copy2pdf(file="Theo-mod_range-contrxn_reg-Bcur-vs-loc-Bcur.pdf")

}

###############################################################
###############################################################
## relative abundance over time CPUE
cpue.ts <- function(emat=envpop) {

  # cpue at the population level:
  fs <- emat$run.info$fishing.start
  catch.ts <- apply(emat$catch.mat[,,fs:ts.max], 3, sum, na.rm=TRUE)
  effort.ts <- apply(emat$F.array[,,fs:ts.max], 3, sum, na.rm=TRUE)   # assuming q = 1
  effort.ts[effort.ts == 0] <- NA
  catch.ts[effort.ts == 0] <- NA
  cpue.ts <- catch.ts/effort.ts
  stand.cpue.ts <- cpue.ts/cpue.ts[1]

  # abundance at the population level
  abund.ts <- apply(emat$mat[,,fs:ts.max], 3, sum, na.rm=TRUE)
  abund.ts <- abund.ts/abund.ts[1]

  # by cell
  cpue.cell <- function(wcell) {

    cell.pos <- arrayInd(wcell, .dim=dim(emat$mat[,,1]))
    catch.ts <- emat$catch.mat[cell.pos[1],cell.pos[2],fs:ts.max]
    effort.ts <- emat$F.array[cell.pos[1],cell.pos[2],fs:ts.max] # using F=qE and q=1
    effort.ts[effort.ts == 0] <- NA
    catch.ts[effort.ts == 0] <- NA
    cpue.ts <- catch.ts/effort.ts
    stand.cpue.ts <- cpue.ts/cpue.ts[1]
    abund.ts <- emat$mat[cell.pos[1],cell.pos[2],fs:ts.max]
    abund.ts <- abund.ts/abund.ts[1]

#    lines(stand.cpue.ts, col=col.mat.transp[wcell])
    #lines(abund.ts/stand.cpue.ts, col=col.mat.transp[wcell])
    lines(abund.ts, col=col.mat.transp[wcell])
  }

  ncell <- emat$run.info$grid.width^2

  plot(stand.cpue.ts, type="n", ylim=c(0,1.1))
  dmm <- sapply(1:ncell, cpue.cell)
  lines(stand.cpue.ts, lwd=2)
  lines(abund.ts, col="green")

}
