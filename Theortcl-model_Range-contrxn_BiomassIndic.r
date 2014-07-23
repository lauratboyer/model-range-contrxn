## Theortcl-model_Range-contrxn_BiomassIndic.r
## Comparing biomass in individual cells compared
## to overall biomass under different scenarios
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 18, 2014
## Time-stamp: <2014-07-23 10:44:10 Laura>

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
                            pref.disp=c(0,1))
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
