## Theortcl-model_Range-contrxn_BiomassIndic.r
## Comparing biomass in individual cells compared
## to overall biomass under different scenarios
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 18, 2014
## Time-stamp: <2014-07-19 12:42:45 Laura>

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

    rg.core.vect <- seq(0.05, 0.25, by=0.05)
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


