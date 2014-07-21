## Theortcl-model_Range-contrxn_Scenarios.r
## Setting up scenarios to compare model runs
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 14, 2014
## Time-stamp: <2014-07-20 16:55:14 Laura>

########################################################
scenario.run <- function(rgcore = 0.2, rgedge=rgcore, Fval=0.5) {

## Constant environment
## fishing everywhere
    evenfish <- list()
    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, emig.base=0)
    evenfish$base <- envpop$mat

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, emig.base=0.1)
    evenfish$emig <- envpop$mat

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, emig.base=0.1, pref.disp=1, add.r.pref=TRUE)
    evenfish$wpref <- envpop$mat

## ... with fishing only in core
    corefish <- list()
    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, fish.fact.edge=0, emig.base=0)
    corefish$base <- envpop$mat

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, fish.fact.edge=0, emig.base=0.1)
    corefish$emig <- envpop$mat

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=Fval, fish.fact.edge=0, emig.base=0.1,
               pref.disp=1, add.r.pref=TRUE)
    corefish$wpref <- envpop$mat

    ## ... with fishing only in edge
    edgefish <- list()
    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=0, fish.fact.edge=Fval, emig.base=0)
    edgefish$base <- envpop$mat

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=0, fish.fact.edge=Fval, emig.base=0.1)
    edgefish$emig <- envpop$mat

    tm.spatial.dyn(r.growth.core=rgcore, r.growth.edge=rgedge,
               fish.fact=0, fish.fact.edge=Fval, emig.base=0.1,
               pref.disp=1, add.r.pref=TRUE)
    edgefish$wpref <- envpop$mat

    scenl <- list(even=evenfish, core=corefish, edge=edgefish)
    attr(scenl, "ts.max") <- ts.max
    attr(scenl, "Fval") <- Fval
    attr(scenl, "K") <- K
    attr(scenl, "rgval") <- c(core=rgcore, edge=rgedge)
    attr(scenl, "grid.width") <- grid.width
    attr(scenl, "core.layout") <- core.layout
    return(scenl)
}

## For each scenario calculate the different in biomass
## expected compared to no dispersal version

if(!exists("ee.lowF")) ee.lowF <- scenario.run()
if(!exists("ee.highF")) ee.highF <- scenario.run(Fval=0.8)
if(!exists("ce.lowF")) ce.lowF <- scenario.run(rgedge=0.1)
if(!exists("ce.highF")) ce.highF <- scenario.run(rgedge=0.1, Fval=0.8)

comp.scen.ec <- function(sim.list=ee.lowF, scen="emig",
                         calc.type="absolute") { #scen is emig or wpref

    tsm <- attr(sim.list, "ts.max")
    var.in <- expand.grid(F.type=c("even","core","edge"),
                          region.focus=c("core","edge"),
                          stringsAsFactors=FALSE)

    getstat <- function(ri) { # wregion is core or edge

        F.type <- var.in[ri,1]
        wregion <- var.in[ri,2]
        sl.scen <- sim.list[[F.type]][[scen]][,,tsm]
        sl.base <- sim.list[[F.type]]$base[,,tsm]

        varbl <- paste(wregion, ".cells", sep="")
        if(calc.type=="absolute") {
            valnow <- ((sl.scen-sl.base)/K)[core.layout[[varbl]]]

        }else{
            valnow <- ((sl.scen/sl.base))[core.layout[[varbl]]]-1
        }

        c(scen.type=paste(paste(F.type,"F",sep=""), varbl, sep="."),
          100*c(mean=mean(valnow), sd=sd(valnow)))
    }

    out <- as.data.frame(t(sapply(1:nrow(var.in), getstat)))
    class(out$mean) <- "numeric"
    class(out$sd) <- "numeric"
    out

}

histo.scenario <- function(scen.list) {

    check.dev.size(7, 6)
    par(family="HersheySans", mfrow=c(2,1),
        omi=c(0.1,0.6,0,0.2), mai=c(0.1,0.25,0.95,0.2))
    Fval <- attr(scen.list[[1]], "Fval")

   draw.panel <- function(wscen) {

     rgval <- attr(wscen, "rgval")
     hab.struct <- ifelse(rgval[1]!=rgval[2], "Core-edge", "Even")

     scen.emig <- comp.scen.ec(wscen, scen="emig")
    scen.wpref <- comp.scen.ec(wscen, scen="wpref")
    comp.mat <- rbind(scen.emig$mean, scen.wpref$mean)
     comp.mat <- rbind(comp.mat[c(1,3,5,2,4,6)],
                       comp.mat[c(7,9,11,8,10,12)])

    spacev <- c(0.1, rep(c(0.1,1),2),0.1,3,rep(c(0.1,1),2),0.1)

    bp <- barplot(comp.mat,ylim=c(-50, 50), beside=TRUE,
                  las=1, col=NA,
                  space=spacev, xlim=c(0,20), border=NA)
   abline(h=0, col="grey")
   bp <- barplot(comp.mat,ylim=c(-50, 50), beside=TRUE, las=1, add=TRUE,,
              col=c(rep(c("tomato1","royalblue1"),3),
              rep(c("tomato2","royalblue2"),3)),
              space=spacev, xlim=c(0,20), density=c(rep(50,6), rep(NA,6)),
              border=NA)

   ## Add labels
   sh <- strheight("l", vfont=c("sans serif", "plain"))
   xpos <- apply(bp,2,mean)
   text(xpos, rep(40, 6), c("even.F", "core.F", "edge.F"))
   text(mean(xpos[1:3]), 40+2.5*sh, "EVEN EMIG",xpd=NA, cex=1.25)
   text(mean(xpos[4:6]), 40+2.5*sh, "PREF EMIG",xpd=NA, cex=1.25)
   #text(mean(xpos[5:6]), 40+2*sh, "F: edges",xpd=NA, cex=1.5)
   mtext(sprintf("Habitat structure: %s", hab.struct),
         adj=0, cex=1.5, line=2)
}

   dmm <- lapply(scen.list, draw.panel)
    mtext("Relative % difference with no dispersal scenario, scaled against K", side=2,
          outer=TRUE, line=2)
   fname <- sprintf("tm-range-dyn_histo-comp-scen_F%s.pdf", Fval*10)
   dev.copy2pdf(file=fname)

}


aor.scenario.rerun <- FALSE
if(aor.scenario.rerun) {

    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.5,
                   fish.fact=0.5, pref.disp=0, grid.width=10)
    epop.aor.emig.noss <- envpop$mat
    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.5,
                   fish.fact=0.5, pref.disp=1, add.r.pref=TRUE, grid.width=10)
    epop.aor.wpref.noss <- envpop$mat

    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.3,
                   fish.fact=0.5, pref.disp=0, grid.width=10)
    epop.aor.emig.lss <- envpop$mat
    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.3,
                   fish.fact=0.5, pref.disp=1, add.r.pref=TRUE, grid.width=10)
    epop.aor.wpref.lss <- envpop$mat

    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.2,
                   fish.fact=0.5, pref.disp=0, grid.width=10)
    epop.aor.emig.mss <- envpop$mat
    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.2,
                   fish.fact=0.5, pref.disp=1, add.r.pref=TRUE, grid.width=10)
    epop.aor.wpref.mss <- envpop$mat

    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.1,
                   fish.fact=0.5, pref.disp=0, grid.width=10)
    epop.aor.emig.vhss <- envpop$mat
    tm.spatial.dyn(emig.base=0.2, r.growth.core=0.5, r.growth.edge=0.1,
                   fish.fact=0.5, pref.disp=1, add.r.pref=TRUE, grid.width=10)
    epop.aor.wpref.vhss <- envpop$mat
}


aor.scenarios <- function() {

    emig.col <- ("dark grey")
    wpref.col <- ("royalblue3")
    par(mfrow=c(2,2), mai=rep(0.25, 4), omi=c(0.5, 0.5, 0.2, 0.1), family="HersheySans")

    plot.tm.aor(popmat=epop.aor.emig.noss, pres.thresh=0.05, col=emig.col)
    plot.tm.aor(popmat=epop.aor.wpref.noss, pres.thresh=0.05, add=TRUE, col=wpref.col)
    mtext("No core-edge spatial structure")

    plot.tm.aor(popmat=epop.aor.emig.lss, pres.thresh=0.05, col=emig.col)
    plot.tm.aor(popmat=epop.aor.wpref.lss, pres.thresh=0.05, add=TRUE, col=wpref.col)
    mtext("Low core-edge spatial structure")

    plot.tm.aor(popmat=epop.aor.emig.mss, pres.thresh=0.05, col=emig.col)
    plot.tm.aor(popmat=epop.aor.wpref.mss, pres.thresh=0.05, add=TRUE, col=wpref.col)
    mtext("Medium core-edge spatial structure")

    plot.tm.aor(popmat=epop.aor.emig.vhss, pres.thresh=0.05, col=emig.col)
    plot.tm.aor(popmat=epop.aor.wpref.vhss, pres.thresh=0.05, add=TRUE, col=wpref.col)
    mtext("High core-edge spatial structure")

    mtext("Range size (# of occupied squares)", outer=TRUE, side=2, line=2)
    mtext("Total population abundance", outer=TRUE, side=1, line=1)

    dev.copy2pdf(file="tm-range-contrxn_aor-scenario_emig-wpref-comp.pdf")
}
