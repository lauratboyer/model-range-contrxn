## Theortcl-model_Range-contrxn_Scenarios.r
## Setting up scenarios to compare model runs
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 14, 2014
## Time-stamp: <2014-07-15 12:50:07 Laura>

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
    attr(scenl, "K") <- K
    attr(scenl, "rgval") <- c(core=rgcore, edge=rgedge)
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

    message("need to loop through all Ftypes")
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

}

run.histo <- FALSE
if(run.histo) {
core.val <- comp.scen.ec("core")
edge.val <- comp.scen.ec("edge")
comp.mat <- rbind(core.val, edge.val)*100

check.dev.size(6, 8)
#par(mfrow=c(2,1), family="HersheySans")
bp <- barplot(comp.mat,ylim=c(-50, 50), beside=TRUE, las=1, col=c("tomato", "turquoise3"))
abline(h=0, col="cornsilk3")
bp <- barplot(comp.mat,ylim=c(-50, 50), beside=TRUE, las=1,
              col=c("tomato", "turquoise3"), add=TRUE)
sh <- strheight("l", vfont=c("sans serif", "plain"))
xpos <- apply(bp,2,mean)
text(xpos, rep(40, 6), c("even", "pref"))
text(mean(xpos[1:2]), 40+2*sh, "even F",xpd=NA, cex=1.5)
text(mean(xpos[3:4]), 40+2*sh, "F: core",xpd=NA, cex=1.5)
text(mean(xpos[5:6]), 40+2*sh, "F: edges",xpd=NA, cex=1.5)
}
