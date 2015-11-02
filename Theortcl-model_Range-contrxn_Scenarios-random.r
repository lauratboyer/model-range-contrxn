## Theortcl-model_Range-contrxn_Scenarios-random.r
## Version of scenario code applied to random distribution
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (lauratb@spc.int)
## Written on: September 16, 2015
## Time-stamp: <2015-10-01 15:00:13 lauratb>

##
theme.lolo <- theme_bw() + theme(text=element_text(family="Segoe UI Light", size=12),
                                 axis.title.x=element_text(vjust=-0.615),
                                 axis.title.y=element_text(vjust=0.95),
                            panel.grid=element_blank(),
                            panel.margin=unit(0.4,"cm"), strip.background = element_rect(fill="slategrey",colour="slategrey"),
                            strip.text=element_text(colour="white", size=12.5, vjust=0.4))
theme_set(theme.lolo)
# using argument list args.ce.low.lowF and args.ce.high.lowF from ..._Scenarios.r
# these are the scenarios used in Fig 3 of original chapter version from Fall 2014
if(!exists("ce.hlF.even.prefdisp")){
start.timer()
ce.llF.even.nodisp <- spadyn.sims(args.ce.low.lowF$evenfish$nodisp)
stop.timer()
ce.llF.even.evdisp <- spadyn.sims(args.ce.low.lowF$evenfish$evdisp)
ce.llF.even.prefdisp <- spadyn.sims(args.ce.low.lowF$evenfish$prefdisp)

ce.hlF.even.nodisp <- spadyn.sims(args.ce.high.lowF$evenfish$nodisp)
ce.hlF.even.evdisp <- spadyn.sims(args.ce.high.lowF$evenfish$evdisp)
ce.hlF.even.prefdisp <- spadyn.sims(args.ce.high.lowF$evenfish$prefdisp)

##### Figure 4 barplot comparisons scenarios
ee.lowF.even.nodisp <- spadyn.sims(args.eelowF$evenfish$nodisp) # reference case
ee.lowF.core.nodisp <- spadyn.sims(args.eelowF$corefish$nodisp) # for even
ee.lowF.edge.nodisp <- spadyn.sims(args.eelowF$edgefish$nodisp) # habitat structure
ee.lowF.even.evdisp <- spadyn.sims(args.eelowF$evenfish$evdisp) # even disp
ee.lowF.core.evdisp <- spadyn.sims(args.eelowF$corefish$evdisp)
ee.lowF.edge.evdisp <- spadyn.sims(args.eelowF$edgefish$evdisp)
ee.lowF.even.prefdisp <- spadyn.sims(args.eelowF$evenfish$prefdisp) # preferential disp
ee.lowF.core.prefdisp <- spadyn.sims(args.eelowF$corefish$prefdisp)
ee.lowF.edge.prefdisp <- spadyn.sims(args.eelowF$edgefish$prefdisp)

######################################
ce.lowF.even.nodisp <- spadyn.sims(args.celowF$evenfish$nodisp) # reference case
ce.lowF.core.nodisp <- spadyn.sims(args.celowF$corefish$nodisp) # for core-edge
ce.lowF.edge.nodisp <- spadyn.sims(args.celowF$edgefish$nodisp) # habitat structure
ce.lowF.even.evdisp <- spadyn.sims(args.celowF$evenfish$evdisp) # even disp
ce.lowF.core.evdisp <- spadyn.sims(args.celowF$corefish$evdisp)
ce.lowF.edge.evdisp <- spadyn.sims(args.celowF$edgefish$evdisp)
ce.lowF.even.prefdisp <- spadyn.sims(args.celowF$evenfish$prefdisp) # preferential disp
ce.lowF.core.prefdisp <- spadyn.sims(args.celowF$corefish$prefdisp)
ce.lowF.edge.prefdisp <- spadyn.sims(args.celowF$edgefish$prefdisp)

########################################
## make ggplot object
llF.nodisp <- ts.avrg.simu(ce.llF.even.nodisp)+guides(colour=FALSE, fill=FALSE)
llF.evdisp <- ts.avrg.simu(ce.llF.even.evdisp)+guides(colour=FALSE, fill=FALSE)
llF.prefdisp <- ts.avrg.simu(ce.llF.even.prefdisp)+guides(colour=FALSE, fill=FALSE)

hlF.nodisp <- ts.avrg.simu(ce.hlF.even.nodisp)+guides(colour=FALSE, fill=FALSE)
hlF.evdisp <- ts.avrg.simu(ce.hlF.even.evdisp)+guides(colour=FALSE, fill=FALSE)
hlF.prefdisp <- ts.avrg.simu(ce.hlF.even.prefdisp)+guides(colour=FALSE, fill=FALSE)

grid.arrange(llF.nodisp, llF.evdisp, llF.prefdisp, hlF.nodisp, hlF.evdisp, hlF.prefdisp, nrow=2)
}

# calculates source-sink metric for one simulation
srcsink.rd <- function(obj.1sim) {

    ts.max <- obj.1sim$run.info$ts
    core.cells <- obj.1sim$core.layout$core.cells
    edge.cells <- obj.1sim$core.layout$edge.cells
    core.emig.n <- obj.1sim$emigmat[core.cells,ts.max]
    edge.emig.n <- obj.1sim$emigmat[edge.cells,ts.max]
    core.immig.n <- obj.1sim$immigmat[core.cells,ts.max]
    edge.immig.n <- obj.1sim$immigmat[edge.cells,ts.max]

    # source-sink ratio is emigrants/immigrants
    list(core.ss=mean(core.emig.n/core.immig.n-1), edge.ss=mean(edge.emig.n/edge.immig.n-1))
}

ss.rd.launch <- function(rg.core=0.1, rg.edge=0.05, nsim=10, ss.metric="ss.cells.core") {

    # define scenarios
  # scenario 1 is dispersal preference
    scen.prefdisp <- seq(0, 1, by=0.1)
  # scenario 2 is increasing fishing effort in the core
  # no fishing in the edges
  scen.fe <- seq(0, 1, by=0.1)

  go.sim <- function(scen.pf, scen.F) {
      simu.obj <- spadyn.sims(list(r.growth.core=rg.core, r.growth.edge=rg.edge,
                                   fish.fact=scen.F, fish.fact.edge=0,
                                   add.r=TRUE, pref.disp=scen.pf), nsim=nsim)
      a1 <- colMeans(sapply(simu.obj$all.sims, function(i) c(i$ss.cells.core, i$ss.regions.core2edge))) # take mean value of metric by cell
  }

    scenrun0 <- lapply(scen.prefdisp, function(pd) lapply(scen.fe, function(ff) go.sim(pd, ff)))
    metr1 <- sapply(1:length(scen.prefdisp), function(pd) sapply(1:length(scen.fe), function(ff) scenrun0[[pd]][[ff]][1]))
    metr2 <- sapply(1:length(scen.prefdisp), function(pd) sapply(1:length(scen.fe), function(ff) scenrun0[[pd]][[ff]][2]))
    scenrun <- data.frame(pref.disp=rep(scen.prefdisp, each=length(scen.fe)), fcore=scen.fe,
                          ss.cells.core=c(metr1), ss.regions.core2edge=c(metr2))

  fname <- sprintf("source-sink-simu_core %s_edge %s_%s.RData", rg.core, rg.edge, ss.metric)
  save(scenrun, file=fname)
  invisible(scenrun)

}

############################################################################
sourcesink.heatplotX3.random <- function(wmetr="ss.regions.core2edge") {

    aa <- load("source-sink-simu_core 0.1_edge 0.1_ss.cells.core.RData")
    scenrun$scen.name <- "None"
    scenrun.none <- scenrun

    aa <- load("source-sink-simu_core 0.1_edge 0.05_ss.cells.core.RData")
    print(head(get("scenrun",.GlobalEnv)))
    scenrun$scen.name <- "Medium"
    scenrun.med <- scenrun

    aa <- load("source-sink-simu_core 0.1_edge 0.01_ss.cells.core.RData")
    scenrun$scen.name <- "High"
    scenrun.high <- scenrun
    scen.df <<- rbind(scenrun.none, scenrun.med, scenrun.high)

    scen.df$var.expl <- scen.df[,wmetr]
    scen.df$scen.name %<>% factor(levels=c("None","Medium","High"), ordered=TRUE)
    breakv <- seq(-0.65,0.75,by=0.1)
    breakv.labs <- round(seq(-0.65,0.75,by=0.1) + 0.05, 2)
    colv <- c(colorRampPalette(c("royalblue1","grey85"))(8),"grey95",colorRampPalette(c("grey85","seagreen2"))(7))
    scen.df$cutv <- cut(scen.df$var.expl, breaks=breakv)
    intervls <- levels(scen.df$cutv)
    #    print(data.frame(intervls, breakv.labs))
    check.dev.size(11.5, 4.5)
    g1 <- ggplot(scen.df, aes(x=pref.disp, y=fcore, fill=cutv)) + facet_wrap(~scen.name)+
        geom_tile(colour="grey") +
            scale_fill_manual(values=colv[-(1:2)], limits=intervls, labels=breakv.labs, name="Source-sink") +
            coord_equal() + xlab("Preferential dispersal") + ylab("F in core (as factor of r)")

    ggsave("tm-range-dyn_srcsink-surface_RND-DISTRIB.pdf", g1)
}

launch.biom.indic.random <- function(nsim=10, pref.disp=0, core.fact=1, edge.fact=1) {

    fname <- sprintf("biom-indic-simu_prefdisp %s_corefact %s_edgefact %s.RData", pref.disp, core.fact, edge.fact)
    ffv <- seq(0.1,1,by=0.1)
    rg.core <- 0.1
    rg.edge <- 0.05

    fish.loop <- function(F) {
        simu.args <- list(r.growth.core=rg.core, r.growth.edge=rg.edge,
                          fish.fact=F*core.fact, fish.fact.edge=F*edge.fact,
                          pref.disp=pref.disp, add.r=TRUE)
        simu.obj <- spadyn.sims(simu.args, nsim=nsim)$all
        bcur.df <- lapply(simu.obj, function(i) with(i$cell.Bcur, tapply(Bcur.diff, cell.type, quantile, c(0.25,0.5,0.75))))
        core.25 <- mean(sapply(bcur.df, function(i) i$core[1])) # mean of 25th quantile for Bcur.diff by simulation
        core.50 <- mean(sapply(bcur.df, function(i) i$core[2]))
        core.75 <- mean(sapply(bcur.df, function(i) i$core[3]))
        edge.25 <- mean(sapply(bcur.df, function(i) i$edge[1])) # mean of 25th quantile for Bcur.diff by simulation
        edge.50 <- mean(sapply(bcur.df, function(i) i$edge[2]))
        edge.75 <- mean(sapply(bcur.df, function(i) i$edge[3]))
        c(core.25=core.25, core.50=core.50, core.75=core.75,
          edge.25=edge.25, edge.50=edge.50, edge.75=edge.75)
    }

    bi.sims <- lapply(ffv, fish.loop)
    bi.df <- data.frame(F=ffv, do.call(rbind, bi.sims), pref.disp=pref.disp,
                        core.F.fact=core.fact, edge.F.fact=edge.fact)
    bi.df$F.label <- paste(bi.df$core.F, bi.df$edge.F, sep="-")

    save(bi.df, file=fname)
    invisible(bi.df)
}

fig.biom.indic.random <- function() {

    pref.disp.v <- rep(c(0, 1, 0.5), 3)
    core.fact.v <- rep(c(1, 0, 1),each=3)
    edge.fact.v <- rep(c(1, 1, 0),each=3)


    obj0 <- lapply(1:9, function(i)
        {aa <- load(sprintf("biom-indic-simu_prefdisp %s_corefact %s_edgefact %s.RData",
                            pref.disp.v[i], core.fact.v[i], edge.fact.v[i]));
                 get(aa)})

    df <- do.call(rbind, obj0) # bind all scenarios together in single data frame
    df[,c("core.50","core.25","core.75","edge.50","edge.25","edge.75")] %<>% "*"(100)
    F.label.vect <- c("1-1"="Fishing in core and edge", "1-0"="Core fishing only", "0-1"="Edge fishing only")
    df$F.label2 <- factor(F.label.vect[df$F.label], levels=F.label.vect, ordered=TRUE)
    pref.disp.vect <- c("0"="No preferential dispersal",
                        "0.5"="Low preferential dispersal", "1"="High preferential dispersal")
    df$pd.lab <- factor(pref.disp.vect[as.character(df$pref.disp)], levels=pref.disp.vect, ordered=TRUE)

    # make ggplot
    ggplot(dat=df) +
        geom_hline(aes(yintercept=seq(-50,50,by=25)), colour="cornsilk4", size=0.2) +
        geom_hline(aes(yintercept=0), size=0.2) +
        geom_ribbon(aes(x=F, y=core.50, ymin=core.25, ymax=core.75), fill="tomato", alpha=0.6) +
            geom_ribbon(aes(x=F, y=edge.50, ymin=edge.25, ymax=edge.75), fill="royalblue", alpha=0.6) +
                geom_line(aes(x=F, y=core.50)) + geom_line(aes(x=F, y=edge.50)) +
                    facet_grid(F.label2 ~ pd.lab) +
                        coord_cartesian(ylim=c(-65,65)) +
                            xlab("Increasing fishing mortality as factor of r_g in core") +
                            ylab("Relative difference, regional to local Bcur")

    ggsave("Theo-mod_range-contrxn_reg-Bcur-vs-loc-Bcur_RND-DISTRIB.png")
}

launch.aor.random <- function(pref.disp=0, rg.core=0.5, ff.all=0.5, nsim=100, rg.cv=0.1) {

    rge.vect <- c(0.5, 0.3, 0.2, 0.1)/0.5*rg.core

    # fishing mortality used in scenarios

    ff.edge <- ff.all

    go.simu <- function(rge) {

        print(rge)
        simu.args <- list(emig.base=0.1, r.growth.core=rg.core, r.growth.edge=rge,
                          fish.fact=ff.all, fish.fact.edge=ff.edge,
                          pref.disp=pref.disp, add.r=TRUE, grid.width=10, rg.cv=rg.cv)
        simu.obj <- spadyn.sims(simu.args, nsim=nsim)$all
        abund.df0 <- lapply(1:length(simu.obj), function(i) data.frame(r.edge=rge,
                                                                       pref.disp=pref.disp,
                                                                       run.id=i,
                                                                       ts=1:length(simu.obj[[i]]$abund.ts),
                                                                       abund=simu.obj[[i]]$abund.ts,
                                                                       area=simu.obj[[i]]$area.ts))
        do.call(rbind, abund.df0)
        }

    rge.sims <- lapply(rge.vect, go.simu)
    names(rge.sims) <- c("noss","lss","mss","hss")

    fname <- sprintf("aor-simu_prefdisp %s_rgcore %s_rgcv %s_F %s.RData", pref.disp, rg.core, rg.cv, ff.all)
    print(fname)
    save(rge.sims, file=fname)
    invisible(rge.sims)
}

fig.aor.random <- function(pref.disp=0, fval=0.8) {

    format.df <- function(pd.val) {
    load(sprintf("aor-simu_prefdisp %s_rgcore 0.5_rgcv 0.3_F %s.RData", pd.val, fval)) # rge.sims
    dat.all <- do.call(rbind, rge.sims)
    dat.all$abund %<>% "/"(1000)
    re.vect <- unique(dat.all$r.edge)
    gdat <- filter(dat.all, ts>50)
    gdat$r.edge <- factor(gdat$r.edge, levels=re.vect)
    gdat$pref.disp <- factor(gdat$pref.disp, levels=c(0,0.5,1,2))
    gdat$run.id.r.edge <- paste(gdat$pref.disp, gdat$run.id, gdat$r.edge)
    gdat$pref.disp.r.edge <- paste(gdat$pref.disp, gdat$r.edge)
    # calculate B0
    B0df <- filter(gdat, ts == 499) %>% mutate(B0=abund) %>% select(run.id.r.edge, B0)
    gdat %<>% inner_join(B0df, by="run.id.r.edge") %>% mutate(Bcur=abund/B0)
    gdat$ts5 <- 1*floor(gdat$ts/1)
    gdat %<>% group_by(pref.disp.r.edge, pref.disp, r.edge, ts5) %>%
        summarize(Bcur=mean(Bcur, na.rm=TRUE),
                  area.25=quantile(area, 0.25, na.rm=TRUE), area.75=quantile(area, 0.75, na.rm=TRUE),
                  area=median(area, na.rm=TRUE))
    gdat
}
    start.timer()
    a0 <- lapply(c(0, 0.5, 1, 2), format.df)
    gdat <- do.call(rbind, a0)
    gdat <<- gdat
    stop.timer()

    bpal <- c("grey",col2transp(brewer.pal(3, "YlGnBu")))
    g1 <- ggplot(dat=gdat,
           #           aes(x=abund, y=area, color=re, fill=re, ymin=area.25, ymax=area.75)) +
           aes(x=Bcur, y=area, group=pref.disp.r.edge, color=r.edge, shape=pref.disp)) +
        scale_color_manual(values=bpal, name="Spatial structure", label=c("None","Low","Med","High")) +
            scale_shape_discrete(name="Preferential dispersal", label=c("None","Low","Medium","High"),
                                 guide=guide_legend(override.aes=list(size=4))) +
#        geom_ribbon(color="grey", alpha=0.5) +
            geom_line(alpha=0.5, size=1.2) +
                geom_point(alpha=0.75, size=3) +
                    coord_cartesian(xlim=c(0,1)) + xlab("Bcur/B0") + ylab("Area (# squares)") +
    theme(legend.key=element_blank())
    check.dev.size(9,7)
    ggsave(sprintf("tm-range-contrxn_aor-scenario_aor-only_F %s_RND-DISTRIB.pdf", fval), g1)
}

calc.aor.metric.random <- function(fval=0.8) {

    aa <- load(sprintf("aor-simu_prefdisp 0_rgcore 0.5_rgcv 0.3_F %s.RData", fval))
    o1 <- na.omit(do.call(rbind, rge.sims))
#    abund.ts <- apply(popmat, 3, sum)
#    area.ts <- apply((floor(popmat))>(pres.thresh*K), 3, sum)

    ## get 1:1 line
    get.metr <- function(area, abund) {
    sl <- max(area)/max(abund)
    diff <- sum(area - sl*abund)
    sign <- ifelse(diff < 0, -1, 1)

    s1 <- sign*1/log(diff) # metric -- retain sign but log to diminish differences
}

    o1 %>% group_by(r.edge) %>% summarize(aor.val=get.metr(area, abund))
}

####################################################################
####################################################################
####################################################################
rerun.ss <- FALSE
if(rerun.ss) {

start.timer()
ss.rd.launch(0.1, 0.1)
stop.timer()
start.timer()
ss.rd.launch(0.1, 0.05)
stop.timer()
start.timer()
ss.rd.launch(0.1, 0.01)
stop.timer()
}

rerun.biomindic <- FALSE
if(rerun.biomindic){

    launch.biom.indic.random(nsim=100, pref.disp=0)
    launch.biom.indic.random(nsim=100, pref.disp=0.5)
    launch.biom.indic.random(nsim=100, pref.disp=1)

    # core fishing only
    launch.biom.indic.random(nsim=100, pref.disp=0, edge.fact=0)
    launch.biom.indic.random(nsim=100, pref.disp=0.5, edge.fact=0)
    launch.biom.indic.random(nsim=100, pref.disp=1, edge.fact=0)

    # edge fishing only
    launch.biom.indic.random(nsim=100, pref.disp=0, core.fact=0)
    launch.biom.indic.random(nsim=100, pref.disp=0.5, core.fact=0)
    launch.biom.indic.random(nsim=100, pref.disp=1, core.fact=0)
}

rerun.aor.rnd <- FALSE
if(rerun.aor.rnd){

    launch.aor.random(0, nsim=10, rg.cv=0.3, ff.all=0.8)
    launch.aor.random(0.5, nsim=10, rg.cv=0.3, ff.all=0.8)
    launch.aor.random(1, nsim=10, rg.cv=0.3, ff.all=0.8)
    launch.aor.random(2, nsim=10, rg.cv=0.3, ff.all=0.8)
}


