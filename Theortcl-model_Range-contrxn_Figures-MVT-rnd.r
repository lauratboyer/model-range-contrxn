
theme_set(theme_bw())
ts.avrg.N <- function(run.obj) {

    nts <- run.obj$run.info$ts.max
    ncells <- run.obj$run.info$grid.width^2
    Kval <- run.obj$run.info$K.core
    df0 <- run.obj$dfN

    check.dev.size(8,8)
    df1 <- df0 %>% group_by(ts, cell.type) %>% summarize(mean.N=mean(y))
    ggplot(dat=df1, aes(x=ts, y=mean.N, colour=cell.type)) +
        geom_line() + geom_line(data=df0, aes(x=ts, y=y, group=id)) + geom_line(size=2) +
        geom_hline(y=Kval) + geom_hline(y=2000) + coord_cartesian(ylim=c(0,12000))

}

ts.avrg.simu <- function(multi.obj) {

    ns <- length(multi.obj[[1]])#$Nmat)
    df0 <- multi.obj$df#do.call(rbind, lapply(multi.obj, "[[", "dfN"))
    df0$simu.id <- rep(1:length(multi.obj), each=ns)

   # df1 <- df0 %>% group_by(simu.id, cell.type, ts) %>%
    #    summarize(mean.n=mean(y)) %>% mutate(simu.cell=paste(simu.id, cell.type))

    df2 <- df0 %>% group_by(cell.type, ts) %>%
        summarize(mean.n=mean(y), qt25=quantile(y, 0.25), qt75=quantile(y,0.75))
    #        geom_line(data=df1, aes(x=ts, y=mean.n, group=simu.cell, color=cell.type))+

    ggplot(data=df2, aes(x=ts, color=cell.type)) + geom_hline(yint=10000) +
            geom_ribbon(aes(ymin=qt25, ymax=qt75, fill=cell.type), colour=NA, alpha=0.25) +
        geom_line(aes(y=mean.n), alpha=0.95, size=2)# +

#                coord_cartesian(ylim=c(0, 1.2*10000))

}

fig3.explo <- function() {

    fa <- tm.spadyn.rd(emig.base=0, r.growth.edge=0.02, fish.fact=0.5, pref.disp=0, use.mvt=FALSE)
    fb <- tm.spadyn.rd(emig.base=0.1, r.growth.edge=0.02, fish.fact=0.5, pref.disp=0, use.mvt=FALSE)
    fc <- tm.spadyn.rd(emig.base=0.1, r.growth.edge=0.02, fish.fact=0.5, pref.disp=1, add.r.pref=TRUE,
                       use.mvt=FALSE)
    grid.arrange(ts.avrg.N(fa), ts.avrg.N(fb), ts.avrg.N(fc), nrow=1)

    list(fa, fb, fc)
}

fig4.explo <- function(which.type="ee",disp.type="evdisp") {

    # need df with
    # no core-edge structure
    ref.dat <- list(ref.even.even = get(paste0(which.type,".lowF.even.nodisp"))$df,
                    ref.even.edge = get(paste0(which.type,".lowF.edge.nodisp"))$df,
                    ref.even.core = get(paste0(which.type,".lowF.core.nodisp"))$df)
    all.dat <- list(evF.even.even = get(paste0(which.type,".lowF.even.",disp.type))$df,
                    cF.even.core = get(paste0(which.type,".lowF.core.",disp.type))$df,
                    eF.even.edge = get(paste0(which.type,".lowF.edge.",disp.type))$df)

    make.smry <- function(xlab, type="ref") {
        x <- get(paste0(type, ".dat"))[[xlab]]
        print(type);
        print(xlab)
        x0 <- x %>% filter(ts %in% 490:499) %>% group_by(cell.type) %>% summarize(K.eff=mean(y))
        x1 <- x %>% filter(ts %in% 990:999) %>% group_by(cell.type) %>%
            summarize(mean.N=mean(y), sd.N=sd(y)) %>% inner_join(x0) %>%
                mutate(BB0=mean.N/K.eff, sd.top=mean.N+sd.N, sd.bot=mean.N-sd.N, type=gsub(".*(....)$","\\1",xlab))
        #                mutate(scen=xlab, prop=(mean.N-ref.N)/K.eff)

        names(x1)[2:7] <- paste0(c("N.","sd.","K.","BB0.","sdtop.","sdbot."),type)
        x1

    }

    op2 <- function(x) x %>% mutate(deltBB0=(BB0.all-BB0.ref),
                                    prop=(N.all-N.ref)/N.ref,
                                    prop.top=(sdtop.all-sdbot.all)/K.ref)

    df1 <- inner_join(make.smry("ref.even.even"), make.smry("evF.even.even","all")) %>% op2
    df2 <- inner_join(make.smry("ref.even.core"), make.smry("cF.even.core","all")) %>% op2
    df3 <- inner_join(make.smry("ref.even.edge"), make.smry("eF.even.edge","all")) %>% op2
    df.even.evdisp <- rbind(df1, df2, df3)

    # define factor order for ggplot:
    df.even.evdisp$type <- factor(df.even.evdisp$type, levels=c("even","core","edge"), ordered=TRUE)

    hs.lab.vect <- c(ee="Even", ce="Core-edge")
    hs.lab.vect <- factor(hs.lab.vect, levels=hs.lab.vect, ordered=TRUE)
    disp.lab.vect <- c(evdisp="Even dispersal", prefdisp="Preferential dispersal")
    disp.lab.vect <- factor(disp.lab.vect, levels=disp.lab.vect, ordered=TRUE)
    df.even.evdisp$hs.lab <- hs.lab.vect[which.type]
    df.even.evdisp$disp.lab <- disp.lab.vect[disp.type]
    df.even.evdisp

}

fig4.ggp <- function(gdat=histofig.dat) {
    # make ggplot
    type.vect <- c(even="Even", core="Core only", edge="Edge only")
    gdat$type.vect <- type.vect[gdat$type]
    gdat$type.vect %<>% factor(levels=type.vect, ordered=TRUE)
    g1 <- ggplot(data=gdat, aes(x=type.vect, y=100*deltBB0, fill=cell.type)) +
        geom_hline(aes(yintercept=c(-0.25,0,0.25)), colour="grey") +
            facet_grid(hs.lab ~ disp.lab) +
             geom_bar(stat="identity", width=0.65, position=position_dodge(0.75)) + geom_point(colour=NA) +
                  scale_fill_manual(values=c("tomato","royalblue"),name="Cell type",label=c("Core","Edge"))+#, guide="none") +
                      theme(text=element_text(family="Segoe UI Light", size=12),
                            panel.grid=element_blank(),
                            axis.title.x=element_text(vjust=-0.5),
                            axis.title.y=element_text(vjust=1),
                            panel.margin=unit(0.4,"cm"), strip.background = element_rect(fill="slategrey",colour="slategrey"),
                            strip.text=element_text(colour="white", size=12.5, vjust=0.55),
                            legend.key=element_rect(colour="white"),
                            legend.key.size=unit(0.75,"cm")) +

                          coord_cartesian(ylim=c(-50,50)) + xlab("F by region") +
                              ylab("% change in biomass status Bcur/B0")
            #                  guides(fill=guide_legend(override.aes=list(fill=c("green","green"),size=5),size=5),
             #                        colour=guide_legend(override.aes=list(colour=c("green","green"), shape=15)))
#                                                                 fill="white",shape=15, size=2)))
    print(g1)
    ggsave("tm-range-dyn_histo-comp-scen_Bdecl_RND-DISTRIB.pdf", g1, width=8, height=7)
}

if(!exists("histofig.dat")) {
a1 <- fig4.explo("ee", "evdisp")
a2 <- fig4.explo("ee", "prefdisp")
a3 <- fig4.explo("ce", "evdisp")
a4 <- fig4.explo("ce", "prefdisp")
histofig.dat <- rbind(a1, a2, a3, a4)
}
