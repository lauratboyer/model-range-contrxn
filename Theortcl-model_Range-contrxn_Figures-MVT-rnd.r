theme.lolo <- theme_bw() + theme(text=element_text(family="Segoe UI Light", size=12),
                                 axis.title.x=element_text(vjust=-0.615),
                                 axis.title.y=element_text(vjust=1.5),
                            panel.grid=element_blank(),
                            panel.margin=unit(0.4,"cm"), strip.background = element_rect(fill="slategrey",colour="slategrey"),
                                 strip.text=element_text(colour="white", size=12.5, vjust=0.4),
                                 legend.key=element_blank())
theme_set(theme.lolo)

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

    ggplot(data=df2, aes(x=ts, color=cell.type)) + geom_hline(aes(yintercept=c(0,10000))) +
        geom_ribbon(aes(ymin=qt25, ymax=qt75, fill=cell.type), colour=NA, alpha=0.25) +
            scale_fill_discrete(guide="none") +             scale_color_discrete(guide="none") +

                geom_line(aes(y=mean.n), alpha=0.95, size=2) + ylab("Abundance (in # individuals)") +
                                                         xlab("Time")

#                coord_cartesian(ylim=c(0, 1.2*10000))

}

fig3.explo <- function(fish.F=0.5) {


    frmfunk <- function(x) x %>% group_by(run.label, cell.type, ts) %>%
        summarize(y.median=median(y), y.25=quantile(y, 0.25), y.75=quantile(y, 0.75))
    ca <- spadyn.sims(list(emig.base=0, r.growth.edge=0.085,
                           fish.fact=fish.F, pref.disp=0, run.label="No movement, low HS"),10)$df %>% frmfunk
    cb <- spadyn.sims(list(emig.base=0.1, r.growth.edge=0.085,
                           fish.fact=fish.F, pref.disp=0, run.label="Uniform movement, low HS"),10)$df %>% frmfunk
    cc <- spadyn.sims(list(emig.base=0.1, r.growth.edge=0.085,
                           fish.fact=fish.F, pref.disp=1, add.r.pref=TRUE, run.label="Preferential movement, low HS"),10)$df %>% frmfunk

    fa <- spadyn.sims(list(emig.base=0, r.growth.edge=0.02,
                           fish.fact=fish.F, pref.disp=0,
                           run.label="No movement, high HS"), 10)$df %>% frmfunk
    fb <- spadyn.sims(list(emig.base=0.1, r.growth.edge=0.02,
                           fish.fact=fish.F, pref.disp=0,
                           run.label="Uniform movement, high HS"), 10)$df %>% frmfunk
    fc <- spadyn.sims(list(emig.base=0.1, r.growth.edge=0.02,
                           fish.fact=fish.F, pref.disp=1, add.r.pref=TRUE,
                           run.label="Preferential movement, high HS"), 10)$df %>% frmfunk

    fig3df <- rbind(ca, cb, cc, fa, fb, fc)
    attr(fig3df, "K") <- K[1]
    attr(fig3df, "r.growth.core") <- 0.1
    attr(fig3df, "r.growth.edge") <- c(0.085, 0.02)
    save(fig3df, file=paste0("tm_simu-data_ts-3-scen-two-envir_F", fish.F, ".RData"))
#    grid.arrange(ts.avrg.simu(ca), ts.avrg.simu(cb), ts.avrg.simu(cc),
 #                ts.avrg.simu(fa), ts.avrg.simu(fb), ts.avrg.simu(fc), nrow=2)
  #  ggsave("tm_ts-3-scen_FMx3_two-envir-scenarios_RND-DISTRIB.pdf")

    #ts.avrg.simu(ca)
    #list(fa, fb, fc)
    invisible(fig3df)
}

fig3.ggp <- function(fish.fact=0.5) {

    aa <- load(paste0("tm_simu-data_ts-3-scen-two-envir_F", fish.fact, ".RData"))
    fig3df$mov <- factor(sapply(strsplit(fig3df$run.label,", "), "[[", 1),
                         levels=c("No movement","Uniform movement","Preferential movement"), ordered=TRUE)
    fig3df$hs <- factor(sapply(strsplit(fig3df$run.label,", "), "[[", 2),
                        levels=c("low HS","high HS"), ordered=TRUE)
    K <- attr(fig3df,"K")
    r.mort <- attr(fig3df,"r.growth.core")
    r.edge <- attr(fig3df,"r.growth.edge")
    edge.line <- K * r.edge/r.mort; names(edge.line) <- c("low HS", "high HS")
    fig3df$edge.line <- edge.line[fig3df$hs]
    fig3df$note.x <- 1000
    fig3df$note.y <- 10000
    fig3df$note.lab <- "Edge K"
    fig3df$note.lab[fig3df$ts != 1 | fig3df$cell.type=="core"] <- ""
    fish.poly <- data.frame(x=c(0,500,500,0), y=c(0,0,12000,12000))

    colv <- c("tomato","dodgerblue")

    check.dev.size(12.1, 7)
    ggplot(data=fig3df, aes(x=ts, y=y.median, color=cell.type, fill=cell.type)) +
        # geom_polygon(data=fish.poly, aes(x=x, y=y), color=NA, fill="cornsilk4",alpha=0.2) +
        annotate("rect", xmin=0, xmax=500, ymin=0, ymax=12000, color=NA, fill="cornsilk4",alpha=0.2) +

        geom_hline(aes(yintercept=edge.line), colour="grey") + geom_hline(aes(yintercept=K), colour="grey") +
        geom_ribbon(aes(ymin=y.25, ymax=y.75), alpha=0.5, color=NA) + geom_line() +
            facet_grid(hs~mov) +
                geom_text(aes(x=note.x, y=edge.line, label=note.lab), family="Segoe UI Light",
                          size=4, hjust=1.2, vjust=-0.2, color="grey") +
         annotate("text", x=1000, y=10000, label="Core K", family="Segoe UI Light",
                          size=4, hjust=1.2, vjust=-0.2, color="grey") +
                scale_color_manual(name="Cell type", values=c("tomato","dodgerblue"),
                                   guide=guide_legend(override.aes=list(size=1.5))) +

                coord_cartesian(xlim=c(0,1000), ylim=c(0,12000)) + scale_fill_manual(values=colv, guide="none") +
                    ylab("Median abundance (+/- 25-75th quantile)") + xlab("Time")
ggsave("tm_ts-3-scen_FMx3_two-envir-scenarios_RND-DISTRIB.pdf")
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
