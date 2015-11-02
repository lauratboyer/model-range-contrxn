
df0 <- data.frame(x = 1:4,
                  y = c(2, 0, 5, 8),
                  z = c(0.5, 1, 0.8, 1.08)) %>% mutate(z.mean= mean(z), z.dev=z-z.mean)

breakv <- seq(-0.5, 0.35, by=0.025)
df0$cutv <- cut(df0$z.dev, breakv)#, labels=FALSE)
colv <- colorRampPalette(brewer.pal(9,"RdYlBu"))(length(breakv))[df0$cutv]
colv <- c("royalblue","azure2","pink","indianred1")
print(ggplot(df0) +
    geom_point(aes(x=x, y=y, size=z, colour=cutv), pch=19) +
        scale_size_continuous(range=c(3, 20), guide="none") +
            scale_colour_manual(values=colv, name="deviation",
                                guide=guide_legend(override.aes=list(size=5))) +
                geom_point(aes(x=x, y=y, size=z.mean), pch=1) + theme(legend.key=element_blank())
      )
