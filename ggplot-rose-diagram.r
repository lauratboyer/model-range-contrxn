df0 <- data.frame(x=rep(1:4, each=4),
                  y=rep(1:4, each=4),
                  id=rep(1:4, each=4),
                  z=runif(16))

df0 <- data.frame(x=c(1,2,2,3,3,3,4,4,4,4,4,4,5,5,5,5,6,6,6,7,7,7,8,8,8,8))
g0 <- ggplot(dat=df0, aes(x=factor(x))) +
          geom_bar(colour="black", fill="grey", width=1) + coord_polar()

df1 <- data.frame(x=1:8, y=c(1,2,3,6,4,3,3,4)/26, reg=rep(c("edge","core"), each=4))
df1 <- rbind(df1, df1, df1, df1, df1, df1, df1, df1, df1, df1, df1, df1)
df1$id <- rep(1:4, each=8)
df1$id2 <- rep(1:3)
g1 <- ggplot(dat=df1, aes(x=factor(x), y=y, fill=reg)) +
    geom_bar(stat="identity", colour="white", width=1) + coord_polar(start=-pi/8) +
        facet_grid(id2~id) + theme(strip.background=element_blank(),
                                   panel.margin=unit(-0.1,"cm"),
                                   text=element_blank(),
                                   axis.ticks=element_blank(),
                                   panel.border=element_blank()) +
            guides(fill="none")

df2 <- data.frame(sender=1:25, id=rep(1:25,each=25), val=c(h0$neighmat))
df2 %<>% filter(val >0)
df2$neigbour <- c(apply(h0$neighmat>0, 2, which))
ind.df <- unique(h0$dfN[,c("cell.x","cell.y","cell.type","id")])
df2 %<>% merge(ind.df[,c("cell.x","cell.y","id")], by.x="sender", by.y="id")
df2 %<>% merge(ind.df[,c("cell.type","id")], by="id")

g1 <- ggplot(dat=df2, aes(x=id, y=val, fill=cell.type)) +
    geom_bar(stat="identity", colour="white", width=1) + #coord_polar(start=-pi/8) +
        facet_grid(cell.x~cell.y) + theme(strip.background=element_blank(),
                                   panel.margin=unit(-0.1,"cm"),
                                   text=element_blank(),
                                   axis.ticks=element_blank(),
                                   panel.border=element_blank()) +
            guides(fill="none")
