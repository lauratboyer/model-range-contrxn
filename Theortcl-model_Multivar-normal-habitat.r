## Theortcl-model_Multivar-normal-habitat.r
## Defining core-edge identity by cell based on bivariate
## normal distribution, assigning r as normal random deviates
## for each cell based on mean r core/edge
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: October 18, 2014
## Time-stamp: <2015-09-29 08:53:58 lauratb>
require(Cairo)
require(ggplot2)
require(extrafont)
require(mvtnorm)
require(colorspace)
require(data.table)

## test change
###########################################
###########################################
make.mvt.layout <- function(grid.width=5, r.growth.core=0.5, r.growth.edge=0.1,
                            rg.cv=0.2, edge.ratio=0.5) {

  dx <- 8/(grid.width-1)
  xv <- seq(-4,4,by=dx)
  yv <- xv

  make.single.dist <- function() {
  # generate random location for core center
  pos.x <- sample(xv[xv %between% c(-2,2)], 1)
  pos.y <- sample(yv[yv %between% c(-2,2)], 1)
  grd <- expand.grid(x=xv,y=yv)
  vl <- length(xv) #length of grid along one dim
  # generate density for multivariate normal
  sign.cor <- sample(c(-1,1),1)
  var1 <- sample(seq(1.5,4,by=0.1),1)
  var2 <- sample(seq(1.5,4,by=0.1),1)
  cov.prm <- sample(seq(0,2,by=0.2),1)
  m1 <- rnorm(1000, sd=sqrt(var1))
  m2 <- cov.prm*m1 + rnorm(1000, sd=sqrt(var2))
  var2 <- var(m2)
  cov.val <- cov(m1, m2)
  cov.val <- cov.val*sample(c(-1,1),1)

  sigma.mat <- matrix(c(var1,cov.val,cov.val,var2), ncol=2)

  a1 <- dmvnorm(grd, mean=c(pos.x, pos.y), sigma=sigma.mat)
  a1[is.infinite(a1)] <- 0
  attr(a1,"dim") <- c(vl,vl)
  d1 <- rmvnorm(5000, mean=c(pos.x, pos.y), sigma=sigma.mat)

  list(grid=a1,rdev=d1)
}

  m1 <- make.single.dist()
  m2 <- make.single.dist()
  mcomb <- m1$grid+m2$grid
  mcomb <- mcomb/max(mcomb)

  mcomb.stand <- mcomb
  ce.thresh <- quantile(mcomb, edge.ratio)
  core.r.vals <- rnorm(sum(mcomb>=ce.thresh),
                r.growth.core, sd=rg.cv*r.growth.core)
  edge.r.vals <- rnorm(sum(mcomb<ce.thresh),
                       r.growth.edge, sd=rg.cv*r.growth.edge)
  mcomb.stand[mcomb>=ce.thresh] <- core.r.vals
  mcomb.stand[mcomb<ce.thresh] <- edge.r.vals

  list(core.cells=which(mcomb>=ce.thresh),
       edge.cells=which(mcomb<ce.thresh),
       core.r.vals=core.r.vals,
       edge.r.vals=edge.r.vals)
}

# show single layout in ggplot
gg.draw.layout <- function(run.obj) { # from tm.spadyn.rd ()

    gw <- run.obj$run.info$grid.width
    ly.mat <- data.frame(id=1:(gw^2), cell.x=rep(1:gw, each=gw), cell.y=1:gw, r=runif(gw*gw))

    ly.mat$r[run.obj$core.layout$core.cells] <- run.obj$core.layout$core.r.vals
    ly.mat$r[run.obj$core.layout$edge.cells] <- run.obj$core.layout$edge.r.vals
    ggplot(dat=ly.mat, aes(x=cell.x, y=cell.y, fill=r)) + geom_tile(colour="black") + coord_equal() +
        scale_fill_continuous(low="gold", high="indianred3") +
            theme_bw() +
                theme(text=element_text(family="Segoe UI Light"), panel.grid=element_blank(), axis.title=element_blank(), axis.text=element_blank())
}

# show single layout in ggplot
gg.draw.layout.ts <- function(run.obj, ts2draw=1000) { # from tm.spadyn.rd ()

    dat0 <- run.obj$dfN
    gdat <- filter(dat0, ts==ts2draw)
    ggplot(dat=gdat, aes(x=cell.x, y=cell.y, fill=y)) + geom_tile(colour="black") + coord_equal() +
        scale_fill_continuous(low="gold", high="indianred3") +
            theme_bw() +
                theme(text=element_text(family="Segoe UI Light"), panel.grid=element_blank(), axis.title=element_blank(), axis.text=element_blank())
}

abund.colpal <- colorRampPalette(c("grey", "slategray2","royalblue","gold","hotpink3"))
gg.draw.layout.multi.ts <- function(run.obj, resp="y",
                                    ts2draw=c(1,10,25,50,75,100,250,499,500,550,600,650,700,750,900,1000)) { # from tm.spadyn.rd ()

    dat0 <- run.obj$dfN
    gdat <- filter(dat0, ts%in%ts2draw)
    gdat$in.range <- gdat$y > gdat$Bthresh.bycell

    gdat$resp <- gdat[,resp]
    g1 <- ggplot(dat=gdat, aes(x=cell.x, y=cell.y, fill=resp)) + geom_tile(colour="black") + coord_equal() +
        facet_wrap(~ts) +
                theme_bw() +
                theme(text=element_text(family="Segoe UI Light"),
                      panel.grid=element_blank(),
                      axis.title=element_blank(), axis.text=element_blank())
    if(resp=="y") return(g1 + scale_fill_gradientn(colours=abund.colpal(10)))
    if(resp=="in.range") return(g1 + scale_fill_manual(values=c("grey","royalblue")))
}

###########################################
###########################################
draw.mvt.layout <- function(core.ratio=0.5,
                            core.r=0.5, edge.r=0.1, cv=0.2,
                            num.centrs=2, grid.width=5) {

  dx <- 8/grid.width
  xv <- seq(-4,4,by=dx)
  yv <- xv

  make.single.dist <- function() {
  # generate random location for core center
  pos.x <- sample(xv[xv %between% c(-2,2)], 1)
  pos.y <- sample(yv[yv %between% c(-2,2)], 1)
  grd <- expand.grid(x=xv,y=yv)
  vl <- length(xv) #length of grid along one dim
  # generate density for multivariate normal
  sign.cor <- sample(c(-1,1),1)
  var1 <- sample(seq(1.5,4,by=0.1),1)
  var2 <- sample(seq(1.5,4,by=0.1),1)
  cov.prm <- sample(seq(0,2,by=0.2),1)
  m1 <- rnorm(1000, sd=sqrt(var1))
  m2 <- cov.prm*m1 + rnorm(1000, sd=sqrt(var2))
  var2 <- var(m2)
  cov.val <- cov(m1, m2)
  cov.val <- cov.val*sign.cor

  sigma.mat <- matrix(c(var1,cov.val,cov.val,var2), ncol=2)

  a1 <- dmvnorm(grd, mean=c(pos.x, pos.y), sigma=sigma.mat)
  a1[is.infinite(a1)] <- 0
  attr(a1,"dim") <- c(vl,vl)
  d1 <- rmvnorm(5000, mean=c(pos.x, pos.y), sigma=sigma.mat)

  list(grid=a1,rdev=d1)
}

  make.plot <- function() {

      m1 <- replicate(num.centrs, "[["(make.single.dist(),"grid"), simplify=FALSE)
      mcomb <- do.call("+", m1[1:min(2,num.centrs)])
      if(num.centrs==3) mcomb <- mcomb + m1[[3]]

      mcomb <- mcomb/max(mcomb)

  mcomb.stand <- mcomb
  ce.thresh <- quantile(mcomb, 1-core.ratio)
  mcomb.stand[mcomb>=ce.thresh] <- rnorm(sum(mcomb>=ce.thresh),
                  core.r,sd=core.r*cv)
  mcomb.stand[mcomb<ce.thresh] <- rnorm(sum(mcomb<ce.thresh),
                                        edge.r, sd=edge.r*cv)
  image(xv,yv,mcomb.stand , col=rev(heat_hcl(12)),asp=1,las=1,
        axes=FALSE)
  abline(h=xv-0.5*dx,lwd=0.5); abline(v=yv-0.5*dx,lwd=0.5)
  box()
}

  ww <- 7.25; hh <- 9.5
  check.dev.size(ww, hh)
  par(family="HersheySans", mfrow=c(5,4), mai=rep(0.15,4),
      omi=c(0.25,0.22,0.65,0.15))
  replicate(20, make.plot())

  mtext(sprintf("Core-edge ratio: %s, core r: %s, edge r: %s, CV: %s",
                core.ratio, core.r, edge.r, cv),
        outer=TRUE, adj=0, cex=1.2, line=1)
  other.figs <- FALSE
  if(other.figs) {
  plot(m1$rdev, asp=1, las=1, xlab="Y_1",ylab="Y_2",
       pch=19, col=col2transp("black",0.2))
  abline(0,1,col="dodgerblue")
  hist(m1$rdev[,1],xlim=c(-8,8),main="",ylim=c(0,1500))
  mtext(sprintf("Var = %s", round(var(m1$rdev[,1]),3)))
  abline(h=0)
  hist(m1$rdev[,2],xlim=c(-8,8),main="",ylim=c(0,1500))
  mtext(sprintf("Var = %s", round(var(m1$rdev[,2]),3)))
  abline(h=0)
}
  dev.copy(CairoPNG, file=sprintf("Sample-rd-map-layout_CE-rtio%s_C%s_E%s_%sx%s.png",
                         core.ratio, core.r, edge.r, grid.width, grid.width), width=ww, height=hh,
           res=100, units="in"); dev.off()

}
