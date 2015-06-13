## Theortcl-model_Multivar-normal-habitat.r
## Defining core-edge identity by cell based on bivariate
## normal distribution, assigning r as normal random deviates
## for each cell based on mean r core/edge
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: October 18, 2014
## Time-stamp: <2015-06-14 09:46:57 lauratb>
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
  cov.val <- cov.val*sample(c(-1,1),1)

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
  dev.copy(CairoPNG, file=sprintf("Sample-rd-map-layout_CE-rtio%s_C%s_E%s.png",
                         core.ratio, core.r, edge.r), width=ww, height=hh,
           res=100, units="in"); dev.off()

}
