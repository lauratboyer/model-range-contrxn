## Theortcl-model_Range-contrxn_CellLayout-setup.r
## Theoretical model of range dynamics. Set-up cell layout for habitat quality.
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 25, 2014
## Time-stamp: <2015-09-10 15:05:16 lauratb>


## This one is for a middle core set-up
## Calculate width of square in the middle
## aiming between 30-40% of cells within core
calc.core.size <- function(n=ncell, r.growth.core=NA, r.growth.edge=NA) {

    gw <- sqrt(n)
    even.width <- (gw %% 2)==0
    core.width <- round(sqrt(0.35*n))
    even.core <-  (core.width %% 2)==0
    # for core to be perfectly centered square should be
    # even if range width is even and vice versa
    if(!identical(even.width, even.core)) core.width <- core.width-1
    edge.width <- 0.5*(gw - core.width)
    core.start <- edge.width + 1
    core.end <- core.start + core.width -1

    mat1 <- matrix(0,nrow=gw,ncol=gw)
    mat1[core.start:core.end,] <- 1
    core.mat <- mat1*t(mat1)
    core.ind <- which(core.mat==1)
    edge.ind <- which(core.mat==0)

    # calculate distance from center to set colours accordingly
    center.val <- median(1:grid.width)
    core.pos <- arrayInd(core.ind, .dim=c(grid.width, grid.width))
    edge.pos <- arrayInd(edge.ind, .dim=c(grid.width, grid.width))

    core.ind <- core.ind[order(-apply(abs(core.pos - center.val), 1, mean))]
    edge.ind <- edge.ind[order(-apply(abs(edge.pos - center.val), 1, mean))]

    return(list(core.cells=core.ind,edge.cells=edge.ind,
                core.r.vals=rep(r.growth.core, length(core.ind)),
                edge.r.vals=rep(r.growth.edge, length(edge.ind)), core.mat=core.mat))
}

## Band core set-up (like latitudes)
## Ideally want the same number of cells in core and bands,
## if that is not possible find band combination that centers it in the middle
## i.e. if grid.width is odd, center band + adjoining bands so long as # core < # edge
## if even grid.width, two middle ones + adjoining...
calc.core.band.size <- function(n=ncell) {

    gw <- sqrt(n) # grid.width
    even.width <- (gw %% 2)==0 # is grid.*width even
    # round to half x nearest multiple of 4 to calculate core width
    core.width <- 2*floor(gw/4) + ifelse(even.width, 0, 1)
    # if grid.width is odd, keep core width smaller than edge width
    if(core.width > 0.5*gw) core.width <- core.width - 2
    edge.width <- gw - core.width

    core.start <- 0.5*edge.width + 1
    core.end <- core.start + core.width - 1

    core.mat <- matrix(0,nrow=gw,ncol=gw)
    core.mat[,core.start:core.end] <- 1
    core.ind <- which(core.mat==1)
    edge.ind <- which(core.mat==0)

    # calculate distance from center to set colours accordingly
    center.val <- median(1:grid.width)
    core.pos <- arrayInd(core.ind, .dim=c(grid.width, grid.width))
    edge.pos <- arrayInd(edge.ind, .dim=c(grid.width, grid.width))

    core.ind <- core.ind[order(abs(core.pos[,2] - center.val))]
    edge.ind <- edge.ind[order(abs(edge.pos[,2] - center.val))]

    return(list(core.mat=core.mat,core.cells=core.ind,edge.cells=edge.ind))
}

