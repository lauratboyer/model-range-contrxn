## Theortcl-model_Range-contrxn_Layout.r
## Starting basic version of range contraction population
## dynamics model
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: June  3, 2014
## Time-stamp: <2014-07-18 17:39:44 Laura>
require(colorspace)
require(Rcpp)
source("Theortcl-model_Range-contrxn_Cpp.r")
source("Theortcl-model_Range-contrxn_Figures-base.r")

tm.spatial.dyn <- function(emig.base=0.1, emig.max=emig.base,
                           r.growth.core=0.05, r.growth.edge=r.growth.core,
                           K.core=10000, K.edge=K.core,
                           fish.fact=0, fish.fact.edge=fish.fact,
                           pref.disp=0, add.r.pref=FALSE,
                           grid.width = 5) {

  if(add.r.pref) pref.disp <- 2

    ## Model parameters
    grid.width <<- grid.width
    ncell <<- grid.width^2 ## Number of cells

    ## Population parameters
    r.growth <<- rep(r.growth.core, ncell) # growth rate
    r.mrt <<- r.growth # mortality rate

    Ni <- 100 # starting population initial
    set.seed(99)

    # if core habitat distribution, calculate
    # width of square in the middle
    # aiming between 30-40% of cells within core
calc.core.size <- function(n=ncell) {

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

    return(list(core.mat=core.mat,core.cells=core.ind,edge.cells=edge.ind))
}

    core.layout <- calc.core.size(ncell)

    ### Set up habitat
    habtype <<- "core" # should be 'core','even', 'random'
    if(habtype=="random") K <- rpois(ncell, 10000)
    if(habtype=="core") {
    K <<- rep(K.core, ncell)
    #[core.layout$edge.cells] <- K.edge
    r.growth[core.layout$edge.cells] <<- r.growth.edge
    #r.mrt <<- r.growth # mortality rate
    #r.mrt[core.layout$edge.cells] <- r.growth[core.layout$edge.cells]
}

if(habtype=="even") K <- rep(10000, ncell)
K.order <- numeric(); K.order[order(K)] <- 1:ncell

col.mat <<- matrix(NA, nrow=grid.width)

    core.colpal <- colorRampPalette(c("tomato1","tomato2","tomato3","tomato4"))
    edge.colpal <- colorRampPalette(c("turquoise1","turquoise3","royalblue3"))

    col.mat[core.layout$core.cells] <<- core.colpal(length(core.layout$core.cells))
    col.mat[core.layout$edge.cells] <<- edge.colpal(length(core.layout$edge.cells))

    col.mat.transp <<- col.mat
    col.mat.transp[] <<- col2transp(col.mat.transp)

    #############
    ts.max <<- 500

    # Set up environment and matrices to hold population values
    # in recursive function call

    envpop <<- new.env() # create environment to store
    envpop$Ntvect <- rep(NA, ts.max)
    envpop$Ntvect[1] <- Ni
    envpop$mat <- array(NA, dim=c(grid.width, grid.width, ts.max))
    #randomNi <- rpois(ncell, Ni)
    randomNi <- rep(Ni, ncell)
    #randomNi[core.layout$edge.cells] <- 0
    envpop$mat[,,1] <- randomNi
    #envpop$mat[,,1] <- c(rep(0,4),rep(Ni,5))
    envpop$emig.mat <- array(NA, dim=dim(envpop$mat))
    envpop$emig.mat[,,1] <- 0

    # record number of immigrants/emigrants to/by cell at every time step
    envpop$immig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$emig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$immig.rate <- array(NA, dim=c(ncell, ncell, ts.max))



###################################################
# DEFINE NEIGHBOURS BY CELL
 cell.index <<- data.frame(index=1:ncell,
                         xx=rep(1:grid.width,each=grid.width),
                         yy=1:grid.width)
boundary.opt <- "periodic"
if(boundary.opt == "periodic") {

######################################################
######################################################
## Find neighbour given lattice spatial structure with
## periodic boundaries
## wi is index of cell for which neighbours are wanted,
## n are matrix dimension
pb.neighbour <- function(wi, n=3) {

    # get position in array
    pos <- arrayInd(wi, .dim=c(n, n))
    nr <- pos[1]; nc <- pos[2]

    # get neighbours on sides, top and corners
    # starting from top left to bottom right
    ypos <- nc + c(-1, 0, 1, -1, 1, -1, 0, 1)
    xpos <- nr + rep(c(-1, 0, 1), c(3, 2, 3))

    # make boundaries periodic -- if pos = 0,
    # then neighbour is last cell (loops to end)
    # and vice versa
    xpos[xpos == 0] <- n
    xpos[xpos > n] <- 1
    ypos[ypos == 0] <- n
    ypos[ypos > n] <- 1

    # get matrix index corresponding to these positions
    aI <- arrayInd.rev(xpos, ypos, .dim=c(n,n))
    neighvect <- numeric(n*n)
    neighvect[aI] <- 1 # set position to 1 if it's a neighbour
    neighvect
}

    # neighbours of each cell (base cells are in row)
    # i.e. neighbours of cell i are in neighbour.mat[i,]
    neighbour.mat <- t(sapply(1:ncell, function(x) pb.neighbour(x, grid.width)))
    num.neigh <- apply(neighbour.mat, 1, sum)
    neighbycell <- t(apply(neighbour.mat==1, 1, which))

}else{

    dmat <- as.matrix(dist(cell.index[,c("xx","yy")]))
    dmat[dmat==0] <- 99 # so that the cell itself isn't counted as a neighbour
    nbycell <- apply(dmat<1.5, 1, which) # neighbours of each cell (base cells are in row)
    num.neigh <- sapply(nbycell, length)
    neighbour.mat <- sapply(1:nrow(dmat), function(x) as.numeric(dmat[x,]<1.5))
}

######################################################
######################################################
# Dispersal component of the model:
#
######################################################
##### EMIGRATION FROM NATAL CELL. START.        ######

# Define emigration rate as a linear function of N/K
# where emigration rate is constant if N < K/2
# and increases linearly to 'emig.max' from N = K/2 to N = K
emig.base <<- emig.base # proportion of individuals that emigrate when N < K/2
emig.max <<- emig.max # maximum proportion of individuals that emigrate from cell at
                # each time step
emig.slope <- 2*(emig.max - emig.base)/K
emig.int <- emig.max-emig.slope*K

# Linear relationship between emig.base and emig.max
calc.emig.straight <<- function(wN, es.cell=emig.slope, K.cell=K) { # emigration rate as a function of N
    emv <- es.cell * wN + emig.int
    emv[emig.slope == 0] <- emig.max # if no gradient set to emig.max
    emv[wN<(K.cell/2)] <- emig.base # if N < K/2, set to emigration baseline
    emv[wN > K.cell] <- emig.max # if N > K set to maximum emigration rate
    return(emv)
}

#########  EMIGRATION FROM NATAL CELL. END. #########
#####################################################

#####################################################
###### IMMIG TO NEIGHBOUR CELLS. START.      ########

# habitat quality defines K
# emigration rates from cells are higher as N -> K (reacting to existing stimuli)
# sensing future conditions:
# (smart immigration) alternatively, immigration rates to cell are lower as N -> K

# Calculate no of immigrants by cell
ibc.nkratio <<- function(nmat, Kcell=K, pref.disp=1, add.r=FALSE,
                        pref.calc="expo") {

    # take 1 - N/K ratio
    nkmat <- (1-nmat/Kcell)[c(neighbycell)]
    if(add.r) nkmat <- (r.growth*(1-nmat/Kcell))[c(neighbycell)]
    nkmat[nkmat<(min(r.growth)*emig.base)] <- min(r.growth)*emig.base # this allows a minimum of flow between cells
    # and stabilises population dynamics
    # but probably not the best way to do this

    dim(nkmat) <- dim(neighbycell) # switch back to matrix

    # if pref.disp=0, divide by identity -> no preferential dispersal
    # if pref.disp=1, divide by 1, weighted preferential dispersal towards
    # cells with lowest N/K ratio
    if(pref.calc=="prod") {
           nkweight <- 1+(nkmat - 1)*(1-pref.disp)
           nkmat <- nkmat/nkweight # apply weight to nk matrix
       }
    if(pref.calc=="expo") nkmat <- nkmat^pref.disp

    # standardize to send neighbours to cell based on
    # relative pref.disp weighted proportion
    # (neighbours in columns)
    relprop <- nkmat/apply(nkmat, 1, sum)
    nbmat <- t(neighbour.mat)
    nbmat[nbmat==1] <- c(t(relprop))
    t(nbmat)
}

# Calculate no of immigrants by cell -- including home cell in
# possible 'destination'
ibc.nkratio.wnatal <- function(nmat, Kcell=K, pref.disp=1, add.r=FALSE,
                        pref.calc="expo") {

    # add home cell to compare against neighbours
    neighmat <- cbind(1:nrow(neighbycell), neighbycell)
    # take 1 - N/K ratio
    nkmat <- (1-nmat/Kcell)[c(neighmat)]
    if(add.r) nkmat <- (r.growth*(1-nmat/Kcell))[c(neighmat)]
    nkmat[nkmat==0] <- r.growth*emig.base # this allows a minimum of flow between cells
    # and stabilises population dynamics
    # but probably not the best way to do this
    dim(nkmat) <- dim(neighmat) # switch back to matrix

    # if pref.disp=0, divide by identity -> no preferential dispersal
    # if pref.disp=1, divide by 1, weighted preferential dispersal towards
    # cells with lowest N/K ratio
    if(pref.calc=="prod") {
           nkweight <- 1+(nkmat - 1)*(1-pref.disp)
           nkmat <- nkmat/nkweight # apply weight to nk matrix
       }
    if(pref.calc=="expo") nkmat <- nkmat^pref.disp

    # add home advantage!
    home.adv <- 1
    nkmat[,1] <- nkmat[,1]*home.adv


    # standardize to send neighbours to cell based on
    # relative pref.disp weighted proportion
    # (neighbours in columns)
    relprop <- nkmat/apply(nkmat, 1, sum)
    relprop <- relprop[,-1] # remove home cell from proportions
    nbmat <- t(neighbour.mat)
    nbmat[nbmat==1] <- c(t(relprop))
    t(nbmat)
}

###### IMMIG TO NEIGHBOUR CELLS. END.      ########
###################################################


###################################################
###################################################
#######              FISHING!              ########

    # define F.array in global environment
    fishing.start <- 300 # when fishing starts
    F.array <<- array(0, dim=c(grid.width, grid.width, ts.max))
    F.array[,,fishing.start:ts.max] <<- fish.fact*max(r.growth)
    mat.ind.edge <- arrayInd(core.layout$edge.cells, .dim=c(grid.width,grid.width))
    array.ind.edge <- arrayInd.rev(rep(mat.ind.edge[,1],ts.max-fishing.start+1),
                               rep(mat.ind.edge[,2],ts.max-fishing.start+1),
                               rep(fishing.start:ts.max, each=nrow(mat.ind.edge)),
                    .dim=c(grid.width, grid.width, ts.max))
    F.array[array.ind.edge] <<- fish.fact.edge*max(r.growth)

    ###################################################
    ###################################################
    ## Launch the simulation and store update in envpop environment
    ffunk <- formals(tm.spatial.dyn)
    fcall <- sys.call()
    ffunk[names(fcall)[-1]] <- fcall[-1]
    envpop$run.info <- ffunk

    cell.dyn(pref.disp=pref.disp, add.r=add.r.pref)
}


