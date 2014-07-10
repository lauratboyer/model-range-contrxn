## Theortcl-model_Range-contrxn_Layout.r
## Starting basic version of range contraction population
## dynamics model
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: June  3, 2014
## Time-stamp: <2014-06-19 14:44:56 Laura>
require(colorspace)
require(Rcpp)
source("Theortcl-model_Range-contrxn_Cpp.r")
source("Theortcl-model_Range-contrxn_Figures-base.r")
## Model parameters
grid.width <- 5
ncell <- grid.width^2 ## Number of cells

## Population parameters
r.growth <- rep(0.05, ncell) # growth rate
r.mrt <- rep(0.05, ncell) # mortality rate

#r.growth[1:5] <- 0.025
#r.mrt[1:5] <- 0.025
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

    return(list(core.mat=core.mat,core.cells=core.ind,edge.cells=edge.ind))
}

core.layout <- calc.core.size(ncell)

habtype <- "core" # should be 'core','even', 'random'
if(habtype=="random") K <- rpois(ncell, 10000)
if(habtype=="core") {
    K <- rep(5000, ncell)
    K[core.layout$core.cells] <- 10000
    r.growth[core.layout$edge.cells] <- 0
}
if(habtype=="even") K <- rep(10000, ncell)
count <- function(x) length(unique(x))
K.order <- numeric(); K.order[order(K)] <- 1:ncell
colv <- c("dodgerblue3", "turquoise1", "tomato")
names(colv) <- sort(unique(K))
col.mat <- colv[as.character(K)]
col.mat <- matrix(col.mat, nrow=grid.width)
if(ncell==25) {
    col.mat[c(7,9,17,19)] <- "turquoise3"
    col.mat[c(1,5,21,25)] <- "dodgerblue1"
    col.mat[c(3,11,15,23)] <- "royalblue4"
}
if(ncell==100) col.mat[c(45,46,55,56)] <- "tomato"
col.mat.transp <- col.mat
col.mat.transp[] <- col2transp(col.mat.transp)
#############
ts.max <- 500

#############
## Defined fishing mortality (mouhaha)
fishing.start <- 250 # when fishing starts
F.array <- array(0, dim=c(grid.width, grid.width, ts.max))
F.array[,,fishing.start:ts.max] <- 0.5*max(r.growth)
mat.ind.edge <- arrayInd(core.layout$edge.cells, .dim=c(grid.width,grid.width))
array.ind.edge <- arrayInd.rev(rep(mat.ind.edge[,1],ts.max),
                               rep(mat.ind.edge[,2],ts.max),
                               rep(fishing.start:ts.max, each=ncell),
                    .dim=c(grid.width, grid.width, ts.max))
F.array[array.ind.edge] <- 0*max(r.growth)
# Set up environment and matrices to hold population values
# in recursive function call
if(!exists("envpop")) {
    envpop <- new.env() # create environment to store
    envpop$Ntvect <- rep(NA, ts.max)
    envpop$Ntvect[1] <- Ni
    envpop$mat <- array(NA, dim=c(grid.width, grid.width, ts.max))
    #randomNi <- rpois(ncell, Ni)
    randomNi <- rep(Ni, ncell)
    randomNi[core.layout$edge.cells] <- 0
    envpop$mat[,,1] <- randomNi
    #envpop$mat[,,1] <- c(rep(0,4),rep(Ni,5))
    envpop$emig.mat <- array(NA, dim=dim(envpop$mat))
    envpop$emig.mat[,,1] <- 0

    # record number of immigrants/emigrants to/by cell at every time step
    envpop$immig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$emig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$immig.rate <- array(NA, dim=dim(envpop$mat))
}

# Dispersal component of the model:
#
# Define emigration rate as a linear function of N/K
# where emigration rate is constant if N < K/2
# and increases linearly to 'emig.max' from N = K/2 to N = K
emig.base <- 0.1                                    # proportion of individuals that emigrate when N < K/2
emig.max <- 0.1 # maximum proportion of individuals that emigrate from cell at
                # each time step
emig.slope <- 2*(emig.max - emig.base)/K
emig.int <- emig.max-emig.slope*K

calc.emig.straight <- function(wN, es.cell=emig.slope, K.cell=K) { # emigration rate as a function of N
    emv <- es.cell * wN + emig.int
    emv[emig.slope == 0] <- emig.max # if no gradient set to emig.max
    emv[wN<(K.cell/2)] <- emig.base # if N < K/2, set to emigration baseline
    emv[wN > K.cell] <- emig.max # if N > K set to maximum emigration rate
    return(emv)
}

calc.emig.poly <- function(curv=2) {

    x.start <- 0
    x.end <- K
    y.base <- 0
    y.max <- 1

    div10 <- 10^curv

}



# set-up neighbours
 cell.index <- data.frame(index=1:ncell,
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

# habitat quality defines K
# emigration rates from cells are higher as N -> K (reacting to existing stimuli)
# sensing future conditions:
# (smart immigration) alternatively, immigration rates to cell are lower as N -> K

# Calculate no of immigrants by cell
ibc.nkratio <- function(nmat, Kcell=K, pref.disp=1, pref.calc="expo") {

    # take 1 - N/K ratio
    nkmat <- (1-nmat/Kcell)[c(neighbycell)]
    nkmat[nkmat<=emig.base] <- emig.base # this allows a minimum of flow between cells
#    nkmat[nkmat<=0.05] <- 0.05 # this allows a minimum of flow between cells
    # and stabilises population dynamics
    dim(nkmat) <- dim(neighbycell)

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

cell.dyn.slow <- function() { # Nt as a function of N_t-1

    # Reset population matrix
    envpop$mat[,,2:ts.max] <- NA

    upd.ts <- function(ts) {

        # calculate propotion of each cell that emigrates given current N
        prop.emig <- calc.emig.straight(envpop$mat[,,ts-1])
        # store it
        envpop$emig.mat[,,ts] <- prop.emig

        # calculate emigration matrix for this time-step
        # for each cell, N_t-1 * prop.emig
        emig.by.cell <- c(envpop$mat[,,ts-1]*prop.emig)
        # assign to each neighbour by multiplying on neighbour mat
        # and spreading evenly over neighbours by dividing by # neighbours
        # senders in rows, receivers in columns
        prop.immig <- ibc.nkratio(envpop$mat[,,ts-1], pref.disp=0)
        immig.by.cell <- apply(emig.by.cell*prop.immig, 2, sum)
        envpop$immig.store[ts,] <- immig.by.cell
        envpop$emig.store[ts,] <- emig.by.cell

        upd.cell <- function(wcell) {

            cp <- arrayInd(wcell, .dim=dim(envpop$mat))
            Ntm1 <- envpop$mat[cp[1],cp[2],ts-1] # N_t-1 for cell
            # calculate N_t from Ntm1, reproduction happens before emigration
            Ntp1 <- Ntm1 + Ntm1*r*(1-Ntm1/K[wcell]) -
                Ntm1*prop.emig[wcell] + immig.by.cell[wcell]
            # set maximum population size to K
            # (i.e. if number of migrants exceeds capacity, 'extra' mortality
            # gets rid of them -- check)
            if(Ntp1 < 0) Ntp1 <- 0
            #if(Ntp1 > K[wcell]) Ntp1 <- K[wcell]
            envpop$mat[cp[1],cp[2],ts] <- Ntp1


        }
        dmm <- sapply(1:ncell, upd.cell)
    }


    dmm <- sapply(2:ts.max, upd.ts)
}
