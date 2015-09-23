## Theortcl-model_Range-contrxn_Layout.r
## Starting basic version of range contraction population
## dynamics model
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: June  3, 2014
## Time-stamp: <2015-09-15 10:23:52 lauratb>
require(colorspace); require(RColorBrewer)
require(Rcpp)

tm.spatial.dyn <- function(emig.base=0.1, emig.max=emig.base,
                           r.growth.core=0.1, r.growth.edge=r.growth.core,
                           K.core=10000, K.edge=K.core,
                           fish.fact=0, fish.fact.edge=fish.fact,
                           fishing.start = 500,
                           pref.disp=0, add.r.pref=TRUE,
                           grid.width = 5) {
message("Setting seed!")
set.seed(999)

    ## Model parameters
    grid.width <<- grid.width
    ncell <<- grid.width^2 ## Number of cells

    ## Population parameters
    r.growth <<- rep(r.growth.core, ncell) # growth rate
    r.mrt <<- r.growth # mortality rate

    Ni <<- 100 # starting population initial

    ## get environment layout
    hab.setup.funk <- "calc.core.size" # or calc.core.size for a central core
core.layout.old <<- do.call(hab.setup.funk, list(ncell))
    core.layout <<- do.call(hab.setup.funk, list(ncell, r.growth.core, r.growth.edge))
#    core.layout <<- make.mvt.layout(grid.width, r.growth.core, r.growth.edge)

    ### Set up habitat
    habtype <<- "core" # should be 'core','even', 'random'
    if(habtype=="random") K <- rpois(ncell, K.core)
    if(habtype=="core") {
    K <<- rep(K.core, ncell)

    #r.growth[core.layout$edge.cells] <<- r.growth.edge
    r.growth[core.layout$core.cells] <<- core.layout$core.r.vals
    r.growth[core.layout$edge.cells] <<- core.layout$edge.r.vals
    #r.mrt <<- r.growth # mortality rate
    #r.mrt[core.layout$edge.cells] <- r.growth[core.layout$edge.cells]
}

if(habtype=="even") K <- rep(10000, ncell)
K.order <- numeric(); K.order[order(K)] <- 1:ncell

col.mat <<- matrix(NA, nrow=grid.width)

    core.colpal <- colorRampPalette(c("tomato1","tomato2","tomato3","tomato4"))
    edge.colpal <- colorRampPalette(c("turquoise1","turquoise3","royalblue1","royalblue3"))
#    core.colpal <- colorRampPalette(rev(c("tomato1","tomato3")))
#    edge.colpal <- colorRampPalette(rev(c("aquamarine1","aquamarine3")))

#    col.mat[core.layout$core.cells] <<- rep(core.colpal(2), each=2*grid.width)#length(core.layout$core.cells))
#    col.mat[core.layout$edge.cells] <<- rep(edge.colpal(2), each=2*grid.width)#length(core.layout$edge.cells))

    col.mat[core.layout$core.cells] <<- core.colpal(length(core.layout$core.cells))
    col.mat[core.layout$edge.cells] <<- edge.colpal(length(core.layout$edge.cells))

    col.mat.transp <<- col.mat
    col.mat.transp[] <<- col2transp(col.mat.transp)

    #############
    ts.max <<- 1000

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
    envpop$catch.mat <- array(NA, dim=dim(envpop$mat))
    envpop$F.array <- array(NA, dim=dim(envpop$mat))

    # record number of immigrants/emigrants to/by cell at every time step
    envpop$immig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$emig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$nkmat <- matrix(NA, nrow=ts.max, ncol=ncell)
    envpop$immig.rate <- array(NA, dim=c(ncell, ncell, ts.max))

###################################################
# DEFINE NEIGHBOURS BY CELL
 cell.index <<- data.frame(index=1:ncell,
                           xx=1:grid.width,
                           yy=rep(1:grid.width,each=grid.width))
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
    neighbour.mat <<- t(sapply(1:ncell, function(x) pb.neighbour(x, grid.width)))
    num.neigh <- apply(neighbour.mat, 1, sum)
    neighbycell <<- t(apply(neighbour.mat==1, 1, which))

}else{

    dmat <- as.matrix(dist(cell.index[,c("xx","yy")]))
    dmat[dmat==0] <- 99 # so that the cell itself isn't counted as a neighbour
    nbycell <- apply(dmat<1.5, 1, which) # neighbours of each cell (base cells are in row)
    num.neigh <- sapply(nbycell, length)
    neighbour.mat <<- sapply(1:nrow(dmat), function(x) as.numeric(dmat[x,]<1.5))
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

  # for cells with lower r.growth, adjust K to effective K value
  use.K.eff <- TRUE
  if(use.K.eff) K.cell <- K.cell*r.growth/r.mrt
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
ibc.nkratio <<- function(ts.now, nmat, Kcell=K, pref.disp=1, add.r=FALSE,
                        pref.calc="expo") {

  # When weighing preference, use equilibrium edge K to scale
  # preference
  use.K.edge.pref <- TRUE
  if(use.K.edge.pref) { Kcell <- K.core*r.growth/r.mrt
                        Kcell[Kcell == 0] <- 0.1*K.core }

  # take 1 - N/K ratio for each neighbour -- nkmat becomes a vector
  # of the same length as object neighbycell
  nkmat <- exp(1-nmat/Kcell)

  # ... if r counts in preferred dispersal, add r in product instead
  if(add.r) nkmat <- (r.growth*exp(1-nmat/Kcell))
  envpop$nkmat[ts.now,] <- nkmat
  nkmat <- nkmat[c(neighbycell)] # add value for each cell ID in neighbour matrix

  # set minimum dispersal rate between cells:
  # ** note that this results in a sharp change once edge cells reach
  # ** x% of their effective carrying capacity
  # if(any(nkmat<(min(r.growth)*emig.base))) {
  #  nkmat[nkmat<(min(r.growth)*emig.base)] <- min(r.growth)*emig.base # this allows a minimum of flow between cells
  #  nkmat[nkmat<(min(r.growth)*0.05)] <- min(r.growth)*0.05 # this allows a minimum of flow between cells
  #  print(ts.now)}
  # and stabilises population dynamics
  # but probably not the best way to do this

  # switch back to matrix
  dim(nkmat) <- dim(neighbycell)

  # if pref.disp=0, divide by identity -> no preferential dispersal
  # if pref.disp=1, divide by 1, weighted preferential dispersal towards
  # cells with lowest N/K ratio
    if(pref.calc=="prod") {
           nkweight <- 1+(nkmat - 1)*(1-pref.disp)
           nkmat <- nkmat/nkweight # apply weight to nk matrix
       }
    if(pref.calc=="expo") nkmat <- nkmat^pref.disp

  # remove dispersal from one cell at random, with higher probability for poor cells
  add.random <- FALSE
  if(add.random) {
    numcell0 <- 1
      sample.1 <- function(x) {
        ss <- sample(1:length(x), numcell0)
        }
    cell.0.disp <- apply(nkmat, 1, sample.1)
    cell.ind.0 <- arrayInd.rev(rep(1:nrow(nkmat),each=numcell0), c(cell.0.disp), .dim=dim(nkmat))
    nkmat[cell.ind.0] <- 0
    }

  # nkmat is ncell x 8 (num neighbours matrix)
  # containing absolute attractiveness index
  # first standardize by dividing by sum index over all cell neighbours
  # (senders in rows, neighbours in columns)
    relprop <- nkmat/apply(nkmat, 1, sum)
  if(ts.now == 2) relprop2 <<- relprop

  # assign these values to 0/1 ncell x ncell matrix
  # neet to transpose both neighbour.mat and relprop
  # so that we can fill in neighbours by cell first
  # (not sure this is necessary but too lazy to check)
    nbmat <- t(neighbour.mat)
    nbmat[nbmat==1] <- c(t(relprop))

    t(nbmat) # transpose back
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

    F.array <<- array(0, dim=c(grid.width, grid.width, ts.max))
    F.array[,,fishing.start:ts.max] <<- fish.fact*max(r.growth)
    mat.ind.edge <- arrayInd(core.layout$edge.cells, .dim=c(grid.width,grid.width))
    array.ind.edge <- arrayInd.rev(rep(mat.ind.edge[,1],ts.max-fishing.start+1),
                               rep(mat.ind.edge[,2],ts.max-fishing.start+1),
                               rep(fishing.start:ts.max, each=nrow(mat.ind.edge)),
                    .dim=c(grid.width, grid.width, ts.max))
    F.array[array.ind.edge] <<- fish.fact.edge*max(r.growth)
    envpop$F.array <- F.array # save in envpop

    ###################################################
    ###################################################
    ## Launch the simulation and store update in envpop environment
    ffunk <- formals(tm.spatial.dyn)
    fcall <- sys.call()
    ffunk[names(fcall)[-1]] <- fcall[-1]
    envpop$run.info <- ffunk

    cell.dyn(pref.disp=pref.disp, add.r=add.r.pref, emig.bef=TRUE)
    invisible(envpop)
}


if(!exists("everything.in")) {

  source("Theortcl-model_Range-contrxn_Figures-base.r")
  source("Theortcl-model_Multivar-normal-habitat.r")
  source("Theortcl-model_Range-contrxn_CellLayout-setup.r")
  source("Theortcl-model_Range-contrxn_Cpp.r")
  source("Theortcl-model_Range-contrxn_Cpp-bloops.r")
#  source("Theortcl-model_Range-contrxn_Scenarios.r")
#  source("Theortcl-model_Range-contrxn_BiomassIndic.r")
  everything.in <- TRUE
}




