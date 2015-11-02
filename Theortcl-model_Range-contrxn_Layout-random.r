tm.spadyn.rd <- function(emig.base=0.1, emig.max=emig.base,
                           r.growth.core=0.1, r.growth.edge=r.growth.core,
                           K.core=10000, K.edge=K.core,
                           fish.fact=0, fish.fact.edge=fish.fact,
                           rg.cv=0.1, ts.max=1000,
                         pref.disp=0, add.r.pref=TRUE,
                         use.mvt.range=TRUE,
                         run.label,
                           grid.width = 5) {

  ## Model parameters
  grid.width <<- grid.width
  ncell <<- grid.width^2 ## Number of cells
  Ni <<- 100 # starting population initial
  ts.max <<- ts.max
  fishing.start <- floor(0.5*ts.max)
  ## Population parameters
  r.growth <<- rep(r.growth.core, ncell) # growth rate
  r.mrt <<- r.growth # mortality rate

  ## get environment layout and set-up habitat by cells
  if(use.mvt.range) {
  core.layout <<- make.mvt.layout(grid.width,
                                  r.growth.core, r.growth.edge,
                                  edge.ratio=0.64, rg.cv=rg.cv)
}else{
    core.layout <<- calc.core.size(ncell, r.growth.core, r.growth.edge) # deterministic
}
  K <<- rep(K.core, ncell)
  r.growth[core.layout$core.cells] <<- core.layout$core.r.vals
  r.growth[core.layout$edge.cells] <<- core.layout$edge.r.vals

  # define colour palettes by region, reds are core, blues are edge
  # ordered to have greatest r with darker red/blue
  col.mat <<- matrix(NA, nrow=grid.width)
  core.colpal <- colorRampPalette(c("tomato1","tomato2","tomato3","tomato4"))
  edge.colpal <- colorRampPalette(c("turquoise1","turquoise3","royalblue1","royalblue3"))
  col.mat[core.layout$core.cells[order(core.layout$core.r.vals)]] <<-
    core.colpal(length(core.layout$core.cells))
  col.mat[core.layout$edge.cells[order(core.layout$edge.r.vals)]] <<-
    edge.colpal(length(core.layout$edge.cells))
  col.mat.transp <<- col.mat
  col.mat.transp[] <<- col2transp(col.mat.transp)

  # Set up environment and matrices to hold population values
  # in recursive function call

 #   envpop <<- new.env() # create environment to store

  #envpop$Ntvect <- rep(NA, ts.max)
   # envpop$Ntvect[1] <- Ni
    #envpop$mat <- array(NA, dim=c(grid.width, grid.width, ts.max))
    #envpop$mat[,,1] <- Ni
    #envpop$emig.mat <- array(NA, dim=dim(envpop$mat))
    #envpop$emig.mat[,,1] <- 0
    #envpop$catch.mat <- array(NA, dim=dim(envpop$mat))
    #envpop$F.array <- array(NA, dim=dim(envpop$mat))

    # record number of immigrants/emigrants to/by cell at every time step
    #envpop$immig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    #envpop$emig.store <- matrix(NA, nrow=ts.max, ncol=ncell)
    #envpop$nkmat <- matrix(NA, nrow=ts.max, ncol=ncell)
    #envpop$immig.rate <- array(NA, dim=c(ncell, ncell, ts.max))

  ###################################################
  # Define neighbours by cell
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

    }

  ######################################################
  ######################################################
  # Dispersal component of the model:
  #
  ######################################################
  ##### EMIGRATION FROM NATAL CELL. START.        ######
  # Note: currently not used in simulations, focusing on
  # DD dispersal at the immigration stage
  # (implemented directly in c++ code)(Nov 10 2014)

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

  ###################################################
  ## Fishing
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
  attr(F.array,"dim") <- c(ncell, ts.max) # switch from array to matrix for cpp

  ###################################################
  ###################################################
  ## Launch the simulation and store update in envpop environment
  ffunk <- formals(tm.spatial.dyn)
  fcall <- sys.call()
  ffunk[names(fcall)[-1]] <- fcall[-1]
  ffunk$ts.max <- ts.max
  pmat <- alldyn(Ni, ncell, numts=ts.max,
                 r.growth, r.mrt, K,
                 emig.max, F.array, add_r=add.r.pref, pref_val=pref.disp,
                 neighbycell)
  pmat$run.info <- ffunk
  pmat$core.layout <- core.layout
  pmat$cell.col <- list(op=col.mat,transp=col.mat.transp)

  map.enviro <- gg.draw.layout(pmat)
  pmat$map.layout <- map.enviro

  # in original envpop: F.array, Ntvect, catch.mat, emig.mat,
  # emig.store, immig.rate, immig.store, mat, nkmat, run.info

  # format dataframe for easy graphics with ggplot of N at each time-step by cell vs. cell.type
  dfN <- data.frame(ts=1:ts.max, y=c(t(pmat$Nmat)), id=rep(1:ncell, each=ts.max))
  if(!missing(run.label)) dfN$run.label <- run.label
  dfN$cell.type <- "edge"
  dfN$cell.type[dfN$id %in% core.layout$core.cells] <- "core"
  dfN <- inner_join(dfN, map.enviro$dat, by="id")

  # format dataframe of individual exchanges 'neighmat' for easy graphics
  exch0 <- format.exch.df(pmat$neighmat, core.layout)
  pmat$ss.regions.core2edge <- exch0$core2edge.ratio-1
  pmat$ss.cells.core <- srcsink.rd(pmat)$core.ss

  # add bcur/b0 metric
  cell.type.vect <- rep("edge", ncell)
  cell.type.vect[core.layout$core.cells] <- "core"
  pmat$pop.Bcur <- sum(pmat$Nmat[,ts.max])/sum(pmat$Nmat[,fishing.start-1])
  pmat$cell.Bcur <- data.frame(id=1:ncell, cell.type=cell.type.vect, pop.Bcur=pmat$pop.Bcur,
                               y=pmat$Nmat[,ts.max], y0=pmat$Nmat[,fishing.start-1],
                               Bcur=pmat$Nmat[,ts.max]/pmat$Nmat[,fishing.start-1])
  pmat$cell.Bcur$Bcur.diff <- with(pmat$cell.Bcur, (pop.Bcur-Bcur)/pop.Bcur)
  thresh.in.range <- 0.05 # threshold of B/B0 to consider cell occupied
  pmat$cell.Bcur$Bthresh.bycell <- with(pmat$cell.Bcur, thresh.in.range * y0) # cell-wise threshold
  pmat$cell.Bcur$Bthresh.all <- with(pmat$cell.Bcur, sum(y0)/ncell * thresh.in.range) # average range density threshold
  dfN <- inner_join(dfN, pmat$cell.Bcur[,c("id","Bthresh.bycell","Bthresh.all")], by="id")
  dfN$in.range <- as.numeric(dfN$Bthresh.bycell < dfN$y)
  pmat$dfN <- dfN
  pmat$ts.map <- gg.draw.layout.multi.ts(pmat) # now add map of abundance over time
  pmat$abund.ts <- tapply(dfN$y, dfN$ts, sum)
  pmat$area.ts <- tapply(dfN$in.range, dfN$ts, sum)
  invisible(pmat)
}

#############################################################
# format exchange dataframes
format.exch.df <- function(neighmat, cell.layout) {

    df.exch <- data.frame(sending.id=1:ncell, receiving.id=rep(1:ncell, each=ncell),
                        Nx=c(neighmat),
                          sending.cell.type="edge", receiving.cell.type="edge")
    df.exch$receiving.cell.type[df.exch$receiving.id %in% cell.layout$core.cells] <- "core"
    df.exch$sending.cell.type[df.exch$sending.id %in% cell.layout$core.cells] <- "core"

    exch.metr <- df.exch %>% group_by(sending.cell.type, receiving.cell.type) %>%
        summarize(N=sum(Nx)) %>% mutate(label=paste(sending.cell.type, receiving.cell.type, sep="2")) %>% data.frame
    rownames(exch.metr) <- exch.metr$label
    list(df=df.exch, core2edge.ratio=exch.metr["core2edge","N"]/exch.metr["edge2core","N"])

}

#############################################################
spadyn.sims <- function(..., nsim=100) {

    require(parallel)
    all.sims <- mclapply(1:nsim, function(i) do.call("tm.spadyn.rd",...))
    all.df0 <- mclapply(all.sims, "[[", "dfN")
    all.df <- do.call("rbind",all.df0)
    cell.out <- with(all.sims[[1]]$run.info, c(ts.max, grid.width))
    num.rows <- cell.out[1]*cell.out[2]^2
    all.df$sim.id <- rep(1:nsim, each = num.rows)
    list(all.sims=all.sims, df=all.df)
}
