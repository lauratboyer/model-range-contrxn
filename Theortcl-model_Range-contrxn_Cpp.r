require(Rcpp)

cppFunction('NumericMatrix updcellcpp(NumericMatrix inipopmat,
NumericVector rval_growth, NumericVector rval_mrt, NumericVector Kval,
NumericVector emigv, NumericVector immigv,
NumericVector fishmort) {
int nr = inipopmat.nrow();
int nc = inipopmat.ncol();
NumericMatrix popmatnow(nr,nc); //empty matrix

for (int cc = 0; cc < nc; cc++) {
for (int rr = 0; rr < nr; rr++) {
double Ntm1 = inipopmat[cc*nc+rr];
// calculate cell abundance before fishing
double Ntnow = Ntm1+rval_growth[cc*nc+rr]*Ntm1 -
rval_mrt[cc*nc+rr]*Ntm1*Ntm1/Kval[cc*nc+rr] -
emigv[cc*nc+rr]*Ntm1 + immigv[cc*nc+rr];
// now remove proportion based on F
popmatnow[cc*nc+rr] = Ntnow*(1-fishmort[cc*nc+rr]);
}
}
return popmatnow; }')


cppFunction('NumericMatrix updcellcpp_eb4r(NumericMatrix inipopmat,
NumericVector rval_growth, NumericVector rval_mrt, NumericVector Kval,
NumericVector emigv, NumericVector immigv,
NumericVector fishmort) {
int nr = inipopmat.nrow();
int nc = inipopmat.ncol();
NumericMatrix popmatnow(nr,nc); //empty matrix

for (int cc = 0; cc < nc; cc++) {
for (int rr = 0; rr < nr; rr++) {
double Ntm1 = inipopmat[cc*nc+rr];
// calculate cell abundance before fishing, emigration occurs before repr
double Ntnow = Ntm1 - emigv[cc*nc+rr]*Ntm1 + immigv[cc*nc+rr];
double Ntnowrep = Ntnow + rval_growth[cc*nc+rr]*Ntnow - rval_mrt[cc*nc+rr]*Ntnow*Ntnow/Kval[cc*nc+rr];
// now remove proportion based on F
popmatnow[cc*nc+rr] = Ntnowrep;
}
}
return popmatnow; }')


cell.dyn <- function(pref.disp=0, add.r=FALSE, use.home=FALSE, emig.before=TRUE) { # Nt as a function of N_t-1

    # Reset population matrix
    envpop$mat[,,2:ts.max] <<- NA

    upd.ts <- function(ts) {

        matnow <- envpop$mat[,,ts-1]
        # first step: fish and chips
        # using fishing mortality as defined in current timestep
        fishmat <- matnow*(1-F.array[,,ts])
        matnow <- fishmat

        # then calculate proportion of each cell that emigrates given
        # current, post-fishing N
        prop.emig <- calc.emig.straight(matnow)

        if(use.home){
            pi <- ibc.nkratio.wnatal(ts, matnow, pref.disp=pref.disp, add.r=add.r)
            prop.emig <- apply(pi, 1, sum) # emigration sum of relative proportions that go to neighbours

        }

        # store it
        envpop$emig.mat[,,ts] <- prop.emig

        # calculate emigration matrix for this time-step
        # for each cell, N_t-1 * prop.emig
        #emig.by.cell <- c(envpop$mat[,,ts-1]*prop.emig)
        emig.by.cell <- c(matnow*prop.emig)


        # assign to each neighbour by multiplying on neighbour mat
        # and spreading evenly over neighbours by dividing by # neighbours
        # senders in rows, receivers in columns
        prop.immig <- ibc.nkratio(ts, matnow, pref.disp=pref.disp, add.r=add.r)
        envpop$immig.rate[,,ts] <- prop.immig
        immig.by.cell <- apply(emig.by.cell*prop.immig, 2, sum)
        envpop$immig.store[ts,] <- immig.by.cell
        envpop$emig.store[ts,] <- emig.by.cell

        if(emig.before) {
        matupd <- updcellcpp_eb4r(matnow, r.growth, r.mrt, K,
                             prop.emig, immig.by.cell, c(F.array[,,ts]))
        }else{
        matupd <- updcellcpp(matnow, r.growth, r.mrt, K,
                             prop.emig, immig.by.cell, c(F.array[,,ts]))
    }

        # update value stored in population array
        envpop$mat[,,ts] <- matupd


    }

    dmm <- sapply(2:ts.max, upd.ts)
}

