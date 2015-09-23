## Theortcl-model_Range-contrxn_Cpp-bloops.r
## Rewriting cpp version of model by running both loops in Cpp
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 25, 2014
## Time-stamp: <2015-09-21 10:25:33 lauratb>
require(inline)
require(Rcpp)
require(RcppArmadillo)
# these functions give number of rows/number of columns
#//int nr = Nmat.nrow();
#//int nc = Nmat.ncol();

#################################################
#################################################
cppFunction('
List alldyn(int inipop, int ncell, int numts,
NumericVector r_growth, NumericVector r_mort, NumericVector Kval,
double prop_emig, NumericMatrix Fval, bool add_r, double pref_val,
NumericMatrix neigh_mat) {

//#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
//using namespace arma;
/////////////////////
// Declare variables
// These store values at each ts
NumericMatrix Nmat(ncell, numts);
NumericMatrix immig_store(ncell, numts);
NumericMatrix emig_store(ncell, numts);
NumericMatrix rnk_vals(ncell, numts);
NumericMatrix catch_store(ncell, numts);

// These are used for dispersal calcs, overwritten at each ts
NumericMatrix neigh_store(ncell, 8);
NumericMatrix neigh_vals(ncell, ncell);
//arma::cube neighbour_mat(ncell, ncell, numts, false); // except for this one which archives the entire thing
//NumericMatrix neighbour_mat(ncell, ncell); // archives the entire thing at equilibrium

double Ntm1;
double n_emig;
double n_immig;
double NK_sum;
int nid;
NumericVector K_eff = Kval * r_growth/r_mort;
NumericVector immig_calc(ncell);
NumericVector emigcheck(8);
// Initialize
for (int i=0; i < ncell; i++) {
Nmat[i] = inipop;
emig_store(i,1) = inipop * prop_emig;
if(K_eff[i] == 0) K_eff[i] = 0.1 * Kval[i];
}

// Start loop through time steps
for(int ts = 1; ts < numts; ts++) {

///////////////////////////////////
// BEFORE updating population dynamics for ts,
// figure out where moving individuals are going

// Get relative distribution of target cells of migrants
// Loop through immig_calc to set NK values to use to
// compare neighbour attractiveness at this time-step
for(int cc = 0; cc < ncell; cc++) {
immig_calc[cc] = exp(1-Nmat(cc,ts-1)/K_eff[cc]);
if(add_r) immig_calc[cc] *= r_growth[cc];
immig_calc[cc] = pow(immig_calc[cc], pref_val);
rnk_vals(cc,ts) = immig_calc[cc];
}

for(int cc = 0; cc < ncell; cc++) { // for each cell

NK_sum=0; // reset sum of neighbour indices

// get NK for each neighbour of the cell and store sum to standardize later
for(int nn = 0; nn < 8; nn++) {
nid = neigh_mat(cc, nn)-1; // get neighbour cell number, conv to C++ indx
neigh_store(cc,nn) = immig_calc[nid]; // get NK
NK_sum += immig_calc[nid]; // sum all NKs
}


for(int nn = 0; nn < 8; nn++) { // for each neighbour of cell cc
nid = neigh_mat(cc, nn)-1 ; // get neighbour cell number
neigh_store(cc,nn) /= NK_sum; // standardize NK index to get proportion of migrants sent by cell cc to neighbours 1 to 8

//... and multiply by number of emigrants from the current cell cc...
// ... storing in column for receiving neighbour cell (to be summed up later)
neigh_vals(cc,nid) = emig_store(cc,ts)*neigh_store(cc,nn); // emig_store in no of indivs, done in previous time step
//neighbour_mat(cc,nid,ts)=neigh_vals(cc,nid); // store for subsequent analyses
if(ts==2) emigcheck(nn) = emig_store(nid,ts);

}
}

for(int cc = 0; cc < ncell; cc++) {

for(int nn = 0; nn < ncell; nn++) {
immig_store(cc,ts) += neigh_vals(nn,cc); // sum across columns to get immigrant per *receiving* cell
}

n_emig = emig_store(cc,ts);
n_immig = immig_store(cc,ts);

// Cell abundance at t - 1:
Ntm1 = Nmat(cc,ts-1);
// Fishing last
catch_store(cc, ts) = Ntm1 * Fval(cc,ts); // store catch for that step
Ntm1 *= 1-Fval(cc,ts);

// Start with emigration/immigration
Ntm1 = Ntm1 - n_emig + n_immig;



// Population grows
Ntm1 += r_growth[cc]*Ntm1 -
r_mort[cc]*Ntm1*Ntm1/Kval[cc];

Nmat(cc,ts) = Ntm1;

// Store the number of individuals emigrating for next step
if(ts < (numts-1)) emig_store(cc,ts+1) = prop_emig  * Nmat(cc,ts);
}
}

return Rcpp::List::create(Rcpp::Named("Nmat") = Nmat,
Rcpp::Named("emigmat") = emig_store,
Rcpp::Named("immigmat") = immig_store,
Rcpp::Named("neighmat") = neigh_vals,
Rcpp::Named("catchmat") = catch_store,
Rcpp::Named("Kcell") = K_eff,
Rcpp::Named("rnk") = rnk_vals,
Rcpp::Named("emigcheck") = emigcheck);
}')

cppFunction('
NumericVector arrayC(NumericVector input, IntegerVector dim) {
  input.attr("dim") = dim;
  return input;
}
')


