## Theortcl-model_Range-contrxn_Cpp-bloops.r
## Rewriting cpp version of model by running both loops in Cpp
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: July 25, 2014
## Time-stamp: <2014-11-05 17:09:58 Laura>

require(Rcpp)
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

#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
/////////////////////
// Declare variables
// These store values at each ts
NumericMatrix Nmat(ncell, numts);
NumericMatrix immig_store(ncell, numts);
NumericMatrix emig_store(ncell, numts);
NumericMatrix rnk_vals(ncell, numts);

// These are used for dispersal calcs, overwritten at each ts
NumericMatrix neigh_store(ncell, 8);
NumericMatrix neigh_store2(ncell, 8);
NumericMatrix neigh_vals(ncell, ncell);

double Ntm1;
double n_emig;
double n_immig;
double NK_sum;
int nid;
NumericVector K_eff = Kval * r_growth/r_mort;
NumericVector immig_calc(ncell);

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

for(int cc = 0; cc < ncell; cc++) {

NK_sum=0; // reset sum of neighbour indices

// get NK for each neighbour and store sum to standardize later
for(int nn = 0; nn < 8; nn++) {
nid = neigh_mat(cc, nn)-1; // get neighbour cell number, conv to C++ indx
neigh_store(cc,nn) = immig_calc[nid]; // get NK
NK_sum += immig_calc[nid]; // sum all NKs
}

for(int nn = 0; nn < 8; nn++) {
nid = neigh_mat(cc, nn)-1 ; // get neighbour cell number
neigh_store(cc,nn) /= NK_sum; // standardize
if(ts == 2) neigh_store2(cc,nn) = neigh_store(cc,nn);
//... and multiply by number of emigrants from native cell...
// ... storing in column for receiving neighbour cell
neigh_vals(cc,nid) = emig_store(nid,ts)*neigh_store(cc,nn);
}
}

for(int cc = 0; cc < ncell; cc++) {

for(int nn = 0; nn < ncell; nn++) {
immig_store(cc,ts) += neigh_vals(cc,nn);
}

n_emig = emig_store(cc,ts);
n_immig = immig_store(cc,ts);

// Cell abundance at t - 1:
Ntm1 = Nmat(cc,ts-1);

// Start with emigration/immigration
Ntm1 = Ntm1 - n_emig + n_immig;

// Fishing last
Ntm1 *= 1-Fval(cc,ts);

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
Rcpp::Named("Kcell") = K_eff,
Rcpp::Named("neighsto") = neigh_store2,
Rcpp::Named("rnk") = rnk_vals);

}')

cppFunction('
NumericVector arrayC(NumericVector input, IntegerVector dim) {
  input.attr("dim") = dim;
  return input;
}
')


