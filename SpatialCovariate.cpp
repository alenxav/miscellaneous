#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector SPC(NumericVector y,
                   NumericVector blk, NumericVector row, NumericVector col,
                   int rN=2, int cN=2){
  int n = y.size(); NumericVector Cov(n), Phe(n), Obs(n);
  for(int i=0; i<n; i++){
    // Searching for neighbors
    for(int j=0; j<n; j++){
        if( (i>j) & (blk[i]==blk[j]) & (abs(row[i]-row[j])<=rN) & (abs(col[i]-col[j])<=cN) ){
          Phe[i] = Phe[i]+y[j]; Obs[i] = Obs[i]+1; Phe[j] = Phe[j]+y[i]; Obs[j] = Obs[j]+1; }}}
  // Averaging the neighbors
  Cov = Phe/Obs; return Cov;}
