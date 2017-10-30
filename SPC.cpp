#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector SPC(NumericVector y,
                    NumericVector blk,
                    NumericVector row,
                    NumericVector col,
                    int rN=3, int cN=1){
  int n = y.size();
  NumericVector Cov(n), Phe(n), Obs(n);
  for(int i=0; i<n; i++){
    // Searching for neighbors
    for(int j=0; j<n; j++){
        if( (i!=j) & (blk[i]==blk[j]) & (abs(row[i]-row[j])<=rN) & (abs(col[i]-col[j])<=cN) ){
          Phe[j] = y[j]; Obs[j] = 1; }else{ Phe[j] = 0; Obs[j] = 0; }}
    // Averaging the neighbors
    Cov[i] = sum(Phe)/sum(Obs);
  }
  return Cov;
}
