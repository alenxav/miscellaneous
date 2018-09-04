#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix SPM(NumericVector blk, NumericVector row, NumericVector col,
                  int rN=2, int cN=2){
  int n = blk.size();
  NumericMatrix X(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if( (blk[i]==blk[j]) & (i>j) & (abs(row[i]-row[j])<=rN) & (abs(col[i]-col[j])<=cN) ){
        X(i,j) = 1; X(j,i) = 1; }else{ X(i,j) = 0; X(j,i) = 0; }}
    X(i,i) = 0;}
  return X;}