#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
void WGT(NumericMatrix X, NumericVector w){
  int p = X.ncol(); int n = X.nrow();
  for(int j=0; j<p; j++){
    for(int i=0; i<n; i++){
    X(i,j) = X(i,j)*w[i];};};}