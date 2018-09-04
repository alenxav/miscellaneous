#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix GRM(NumericMatrix X){
  int n = X.nrow(), p = X.ncol();
  NumericMatrix K(n,n); NumericVector xx(p); double zz;
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
        zz = sum( (X(i,_)-xx(i))*(X(j,_)-xx(j)) );
        K(i,j)=zz; K(j,i)=zz;}}
  for(int i=0; i<p; i++){xx[i] = K(i,i);}
  return K/mean(xx);}
