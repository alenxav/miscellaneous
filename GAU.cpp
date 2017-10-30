#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix GAU(NumericMatrix X){
  int n = X.nrow(); NumericVector D; NumericMatrix K(n,n); double d2, md;
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
  if(i==j){ K(i,j)=0; }else if(j>i){; D = X(i,_)-X(j,_);
  d2 = sum(D*D); d2 = d2*d2; K(i,j)=d2; K(j,i)=d2; }}}; md = mean(K);
  for(int i=0; i<n; i++){K(i,_) = exp(-K(i,_)/md);} return K;}