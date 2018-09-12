#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix GRM(NumericMatrix X, bool Code012 = false){
  int n = X.nrow(), p = X.ncol();
  NumericMatrix K(n,n); NumericVector xx(p); double zz, Sum2pq=0.0;
  for(int i=0; i<p; i++){ xx[i] = mean(X(_,i)); }
  if(Code012){
    for(int i=0; i<p; i++){ Sum2pq = Sum2pq + xx[i]*xx[i]/2;}
  }else{
    for(int i=0; i<p; i++){ Sum2pq = Sum2pq + var(X(_,i));}}
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(i<=j ){
   zz = sum( (X(i,_)-xx(i))*(X(j,_)-xx(j)) );
   K(i,j)=zz; K(j,i)=zz;}
    }
  }
  return K/Sum2pq;
}
