#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix ArcCosKernel(NumericMatrix X){
  int n = X.nrow(), p = X.ncol();
  NumericMatrix K(n,n);
  NumericVector xx(p), DiagZZ(p), DiagK(n); double zz, alpha,theta,tmp;
  for(int i=0; i<p; i++){ xx[i] = mean(X(_,i)); }
  for(int j=0; j<n; j++){ DiagZZ[j] = sum( (X(j,_)-xx)*(X(j,_)-xx) ) + 0.1;}
  for(int i=0; i<n; i++){
        zz = sum((X(i,_)-xx)*(X(i,_)-xx));
        theta = acos( zz  / sqrt(DiagZZ[i]*DiagZZ[i]) );
        tmp = (DiagZZ[i]*DiagZZ[i])/3.1416 * (sin(theta)+(3.1416-theta)*cos(theta));
        DiagK(i) = tmp;}
  alpha = 1/mean(DiagK);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(i==j ){ K(i,i) = DiagK(i)*alpha; }else if(i<j ){
          zz = sum((X(i,_)-xx)*(X(j,_)-xx));
          theta = acos( zz  / sqrt(DiagZZ[i]*DiagZZ[j]) );
          tmp = alpha * (DiagZZ[i]*DiagZZ[j])/3.1416 * (sin(theta)+(3.1416-theta)*cos(theta));
          K(i,j) = tmp; K(j,i) = tmp; }}}
  return K;}