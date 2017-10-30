#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP MSX(NumericMatrix X){
  int p = X.ncol(); int n = X.nrow(); double m;
  NumericVector xx(p); NumericVector sx(p);
  for(int k=0; k<p; k++){ xx[k] = sum(X(_,k)*X(_,k));
    m = sum(X(_,k)); sx[k] = m*m/n; }
  double cxx = sum(xx-sx)/(n-1);
  return List::create(Named("MSx")=cxx,Named("xx")=xx);}