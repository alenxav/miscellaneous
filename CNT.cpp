#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
void CNT(NumericMatrix X){for(int j=0;j<X.ncol();j++){X(_,j)=X(_,j)-mean(X(_,j));}}