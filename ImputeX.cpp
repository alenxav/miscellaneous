#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix IMP(NumericMatrix X){
  int p = X.ncol(); int n = X.nrow();
  LogicalVector MIS(n); NumericVector x(n);
  NumericVector z; double EXP;
  for(int j=0; j<p; j++){
    if(is_true(any(is_na(X(_,j))))){
      x = X(_,j); MIS = is_na(x);
      z = x[!MIS]; EXP = mean(z);
      X(_,j) = ifelse(MIS,EXP,x);}
  };return(X);};
