#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP DER(NumericVector y, NumericMatrix gen, double cxx, NumericVector xx, int maxit = 40, double tol = 10e-6){
  int p = gen.ncol(); int n = gen.nrow(); double Ve,Va,b0,b1,eM; double mu = mean(y);
  NumericVector b(p), bc(p), Vb(p), fit(n), E(n); NumericVector e = y-mu;
  NumericVector Lmb = rep(cxx,p); int numit = 0; double cnv = 1;
  while(numit<maxit){ bc = b+0; for(int j=0; j<p; j++){ b0 = b[j];
    b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j)); b[j] = b1; e = e-gen(_,j)*(b1-b0);}
  eM = mean(e); mu = mu+eM; e = e-eM; Ve = sum(e*y)/(n-1); Vb = b*b+(Ve/(xx+Lmb));
  Lmb = sqrt(cxx*Ve/Vb); ++numit; cnv = sum(abs(bc-b)); if( cnv<tol ){break;};}
  for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b); E[k] = y[k]-fit[k];}; Va=sum(Vb);
  return List::create(Named("b")=b,Named("v")=Va,Named("h")=fit,Named("e")=E);}
