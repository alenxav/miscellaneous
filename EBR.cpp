#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP EBR(NumericVector y, NumericMatrix gen){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  NumericVector d(p);
  NumericVector b(p);
  double vy = var(y);
  NumericVector xx(p);
  for(int i=0; i<p; i++){xx[i] = sum(gen(_,i)*gen(_,i));}
  double MSx = mean(xx);
  double Sb = 0.5*vy/MSx;
  double ve = vy;
  NumericVector vb = b+Sb;
  NumericVector Lmb = ve/vb;
  double mu = mean(y);
  NumericVector e = y-mu;
  double b0,b1,eM,h2;
  for(int i=0; i<it; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e = e-gen(_,j)*(b1-b0);
      b[j] = b1;
      vb[j] = (b[j]*b[j]+Sb+ve/(Lmb[j]+xx[j]))/2;
      e = e - gen(_,j)*(b1-b0);}
    ve = sum(e*y)/(n-1);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;}
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(gen(k,_)*b)+mu;}
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Vb") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);}
