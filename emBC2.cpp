#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emBC(NumericVector y, NumericMatrix gen, double df = 10, double R2 = 0.5, double Pi = 0.75){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  NumericVector d(p);
  NumericVector b(p);
  double vy = var(y);
  if(Pi>0.5){Pi = 1-Pi;} 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));
  }
  double MSx = sum(vx)*Pi*(1-Pi);
  double Sa = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double ve = Sa;
  double va = Se;
  double Lmb = ve/va;
  double b0,b1,LR,eM,h2,C;
  double Pi0 = (1-Pi)/Pi;
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      e = e - gen(_,j)*(b[j]-b0);
    }
    ve = (sum(e*e)+Se)/(n+df);
    va = (sum(b*b)+Sa)/(p+df)/(Pi*(1-Pi));
    Lmb = ve/va;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(gen(k,_)*b)+mu;}
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("d") = d,
                      Named("hat") = fit,
                      Named("Vg") = va*MSx,
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}
