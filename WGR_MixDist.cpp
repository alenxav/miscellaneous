#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP MixDist(NumericVector y, NumericMatrix X, double pi=0.7){
  // Regularization settings
  double alpha=0.05;
  double iter=350;
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and variances
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx= sum(vx);
  // Get regularization parameters
  double mu = mean(y);
  double L1 = 0.5*MSx*(sd(y)*alpha);
  double L2 = 0.5*MSx;
  // Create empty objects
  double b0,b1,b2,eM,C,cj,dj,pj;
  NumericVector b(p),d(p),fit(n);
  NumericVector e=y-mu,e1(n),e2(n);
  double ve=var(y);
  // Fitting loop
  for(int i=0; i<iter; i++){
    C = -0.5/ve;
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // OLS marker effect
      b1 = (sum(X(_,j)*e)+xx[j]*b0)/xx[j];
      // EN marker effect
      if(b1>0){
        // Regularization of positive OLS
        b2 = (sum(X(_,j)*e)+xx[j]*b0-L1)/(xx[j]+L2);
        if(b2<0){b2=0;};
      }else{
        // Regularization of negative OLS
        b2 = (sum(X(_,j)*e)+xx[j]*b0+L1)/(xx[j]+L2);
        if(b2>0){b2=0;};
      }
      e1 = e-X(_,j)*(b1-b0);
      e2 = e-X(_,j)*(b2-b0); 
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj); d[j] = pj;
      // Smple from Bernoulli
      if(pj>0.5){
        b[j] = b1;
      }else{
        b[j] = b2;
      }
      // Update residuals
      e = e - X(_,j)*(b[j]-b0);
    }
    // Update intercept and variance
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    ve = sum(e*y)/(n-1);
    L1 = 0.5*MSx*(sqrt(ve)*alpha);
  }
  // Get fitted values
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*b)+mu;}
  double h2 = 1-ve/var(y);
  // Return output
  return List::create(Named("mu") = mu, Named("b") = b,
                      Named("PrQTL") = d, Named("h2") = h2,
                      Named("Lmb1") = L1, Named("Lmb2") = L2,
                      Named("MSx") = MSx, Named("hat") = fit);}
