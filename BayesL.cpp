#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP BayesL(NumericVector y, NumericMatrix X,
            double it = 1500, double bi = 500,
            double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  double Phi = MSx*(1-R2)/R2;
  // Get priors
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,eM,h2,MU,VE,vg,ve=vy;
  NumericVector b(p),D(p),B(p),VB(p),fit(n);
  NumericVector vb=b+Sb,Lmb=ve/vb,e=y-mu;
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]),
                    sqrt(ve/(xx[j]+Lmb[j])));
      e = e-X(_,j)*(b1-b0);
      b[j] = b1;
      // Update marker variance and residuals
      vb[j] = (Sb+b1*b1)/R::rchisq(df+1);
      e = e - X(_,j)*(b1-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = sqrt(Phi*ve/vb);
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; VB=VB+vb; VE=VE+ve;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC;
  VB = VB/MCMC; VE = VE/MCMC;
  // Get fitted values and h2
  vg = sum(VB); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);}