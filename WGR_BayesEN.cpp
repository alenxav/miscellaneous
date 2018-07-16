#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP BayesEN(NumericVector y, NumericMatrix X,
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
  // Get priors
  double alpha = 0.05;
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx/(1-alpha);
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  double L1 = 0.5*MSx*(sd(y)*alpha);
  double L2 = vy/Sb;
  // Create empty objects
  double b0,b1,b2,eM,h2,MU,VE,VB,vg,ve=vy,vb=Sb;
  NumericVector b(p),B(p),p0(p),fit(n),e=y-mu;
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // RHS
      b1 = sum(X(_,j)*e)+xx[j]*b0;
      if(b1>0){
        // Regularization of positive OLS
        b2 = (b1-L1)/(xx[j]+L2);
        if(b2<0){b2=0;p0[j]=0;}else{p0[j]=1;}
      }else{
        // Regularization of negative OLS
        b2 = (b1+L1)/(xx[j]+L2);
        if(b2>0){b2=0;p0[j]=0;}else{p0[j]=1;}
      }
      b[j] = R::rnorm(b2,sqrt(ve/(xx[j]+L2)));
      e = e-X(_,j)*(b[j]-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update variance components
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    vb = (sum(b*b)+Sb)/R::rchisq(sum(p0)+df);
    // Update regularization coefficients
    alpha = R::rbeta(.2,2);
    L1 = 0.5*(ve/vb)*sqrt(ve)*alpha;
    L2 = (ve/vb)/(1-alpha);
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; VB=VB+vb;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC;
  B = B/MCMC;
  VB = VB/MCMC;
  // Get fitted values and h2
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  VE = sum((y-fit)*(y-fit))/n;
  vg = VB*MSx; h2 = vg/(vg+VE);
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit, Named("h2") = h2);}
