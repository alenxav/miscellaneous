#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP BDpi(NumericVector y, NumericMatrix X,
              double it = 2500, double bi = 500,
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
  double priorA = 1;
  double priorB = 1;
  double pi = 0.5;
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,b2,eM,h2,C,MU,VE,Pi,cj,dj,pj,vg,ve=vy;
  double PiAlpha,PiBeta,PiMean,PiVar;
  NumericVector d(p),b(p),D(p),B(p),VB(p),fit(n);
  NumericVector vb=b+Sb,Lmb=ve/vb,e=y-mu,e1(n),e2(n);
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]),sqrt(ve/(xx[j]+Lmb[j])));
      b2 = R::rnorm(0,sqrt(ve/(xx[j]+Lmb[j])));        
      e1 = e-X(_,j)*(b1-b0); // Pr(with marker)
      e2 = e-X(_,j)*(b2-b0); // Pr(without marker)
      // Pr(marker included)
      cj = exp(C*sum(e1*e1)); // Likelihood(with marker)
      dj = exp(C*sum(e2*e2)); // Likelihood(without marker)
      pj = (1-pi)*cj/dj;
      if(pj>1) pj = 1;
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = b2; d[j] = 0;
      }
      // Update marker variance and residuals
      vb[j] = (Sb+b[j]*b[j])/R::rchisq(df+1);
      e = e - X(_,j)*(b[j]-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    // Update Pi from beta
    PiMean = mean(d); PiVar = var(d);
    PiAlpha = priorA+((1-PiMean)/PiVar-1/PiMean)*(PiMean*PiMean);
    PiBeta = priorB+PiAlpha*(1/PiMean-1);
    pi = R::rbeta(PiAlpha,PiBeta);
    // Store posterior sums
    if(i>bi){
      MU=MU+mu; B=B+b; D=D+d;
      VB=VB+vb; VE=VE+ve; Pi = Pi+pi;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC; D = D/MCMC;
  VB = VB/MCMC; VE = VE/MCMC; Pi = Pi/MCMC;
  // Getting GWAS results
  NumericVector PVAL = -log(1-D);
  // Get fitted values and h2
  vg = sum(VB); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("pi") = Pi, 
                      Named("hat") = fit, Named("PVAL") = PVAL,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);}