#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP BayesCpi(NumericVector y, NumericMatrix X,
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
  double pi = 0.5;
  double vy = var(y);
  double Sb = df*(R2)*vy/MSx/(1-pi);
  double Se = df*(1-R2)*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,eM,h2,C,MU,VB,VE,Pi,cj,dj,pj,vg,ve=vy,vb=Sb;
  double PiAlpha,PiBeta,PiMean,PiVar;
  NumericVector d(p),b(p),D(p),B(p),fit(n);
  NumericVector e=y-mu,e1(n),e2(n);
  double Lmb=ve/vb;
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb),
                    sqrt(ve/(xx[j]+Lmb)));
      e1 = e-X(_,j)*(b1-b0); // Pr(with marker)
      e2 = e-X(_,j)*(0-b0); // Pr(without marker)
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = R::rnorm(0,sqrt(ve/(xx[j]+Lmb))); d[j] = 0;
      }
      // Update residuals
      e = e - X(_,j)*(b[j]-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update variance components and lambda
    vb = (sum(b*b)+Sb)/R::rchisq(df+p);
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    // Update Pi from beta
    PiMean = mean(d); PiVar = var(d);
    PiAlpha = 4+((1-PiMean)/PiVar-1/PiMean)*(PiMean*PiMean);
    PiBeta = 1+PiAlpha*(1/PiMean-1);
    pi = R::rbeta(PiAlpha,PiBeta);
    Sb = df*(R2)*vy/MSx/(1-pi);
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; D=D+d; VB=VB+vb; VE=VE+ve; Pi = Pi+pi;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC; D = D/MCMC;
  VB = VB/MCMC; VE = VE/MCMC; Pi = Pi/MCMC;
  // Get fitted values and h2
  vg = VB*MSx; h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("pi") = Pi,
                      Named("hat") = fit, Named("h2") = h2
                      Named("vb") = VB, Named("ve") = VE);}