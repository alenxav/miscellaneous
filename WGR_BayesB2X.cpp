#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP BayesB2(NumericVector y, NumericMatrix X1, NumericMatrix X2,
             double it = 1500, double bi = 500,
             double pi = 0.95, double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int n = X1.nrow();
  int p1 = X1.ncol();
  int p2 = X2.ncol();
  // Estimate crossproducts and MSx
  NumericVector xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = sum(X1(_,i)*X1(_,i));
    vx1[i] = var(X1(_,i));}
  double MSx1 = sum(vx1);
  NumericVector xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = sum(X2(_,i)*X2(_,i));
    vx2[i] = var(X2(_,i));}
  double MSx2 = sum(vx2);
  // Get priors
  double vy = var(y);
  double Sb1 = (R2)*df*vy/MSx1;
  double Sb2 = (R2)*df*vy/MSx2;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b_t0,b_t1,b_t2,eM,h2,C,MU,VE,cj,dj,pj,vg,ve=vy;
  NumericVector d1(p1),b1(p1),D1(p1),B1(p1),VB1(p1),fit(n);
  NumericVector d2(p2),b2(p2),D2(p2),B2(p2),VB2(p2);
  NumericVector vb1=b1+Sb1,vb2=b2+Sb2,Lmb1=ve/vb1,Lmb2=ve/vb2,e=y-mu,e1(n),e2(n);
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    // Update marker effects 1
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      // Sample marker effect
      b_t1 = R::rnorm((sum(X1(_,j)*e)+xx1[j]*b_t0)/(xx1[j]+Lmb1[j]),sqrt(ve/(xx1[j]+Lmb1[j])));
      b_t2 = R::rnorm(0,sqrt(ve/(xx1[j]+Lmb1[j])));
      e1 = e-X1(_,j)*(b_t1-b_t0); // Pr(with marker)
      e2 = e-X1(_,j)*(b_t2-b_t0); // Pr(without marker)
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b1[j] = b_t1; d1[j] = 1;
      }else{
        b1[j] = b_t2; d1[j] = 0;
      }
      // Update marker variance and residuals
      vb1[j] = (Sb1+b1[j]*b1[j])/R::rchisq(df+1);
      e = e - X1(_,j)*(b1[j]-b_t0);
    }
    // Update marker effects 1
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      // Sample marker effect
      b_t1 = R::rnorm((sum(X2(_,j)*e)+xx2[j]*b_t0)/(xx2[j]+Lmb2[j]),sqrt(ve/(xx2[j]+Lmb2[j])));
      b_t2 = R::rnorm(0,sqrt(ve/(xx2[j]+Lmb2[j])));
      e1 = e-X2(_,j)*(b_t1-b_t0); // Pr(with marker)
      e2 = e-X2(_,j)*(b_t2-b_t0); // Pr(without marker)
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b2[j] = b_t1; d2[j] = 1;
      }else{
        b2[j] = b_t2; d2[j] = 0;
      }
      // Update marker variance and residuals
      vb2[j] = (Sb2+b2[j]*b2[j])/R::rchisq(df+1);
      e = e - X2(_,j)*(b2[j]-b_t0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    // Store posterior sums
    if(i>bi){
      MU=MU+mu; VE=VE+ve;
      B1=B1+b1; D1=D1+d1; VB1=VB1+vb1; 
      B2=B2+b2; D2=D2+d2; VB2=VB2+vb2; 
    }
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; VE = VE/MCMC;
  B1 = B1/MCMC; D1 = D1/MCMC; VB1 = VB1/MCMC; 
  B2 = B2/MCMC; D2 = D2/MCMC; VB2 = VB2/MCMC; 
  // Get fitted values and h2
  vg = sum(VB1)+sum(VB2); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X1(k,_)*B1)+sum(X2(k,_)*B2)+MU;}
  // Return output
  return List::create(Named("mu") = MU,
                      Named("b1") = B1, Named("d1") = D1, Named("vb1") = VB1,
                      Named("b2") = B2, Named("d2") = D2, Named("vb2") = VB2,
                      Named("ve") = VE, Named("hat") = fit, Named("h2") = h2);}
