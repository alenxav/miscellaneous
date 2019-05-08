#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP gREML2(NumericVector y, NumericMatrix X){
  int maxit = 500;
  int VCup = 20;
  double tol = 10e-8;
  // Functions starts here
  int p = X.ncol();
  int n = X.nrow();
  // Beta, mu and epsilon
  double b0, eM, ve, vb, h2, mu = mean(y), vy = var(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), ze(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));}
  double Lmb = mean(xx);
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      b0 = b[j];
      ze[j] = sum(X(_,j)*e);
      b[j] = (ze[j]+xx[j]*b0)/(xx[j]+Lmb);
      e = e-X(_,j)*(b[j]-b0);
    }
    // Variance components update
    if(numit%VCup==0){
      ve = sum(e*e)/n;
      //ve = sum(e*y)/(n-1); 
      vb = mean(b/ze)*ve;
      Lmb = ve/vb;
    }
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit(n);
  fit = y-e;
  h2 = 1-ve/vy;
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Vb")=vb, Named("Ve")=ve);}
