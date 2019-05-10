#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP ElsNet(NumericVector y, NumericMatrix X){
  int maxit = 500;
  double tol = 10e-12;
  // Functions starts here
  int p = X.ncol();
  int n = X.nrow();
  // Beta and mu
  double eM, mu = mean(y), bL1, bL2;
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int i=0; i<p; i++){xx[i] = sum(X(_,i)*X(_,i));}
  double Lmb2=mean(xx)/p, Lmb1=Lmb2/(p*p);
  // Convergence control
  NumericVector bc(p), yx(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      e = e+X(_,j)*b[j];
      yx[j] = sum(e*X(_,j));
      bL2 = (yx[j])/(xx[j]+Lmb2);
      if(yx[j]>0){ bL1 = (yx[j]-Lmb1)/(xx[j]); if(bL1<0){bL1=0;}
            }else{ bL1 = (yx[j]+Lmb1)/(xx[j]); if(bL1>0){bL1=0;} }
      b[j] = (bL1+bL2)/2;
      e = e-X(_,j)*b[j];}
    // Update regularization and intercept
    Lmb1 = mean(abs(yx)-b*xx);
    Lmb2 = 0.5*mean((yx-b*xx)/b);
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit = y-e;
  double h2 = 1-(sum(e*y)/(n-1))/var(y);
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Lmb1")=Lmb1, Named("Lmb2")=Lmb2);}