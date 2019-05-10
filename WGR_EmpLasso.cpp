#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP Lasso(NumericVector y, NumericMatrix X){
  int maxit = 300;
  double tol = 10e-12;
  // Functions starts here
  int p = X.ncol();
  int n = X.nrow();
  // Beta, mu and epsilon
  double eM, mu = mean(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int i=0; i<p; i++){xx[i] = sum(X(_,i)*X(_,i));}
  double Lmb=mean(xx)/p;
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
      if(yx[j]>0){
        b[j] = (yx[j]-Lmb)/xx[j];
        if(b[j]<0){b[j]=0;}
      }else{
        b[j] = (yx[j]+Lmb)/xx[j];
        if(b[j]>0){b[j]=0;}
      }
      e = e-X(_,j)*b[j];}
    // Update regularization and intercept
    Lmb = 2*mean(abs(yx)-b*xx);
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
                      Named("Lmb")=Lmb);}
