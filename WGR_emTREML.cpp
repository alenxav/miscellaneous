#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP tREML(NumericVector y, NumericMatrix X,
           int maxit = 500, double tol = 10e-8){
  // Functions starts here
  int p = X.ncol();
  int n = X.nrow();
  // Coefficients and marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx=sum(vx),b0,eM,Ve,b1,Vy=var(y),mu=mean(y),Va,h2;
  NumericVector Vb(p),b(p),bc(p),d(p),fit(n),e=y-mu,Lmb=p+MSx;
  // Convergence control
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      b1 = (sum(X(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j));
      b[j] = b1;
      // Residuals update
      e = e-X(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Vb = b*b+(Ve/(xx+Lmb));
    Lmb = sqrt(MSx*Ve/Vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  for(int k=0; k<n; k++){ fit[k] = sum(X(k,_)*b)+mu; }
  // Output
  Va=(Vy-Ve)/MSx;
  d=Ve/(Va*Lmb);
  h2=Va*MSx/(Va*MSx+Ve);
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("hat")=fit, Named("h2")=h2, Named("D")=d,
                      Named("Vb")=Va, Named("Ve")=Ve);}
