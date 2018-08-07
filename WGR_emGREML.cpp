#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emREML(NumericVector y, NumericMatrix X,
           Rcpp::Nullable<Rcpp::NumericVector> D = R_NilValue,
           int maxit = 500, double tol = 10e-8){
  // Functions starts here
  int p = X.ncol();
  int n = X.nrow();
  // Weights
  bool P_WEIGHTS = FALSE;
  NumericVector d(p);
  if (D.isNotNull()){P_WEIGHTS = TRUE; d=D;}
  // Beta, mu and epsilon
  double b0, eM, ve, vb, h2, mu = mean(y), vy = var(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx), Lmb=MSx;
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
      if(P_WEIGHTS){
        b[j] = (sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb/d[j]);
          }else{
        b[j] = (sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);}
      e = e-X(_,j)*(b[j]-b0);}
    // Variance components update
    ve = sum(e*y)/(n-1);
    vb = (vy-ve)/MSx;
    Lmb = ve/vb;
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
  for(int k=0; k<n; k++){ fit[k] = sum(X(k,_)*b)+mu; }
  h2 = vb*MSx/(vb*MSx+ve);
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Vb")=vb, Named("Ve")=ve);}
