#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP BSR(NumericVector y, NumericMatrix gen){
  int maxit = 250;
  double tol = 10e-8;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0, eM, ve, mu = mean(y);
  NumericVector b(p), vb(p), e=y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  double cxx = mean(xx), mxx = 4*cxx;
  NumericVector Lmb = vb+cxx;
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
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e = e-gen(_,j)*(b[j]-b0);}
    // Variance components update
    vb = b*b+ve/(xx+Lmb);
    ve = sum(e*y)/(n-1);
    Lmb = ve/vb;
    Lmb = ifelse(Lmb>mxx,mxx,Lmb);
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  bool CONV;
  if(numit<maxit){CONV = true;}else{CONV = false;}
  // Fitting the model
  NumericVector fit(n);
  for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  double h2 = 1-ve/var(y);
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("h2")=h2,
                      Named("Vb")=vb,
                      Named("Ve")=ve,
                      Named("conv")=CONV);}
