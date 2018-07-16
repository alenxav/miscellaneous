#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP RMT(NumericVector y, NumericMatrix gen, int maxit = 1000, double tol = 10e-5){
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0, eM, ve, mu = mean(y);
  NumericVector b(p), vb(p), e=y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  double cxx = mean(xx);
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
  // Genome-wide screening
  NumericVector LRT(p),PVAL(p),y0(n),e0(n),e1(n),b_ols(p);
  double ve0, ve1, L0, L1;
  for(int j=0; j<p; j++){
    // Full conditional phenotype
    y0 = e+gen(_,j)*b[j];
    // Fixed effect marker
    b_ols[j] = sum(gen(_,j)*y0)/xx[j];
    // Null model
    e0 = y0-mean(y0);
    // Alternative model
    e1 = y0-gen(_,j)*b_ols[j];
    e1 = e1-mean(e1);
    // ReML variance
    ve0 = sum(y0*e0)/(n-1);
    ve1 = sum(y0*e1)/(n-1);
    // Likelihood ratio
    L0 = -sum(e0*e0)/(2*ve0)-0.5*n*log(6.28*ve0);
    L1 = -sum(e1*e1)/(2*ve1)-0.5*n*log(6.28*ve1);
    LRT[j] = 2*(L1-L0);}
  PVAL = -log10(1-pchisq(LRT,1,true));
  // Fitting the model
  NumericVector fit(n);
  for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  double h2 = sum(vb)/(sum(vb)+ve);
  // Output
  return List::create(Named("LRT")=LRT,
                      Named("PVAL")=PVAL,
                      Named("h2")=h2,
                      Named("mu")=mu,
                      Named("b")=b,
                      Named("b_ols")=b_ols,
                      Named("hat")=fit,
                      Named("Vb")=vb,
                      Named("Ve")=ve,
                      Named("conv")=CONV);}
