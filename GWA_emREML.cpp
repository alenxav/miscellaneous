#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emREML(NumericVector y, NumericMatrix gen, int maxit = 500, double tol = 10e-6){
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0, eM, ve, vb, h2, mu = mean(y), vy = var(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));}
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
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e = e-gen(_,j)*(b[j]-b0);}
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
  for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  h2 = vb*MSx/(vb*MSx+ve);
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
  // Output
  return List::create(Named("mu")=mu, Named("b")=b, Named("b_LS")=b_ols,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Vb")=vb, Named("Ve")=ve,
                      Named("LRT")=LRT, Named("PVAL")=PVAL);}