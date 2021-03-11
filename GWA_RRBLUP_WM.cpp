#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emREML(NumericVector Y, NumericMatrix X,
            Rcpp::Nullable<Rcpp::NumericVector> D = R_NilValue,
            int maxit = 500, double tol = 10e-8){
  // Functions starts here
  int p = X.ncol();
  int n0 = X.nrow();
  // Handle missings y
  NumericVector fit(n0);
  NumericVector o(n0), y(n0);
  o = ifelse(is_na(Y),0,1);
  y = ifelse(is_na(Y),0,Y);
  int n = sum(o);
  // Weights
  bool P_WEIGHTS = FALSE;
  NumericVector d(p);
  if (D.isNotNull()){P_WEIGHTS = TRUE; d=D;}
  // Beta, mu and epsilon
  double b0, eM, ve=var(y), vb=ve/p, h2, mu=sum(y)/n;
  double vy = sum((y-mu)*(y-mu))/(n-1);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i)*o);
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
        b[j] = (sum(X(_,j)*e*o)+xx[j]*b0)/(xx[j]+Lmb/d[j]);
      }else{
        b[j] = (sum(X(_,j)*e*o)+xx[j]*b0)/(xx[j]+Lmb);}
      e = (e-X(_,j)*(b[j]-b0))*o;}
    // Variance components update
    ve = sum(e*y*o)/(n-1);
    vb = (vy-ve)/MSx;
    Lmb = ve/vb;
    // Intercept update
    eM = sum(e*o)/n;
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  for(int k=0; k<n0; k++){ fit[k] = sum(X(k,_)*b)+mu; }
  e = (y-fit)*o;
  h2 = vb*MSx/(vb*MSx+ve);
  // Genome-wide screening
  NumericVector LRT(p),y0(n0),e0(n0),e1(n0),b_ols(p);
  double ve0, ve1, L0, L1;
  for(int j=0; j<p; j++){
    // Full conditional phenotype
    y0 = e+X(_,j)*b[j];
    // Fixed effect marker
    b_ols[j] = sum(X(_,j)*y0*o)/xx[j];
    // Null model
    e0 = y0-mean(y0);
    // Alternative model
    e1 = y0-X(_,j)*b_ols[j];
    e1 = e1-mean(e1);
    // ReML variance
    ve0 = sum(y0*e0)/(n-1);
    ve1 = sum(y0*e1)/(n-1);
    // Likelihood ratio
    L0 = -sum(e0*e0)/(2*ve0)-0.5*n*log(6.28*ve0);
    L1 = -sum(e1*e1)/(2*ve1)-0.5*n*log(6.28*ve1);
    LRT[j] = 2*(L1-L0);}
  // Output
  return List::create(Named("mu")=mu, Named("b")=b, Named("hat")=fit,
                      Named("h2")=h2, Named("LRT")=LRT);}
