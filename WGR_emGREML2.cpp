#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emML2(NumericVector y, NumericMatrix X1, NumericMatrix X2,
           Rcpp::Nullable<Rcpp::NumericVector> D1 = R_NilValue,
           Rcpp::Nullable<Rcpp::NumericVector> D2 = R_NilValue,
           int maxit = 350, double tol = 10e-8){
  // Functions starts here
  int p1 = X1.ncol();
  int p2 = X2.ncol();
  int n = X1.nrow();
  // Weights
  bool P1_WEIGHTS = FALSE;
  bool P2_WEIGHTS = FALSE;
  NumericVector d1(p1), d2(p2);
  if (D1.isNotNull()){P1_WEIGHTS = TRUE; d1=D1;}
  if (D2.isNotNull()){P2_WEIGHTS = TRUE; d2=D2;}
  // Beta, mu and epsilon
  double b0, eM, ve, vb1, vb2, h2, mu = mean(y);
  NumericVector b1(p1), b2(p2), u1(n), u2(n), cY(n), e = y-mu;
  // Marker variance
  NumericVector x1x1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    x1x1[i] = sum(X1(_,i)*X1(_,i));
     vx1[i] = var(X1(_,i));}
  double MSx1 = sum(vx1), Lmb1=MSx1;
  NumericVector x2x2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    x2x2[i] = sum(X2(_,i)*X2(_,i));
     vx2[i] = var(X2(_,i));}
  double MSx2 = sum(vx2), Lmb2=MSx2;
  // Convergence control
  NumericVector bc1(p1), bc2(p2);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Save b(t0)
    bc1 = b1+0;
    bc2 = b2+0;
    // Regression coefficients loop 1
    for(int j=0; j<p1; j++){
      b0 = b1[j];
      if(P1_WEIGHTS){
        b1[j] = (sum(X1(_,j)*e)+x1x1[j]*b0)/(x1x1[j]+Lmb1/d1[j]);
      }else{
        b1[j] = (sum(X1(_,j)*e)+x1x1[j]*b0)/(x1x1[j]+Lmb1);}
      e = e-X1(_,j)*(b1[j]-b0);}
    // Regression coefficients loop 2
    for(int j=0; j<p2; j++){
      b0 = b2[j];
      if(P2_WEIGHTS){
        b2[j] = (sum(X2(_,j)*e)+x2x2[j]*b0)/(x2x2[j]+Lmb2/d2[j]);
      }else{
        b2[j] = (sum(X2(_,j)*e)+x2x2[j]*b0)/(x2x2[j]+Lmb2);}
      e = e-X2(_,j)*(b2[j]-b0);}
    // Fitting the model
    for(int k=0; k<n; k++){
      u1[k] = sum(X1(k,_)*b1);
      u2[k] = sum(X2(k,_)*b2);
    }
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components update
    cY = u1+u2+e;
    ve = sum(e*cY)/n;
    vb1 = (sum(u1*cY)/n)/MSx1;
    vb2 = (sum(u2*cY)/n)/MSx2;
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    // Convergence
    ++numit;
    cnv = sum(abs(bc1-b1))+sum(abs(bc2-b2));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit = mu+u1+u2;
  h2 = 1-ve/var(y);
  // Output
  return List::create(Named("mu")=mu,
                      Named("b1")=b1, Named("b2")=b2, 
                      Named("Vb1")=vb1, Named("Vb2")=vb2, Named("Ve")=ve,
                      Named("u1")=u1, Named("u2")=u2,
                      Named("MSx1")=MSx1, Named("MSx2")=MSx2,
                      Named("h2")=h2, Named("hat")=fit);}