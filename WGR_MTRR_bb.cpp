#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP MVRR(NumericMatrix Y, NumericMatrix X){
  int maxit = 200; double tol = 10e-8;
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  // Functions starts here
  int k = Y.ncol(), p = X.ncol(), n0 = X.nrow();
  // Handle missings Y's
  NumericMatrix fit(n0,k),o(n0,k),y(n0,k),e(n0,k);
  for(int i=0; i<k; i++){
    o(_,i) = ifelse(is_na(Y(_,i)),0,1);
    y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));}
  NumericVector n = colSums(o);
  // Marker variance
  NumericMatrix xx(p,k), vx(p,k);
  double tmp;
  for(int i=0; i<p; i++){
    for(int j=0; j<k; j++){
     xx(i,j) = sum(X(_,i)*X(_,i)*o(_,j));
      tmp = sum(X(_,i)*o(_,j))/n(j);
      vx(i,j) = xx(i,j)/n(j)-tmp*tmp;}}
  NumericVector MSx = colSums(vx);
  // Beta, intersept and residuals
  NumericMatrix b(p,k),vb(k,k),iG(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),eM(k),mu(k),vy(k),ve(k),RHS(k);
  mu = colSums(y)/n;
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = (y(_,i)-mu(i))*o(_,i);
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx(i);}
  iG = solve(vb);
  // Convergence control
  NumericMatrix bc(p,k);
  int numit = 0; double cnv = 1;
  // Loop
  while(numit<maxit){
    // Gauss-Seidel loop
    bc = b+0;
    for(int j=0; j<p; j++){
      b0 = b(j,_);
      LHS = iG+0;
      for(int i=0; i<k; i++){
        LHS(i,i) = iG(i,i)+(xx(j,i)/ve(i));
        RHS(i) = (sum(e(_,i)*X(_,j))+xx(j,i)*b0(i))/ve(i);}
      // Update effects
      b1 = solve(LHS, RHS);
      b(j,_) = b1;
      // Update residuals
      for(int i=0; i<k; i++){
        e(_,i) = (e(_,i)-X(_,j)*(b1(i)-b0(i)))*o(_,i);}}
    // Intercept update
    eM = colSums(e)/n;
    mu = mu+eM;
    for(int j=0; j<k; j++){e(_,j) = (e(_,j)-eM(j))*o(_,j);}
    // Variance components update
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = sum(b(_,i)*b(_,j))/(p-1);}}
    for(int i=0; i<k; i++){ve(i) = sum(e(_,i)*y(_,i))/(n(i)-1);vb(i,i) = (1.001*vy(i)-ve(i))/MSx(i);}
    iG = solve(vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model and heritability
  NumericVector h2(k);
  for(int i=0; i<n0; i++){for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j))+mu(j);}}
  for(int i=0; i<k; i++){ h2 = (vb(i,i)*MSx(i))/((vb(i,i)*MSx(i))+ve); }
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("hat")=fit, Named("h2")=h2,
                      Named("Vb")=vb, Named("Ve")=ve);}