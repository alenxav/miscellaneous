#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using arma::pinv;

// [[Rcpp::export]]
SEXP MV(NumericMatrix Y, NumericMatrix X){
  // Convergence criteria
  int maxit = 250;
  double tol = 10e-8;
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  // Functions starts here
  int k = Y.ncol(), p = X.ncol(), n0 = X.nrow();
  NumericMatrix fit(n0,k),o(n0,k),y(n0,k),eps(n0,k),e(n0,k);
  for(int i=0; i<k; i++){
    o(_,i) = ifelse(is_na(Y(_,i)),0,1);
    y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));}
  NumericVector n = colSums(o);
  // Starting point for overall mean
  NumericVector mu = colSums(y)/n, mu0(k);
  for(int j=0; j<k; j++){eps(_,j) = (y(_,j)-mu(j))*o(_,j);}
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
  NumericMatrix b(p,k),vb(k,k),rho(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),vy(k),ve(k),RHS(k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx(i);
    rho(i,i) = 1;}
  // Inverse G
  arma::mat G = as<arma::mat>(vb);
  arma::mat invG = pinv(G);
  NumericMatrix iG = wrap(invG);
  // Convergence control
  NumericMatrix bc(p,k);
  int numit = 0;
  double cnv = 1;
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
    // Update mu and epsilon
    mu0 = colSums(e)/n; 
    mu = mu+mu0;
    for(int j=0; j<k; j++){eps(_,j) = (y(_,j)-mu(j))*o(_,j);}
    // Residual variance components update
    for(int i=0; i<k; i++){ ve(i) = (sum(e(_,i)*y(_,i)))/n(i);}
    // Genetic covariance components update
    for(int i=0; i<n0; i++){
      for(int j=0; j<k; j++){fit(i,j)=sum(X(i,_)*b(_,j));}}
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
      tmp = (sum(fit(_,i)*eps(_,j))+sum(fit(_,j)*eps(_,i)));
      vb(i,j) = tmp/((n(i)*MSx(i))+(n(j)*MSx(j)));}}
    // Generalized inverse of G
    G = as<arma::mat>(vb);
    invG = pinv(G,0.01);
    iG = wrap(invG);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){fit(i,j)=sum(X(i,_)*b(_,j))+mu(j);}}
  // Heritability and genetic correlations
  NumericVector h2 = 1-ve/vy;
  NumericMatrix GC(k,k);
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
    GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  // Output list
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("h2")=h2,
                      Named("Vb")=vb,
                      Named("Ve")=ve,
                      Named("GC")=GC);}
