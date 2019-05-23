#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP mkr(NumericMatrix Y, NumericMatrix K){
  // Convergence parameters
  int maxit = 200; double tol = 10e-8;
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  Rcpp::Function eigen = base["eigen"];
  // Eigendecomposition
  List EIG = eigen(K,true);
  NumericMatrix X = EIG[1];
  NumericVector d = EIG[0];
  // Functions starts here
  int k = Y.ncol(), p = X.ncol(), n0 = X.nrow();
  NumericMatrix fit(n0,k),o(n0,k),y(n0,k),e(n0,k);
  for(int i=0; i<k; i++){
    o(_,i) = ifelse(is_na(Y(_,i)),0,1);
    y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));}
  NumericVector n = colSums(o);
  // Mu
  NumericVector mu = colSums(y)/n;
  for(int j=0; j<k; j++){y(_,j) = (y(_,j)-mu(j))*o(_,j);}
  // Marker variance
  NumericMatrix xx(p,k);
  double MSx = mean(d);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){
      xx(i,j) = sum(X(_,i)*X(_,i)*o(_,j)); }}
  // Beta, intersept and residuals
  NumericMatrix b(p,k),vb(k,k),iG(k,k),rho(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),vy(k),ve(k),RHS(k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx;
    rho(i,i) = 1;}
  iG = solve(vb);
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
        LHS(i,i) = iG(i,i)+(xx(j,i)/(ve(i)/d(j)) );
        RHS(i) = (sum(e(_,i)*X(_,j))+xx(j,i)*b0(i))/(ve(i)/d(j));}
      // Update effects
      b1 = solve(LHS, RHS);
      b(j,_) = b1;
      // Update residuals
      for(int i=0; i<k; i++){
        e(_,i) = (e(_,i)-X(_,j)*(b1(i)-b0(i)))*o(_,i);}
    }
    // Residual variance components update
    for(int i=0; i<k; i++){ ve(i) = (sum(e(_,i)*y(_,i)))/n(i);}
    // Genetic covariance components update
    for(int i=0; i<n0; i++){ for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j));}}
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
      // Diag VC
      if(i==j){
        vb(i,j) = (1.01*vy(i)-ve(i))/MSx;
      }else{
        vb(i,j) = (sum(fit(_,i)*y(_,j))+sum(fit(_,j)*y(_,i))) / ((n(i)+n(j))*MSx);
      }}}
    iG = solve(vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector h2(k); 
  for(int i=0; i<n0; i++){for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j))+mu(j);}}
  h2 = 1-ve/vy;
  // Genetic correlations
  NumericMatrix GC(k,k);
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("hat")=fit, Named("h2")=h2,
                      Named("Vb")=vb, Named("Ve")=ve,
                      Named("GC")=GC);}