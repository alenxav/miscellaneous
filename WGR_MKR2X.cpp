#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP mkr2X(NumericMatrix Y, NumericMatrix K1, NumericMatrix K2){
  // Convergence parameters
  int maxit = 200; double tol = 10e-8;
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  Rcpp::Function eigen = base["eigen"];
  // Eigendecomposition
  List EIG1 = eigen(K1,true);
  NumericMatrix X1 = EIG1[1];
  NumericVector d1 = EIG1[0];
  List EIG2 = eigen(K2,true);
  NumericMatrix X2 = EIG2[1];
  NumericVector d2 = EIG2[0];
  // Functions starts here
  int k = Y.ncol(), p1 = X1.ncol(), p2 = X2.ncol(), n0 = X1.nrow();
  NumericMatrix fit(n0,k),g1(n0,k),g2(n0,k),o(n0,k),y(n0,k),e(n0,k);
  for(int i=0; i<k; i++){
    o(_,i) = ifelse(is_na(Y(_,i)),0,1);
    y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));}
  NumericVector n = colSums(o);
  // Mu
  NumericVector mu = colSums(y)/n;
  for(int j=0; j<k; j++){y(_,j) = (y(_,j)-mu(j))*o(_,j);}
  // Marker variance
  NumericMatrix xx1(p1,k), xx2(p2,k);
  for(int i=0; i<p1; i++){
    for(int j=0; j<k; j++){
      xx1(i,j) = sum(X1(_,i)*X1(_,i)*o(_,j));}}
  for(int i=0; i<p2; i++){
    for(int j=0; j<k; j++){
      xx2(i,j) = sum(X2(_,i)*X2(_,i)*o(_,j));}}
  double MSx1 = mean(d1), MSx2 = mean(d2);
  // Beta, intersept and residuals
  NumericMatrix b1(p1,k),vb1(k,k),iG1(k,k),LHS1(k,k);
  NumericMatrix b2(p2,k),vb2(k,k),iG2(k,k),LHS2(k,k);
  NumericVector b_0(k),b_1(k),vy(k),ve(k),RHS1(k),RHS2(k);
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      vb2(i,j) = 0; vb1(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb1(i,i) = ve(i)/MSx1;
    vb2(i,i) = ve(i)/MSx2;
  }
  iG1 = solve(vb1);
  iG2 = solve(vb2);
  // Convergence control
  NumericMatrix bc1(p1,k), bc2(p2,k);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Convergence parameter (coef at time zero)
    bc1 = b1+0; bc2 = b2+0;
    // Gauss-Seidel loop 1
    for(int j=0; j<p1; j++){
      b_0 = b1(j,_);
      LHS1 = iG1+0;
      for(int i=0; i<k; i++){
        LHS1(i,i) = iG1(i,i)+(xx1(j,i)/(ve(i)/d1(j)));
        RHS1(i) = (sum(e(_,i)*X1(_,j))+xx1(j,i)*b_0(i))/(ve(i)/d1(j));}
      b_1 = solve(LHS1, RHS1);
      b1(j,_) = b_1;
      for(int i=0; i<k; i++){
        e(_,i) = (e(_,i)-X1(_,j)*(b_1(i)-b_0(i)))*o(_,i);}
    }
    // Gauss-Seidel loop 2
    for(int j=0; j<p2; j++){
      b_0 = b2(j,_);
      LHS2 = iG2+0;
      for(int i=0; i<k; i++){
        LHS2(i,i) = iG2(i,i)+(xx2(j,i)/(ve(i)/d2(j)));
        RHS2(i) = (sum(e(_,i)*X2(_,j))+xx2(j,i)*b_0(i))/(ve(i)/d2(j));}
      b_1 = solve(LHS2, RHS2);
      b2(j,_) = b_1;
      for(int i=0; i<k; i++){
        e(_,i) = (e(_,i)-X2(_,j)*(b_1(i)-b_0(i)))*o(_,i);}
    }
    // Residual variance components update
    for(int i=0; i<k; i++){ ve(i) = (sum(e(_,i)*y(_,i)))/n(i);}
    // Genetic covariance components update
    for(int i=0; i<n0; i++){
      for(int j=0; j<k; j++){
        g1(i,j) = sum(X1(i,_)*b1(_,j));
        g2(i,j) = sum(X2(i,_)*b2(_,j));  
      }}
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        vb1(i,j) = (sum(g1(_,i)*y(_,j))+sum(g1(_,j)*y(_,i))) / ((n(i)+n(j))*MSx1);
        vb2(i,j) = (sum(g2(_,i)*y(_,j))+sum(g2(_,j)*y(_,i))) / ((n(i)+n(j))*MSx2);
      }}
    iG1 = solve(vb1);
    iG2 = solve(vb2);
    // Convergence
    ++numit;
    cnv = sum(abs(bc1-b1))+sum(abs(bc2-b2));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector h2(k); 
  for(int j=0; j<k; j++){ fit(_,j) = mu(j)+g1(_,j)+g2(_,j); }
  h2 = 1-ve/vy;
  // Genetic correlations
  NumericMatrix GC1(k,k),GC2(k,k);
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
    GC1(i,j)=vb1(i,j)/(sqrt(vb1(i,i)*vb1(j,j)));
    GC2(i,j)=vb2(i,j)/(sqrt(vb2(i,i)*vb2(j,j)));
    }}
  // Output
  return List::create(Named("mu")=mu,
                      Named("b1")=b1, Named("b2")=b2,
                      Named("hat")=fit,
                      Named("g1")=g1, Named("g2")=g2, 
                      Named("Vb1")=vb1, Named("Vb2")=vb2,
                      Named("GC1")=GC1, Named("GC2")=GC2,
                      Named("Ve")=ve, Named("h2")=h2);}