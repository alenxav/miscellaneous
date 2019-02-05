#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP MVRR(NumericMatrix Y, NumericMatrix X, bool Choleski = false){
  
  // Convergence parameters
  int maxit = 350;
  double tol = 10e-8;
  
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  Rcpp::Function chol2inv = base["chol2inv"];
  
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
  
  //NumericVector MSx = colSums(xx);
  NumericVector MSx = colSums(vx);
  
  // Beta, intersept and residuals
  NumericMatrix b(p,k),vb(k,k),iG(k,k),rho(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),eM(k),mu(k),vy(k),ve(k),RHS(k);
  mu = colSums(y)/n;
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = (y(_,i)-mu(i))*o(_,i);
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx(i);
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
    for(int i=0; i<k; i++){
      ve(i) = sum(e(_,i)*y(_,i))/(n(i)-1);
      if(Choleski){
        vb(i,i) = (1.1*vy(i)-ve(i))/MSx(i);
      }else{ vb(i,i) = (vy(i)-ve(i))/MSx(i);}
    }
    
    // Approximate genetic correlation
    for(int i=0; i<n0; i++){ 
      for(int j=0; j<k; j++){
        fit(i,j) = sum(X(i,_)*b(_,j));}}
    for(int i=0; i<k; i++){ 
      for(int j=0; j<k; j++){
        rho(i,j) = sum(fit(_,i)*fit(_,j))/sqrt(sum(fit(_,i)*fit(_,i))*sum(fit(_,j)*fit(_,j)));
      }}
    
    // Covariance components
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i>j){
          vb(i,j) = rho(i,j)*sqrt(vb(i,i)*vb(j,j));
          vb(j,i) = vb(i,j); }}}
    for(int i=0; i<k; i++){vb(i,i)=vb(i,i)*1.01;}
    if(Choleski){iG = chol2inv(vb);}else{iG = solve(vb);}
    
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  
  // Fitting the model
  for(int i=0; i<n0; i++){ 
    for(int j=0; j<k; j++){
      fit(i,j) = sum(X(i,_)*b(_,j))+mu(j);}}
  
  // Heritability
  NumericVector h2(k); 
  for(int i=0; i<k; i++){ h2 = (vb(i,i)*MSx(i))/((vb(i,i)*MSx(i))+ve); }
  
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("hat")=fit, Named("h2")=h2,
                      Named("Vb")=vb, Named("Ve")=ve,
                      Named("Vy")=vy, Named("MSx")=MSx);}