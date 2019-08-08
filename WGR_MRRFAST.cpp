#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP mrrFast(NumericMatrix Y, NumericMatrix X){
  int maxit = 200; double tol = 10e-8;  
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];  
  // Functions starts here
  int k = Y.ncol(), p = X.ncol(), n0 = X.nrow();
  NumericMatrix fit(n0,k),o(n0,k),y(n0,k),e(n0,k);
  for(int i=0; i<k; i++){
    o(_,i) = ifelse(is_na(Y(_,i)),0,1);
    y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));  }
  NumericVector n = colSums(o);  
  // Mu
  NumericVector mu = colSums(y)/n;
  for(int j=0; j<k; j++){
    y(_,j) = (y(_,j)-mu(j))*o(_,j);  }  
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
  NumericMatrix b(p,k),vb(k,k),rho(k,k);
  NumericVector b0(k),vy(k),ve(k);
  double b1,lhs,rhs,covs,vb0,lambda;
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx(i);
    rho(i,i) = 1;}
  NumericMatrix iG = solve(vb);  
  // Convergence control
  NumericMatrix bc(p,k);
  int numit = 0;
  double cnv = 1;  
  // Loop
  while(numit<maxit){    
    // Gauss-Seidel loop
    bc = b+0;
    for(int j=0; j<p; j++){      
      // Current marker effect
      b0 = b(j,_);      
      // Update coefficient one trait at a time
      for(int i=0; i<k; i++){        
        // Compute coefficient covariances among traits
        covs = 0.0;
        for(int ii=0; ii<k; ii++){
          if(i!=ii){ covs = covs+iG(i,ii)*b0(ii); } }
        //covs = covs/ve(i);
        covs = covs*ve(i);
        lambda = iG(i,i)*ve(i);        
        // Set LHS RHS
        lhs = xx(j,i) + lambda;
        rhs = sum(e(_,i)*X(_,j)) + xx(j,i)*b0(i) - covs;        
        // Update effects
        b1 = rhs/lhs;
        b(j,i) = b1;        
        // Update residuals
        //e(_,i) = ifelse( o(_,i)==1, (e(_,i)-X(_,j)*(b1-b0(i))), 0 );
        e(_,i) = (e(_,i)-X(_,j)*(b1-b0(i)))*o(_,i); }}    
    // Update variance components every few iterations
    if(numit%3==0){      
      // Get fitted values
      for(int i=0; i<n0; i++){
        for(int j=0; j<k; j++){
          fit(i,j) = sum(X(i,_)*b(_,j));}}      
      // Residual variance components update
      for(int i=0; i<k; i++){ ve(i) = (sum(e(_,i)*y(_,i)))/(n(i)-1); }      
      // Genetic covariance components update
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i==j){
                //vb(i,i) = sum(fit(_,i)*y(_,i))/(n(i)*MSx(i));
                vb(i,i) = (1.1*vy(i)-ve(i))/MSx(i);
              }else{
            if(i>j){
                vb0 = (sum(fit(_,i)*y(_,j))+sum(fit(_,j)*y(_,i))) / ((n(i)*MSx(i))+(n(j)*MSx(j)));
                vb(i,j) = vb0; vb(j,i) = vb0;}}}}      
      // Generalized inverse of G2A_Kernels
      iG = solve(vb);}    
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
                      Named("Vb")=vb, Named("GC")=GC,
                      Named("Ve")=ve);}