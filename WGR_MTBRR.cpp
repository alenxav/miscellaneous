#include <iostream>
using std::cout;
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP MRR(NumericMatrix Y, NumericMatrix X,
         int it = 1500, int bi = 500){
  
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
  
  cout << "\n Model setup: \n\n";
  cout << "Observations: " << n << "\n";
  cout << "Paramaters: " << p << "\n";
  cout << "Responses: " << k << "\n\n";
  
  // Pre-centralization
  NumericVector mu = colSums(y)/n;
  for(int j=0; j<k; j++){
    y(_,j) = (y(_,j)-mu(j))*o(_,j);  }  
  
  // Marker variance
  NumericMatrix GC(k,k);
  NumericMatrix xx(p,k), vx(p,k);
  double tmp;
  for(int i=0; i<p; i++){
    for(int j=0; j<k; j++){
      xx(i,j) = sum(X(_,i)*X(_,i)*o(_,j));
      tmp = sum(X(_,i)*o(_,j))/n(j);
      vx(i,j) = xx(i,j)/n(j)-tmp*tmp;}}
  NumericVector MSx = colSums(vx);  

  // Beta, intercept and residuals
  NumericMatrix b(p,k),vb(k,k),rho(k,k);
  NumericVector b0(k),vy(k),ve(k);
  double b1,lhs,rhs,covs,lambda;
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){vb(i,j) = 0.0;}}
  for(int i=0; i<k; i++){
    e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx(i);
    rho(i,i) = 1.0;}
  
  // Store Posterior
  NumericMatrix iG(k,k); iG = solve(vb);
  NumericMatrix B(p,k),VB(k,k);
  NumericVector VE(k);
  int mcmc;
  
  // Loop
  for(int numit=0; numit<it; numit++){    
    
    // Print iteration
    if ((numit+1) % 100 == 0){
      cout << "Iteration: " << numit+1 << "\n";}
    
    // Gauss-Seidel loop
    for(int j=0; j<p; j++){      
      
      // Current marker effect
      b0 = b(j,_);      
      
      // Update coefficient one trait at a time
      for(int i=0; i<k; i++){        
        
        // Compute coefficient covariances among traits
        covs = 0.0;
        for(int ii=0; ii<k; ii++){
          if(i!=ii){
            covs = covs+iG(i,ii)*b0(ii);
          }
        }
        
        //covs = covs/ve(i);
        covs = covs*ve(i);
        lambda = iG(i,i)*ve(i);        
        
        // Set LHS RHS
        lhs = xx(j,i) + lambda;
        rhs = sum(e(_,i)*X(_,j)) + xx(j,i)*b0(i) - covs;        
        
        // Update effects
        b1 = R::rnorm(rhs/lhs,sqrt(ve(i)/lhs));
        b(j,i) = b1;        
        
        // Update residuals
        e(_,i) = (e(_,i)-X(_,j)*(b1-b0(i)))*o(_,i);}}   
    
    // Residual variance components update
    for(int i=0; i<k; i++){
      ve(i) = sum(e(_,i)*e(_,i))/R::rchisq(n(i)-2);}
    
    // Get fitted values
    for(int i=0; i<n0; i++){for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j));}}      
    
    // GC correlations
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
        GC(i,j)= sum(fit(_,i)*fit(_,j))/sqrt( sum(fit(_,i)*fit(_,i)) * sum(fit(_,j)*fit(_,j))  );
      }
    }
    
    // Genetic variance
    for(int i=0; i<k; i++){ vb(i,i) = sum(b(_,i)*b(_,i))/R::rchisq(p-2);}
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}
    
    // Ridging diagonal and inverse variance
    for(int i=0; i<k; i++){ vb(i,i) = vb(i,i)*1.1; }
    iG = solve(vb);
    
    // Store posterior sums
    if(numit>bi){
      for(int i=0; i<k; i++){
        VB(_,i)=VB(_,i)+vb(_,i);
        B(_,i)=B(_,i)+b(_,i);}
      VE=VE+ve;
      mcmc=mcmc+1;
    }
    
    // End of MCMC iteration
  }
  
  // Compute posterior means
  VB = VB/mcmc;
  VE = VE/mcmc;
  B = B/mcmc;
  
  // Fitting the model
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      fit(i,j) = sum(X(i,_)*b(_,j))+mu(j);
    }}
  
  // Heritability
  NumericVector h2(k); 
  for(int j=0; j<k; j++){
    h2(j) = (VB(j,j)*MSx(j)) / (VB(j,j)*MSx(j) + VE(j) ); 
  }
  
  // Genetic correlations
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=VB(i,j)/(sqrt(VB(i,i)*VB(j,j)));
    }}  
  
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("h2")=h2,
                      Named("GC")=GC,
                      Named("Vb")=VB,
                      Named("Ve")=VE);}
