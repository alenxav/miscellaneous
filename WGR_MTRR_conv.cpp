#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using std::cout;
// [[Rcpp::export]]
SEXP MV2(NumericMatrix Y,
         NumericMatrix X,
         int maxit = 200,
         double tol = 10e-10,
         double SOR = 1.0, // 0.75 is a good value
         double MultiplyOffDiag = 1.0, // 0.97 is a good value
         double MultiplyDiag = 1.0, // 1.03 is a good value
         double AddToDiag = 0.0, // 0.01 is a good value
         bool DecayAddToDiag = false){
  // Obtain environment containing function
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
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
  NumericMatrix xx(p,k), vx(p,k);
  double tmp;
  for(int i=0; i<p; i++){
    for(int j=0; j<k; j++){
      xx(i,j) = sum(X(_,i)*X(_,i)*o(_,j));
      tmp = sum(X(_,i)*o(_,j))/n(j);
      vx(i,j) = xx(i,j)/n(j)-tmp*tmp;}}
  NumericVector MSx = colSums(vx);
  // Beta, intersept and residuals
  NumericMatrix b(p,k),vb(k,k),iG(k,k),rho(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),vy(k),ve(k),RHS(k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){vb(i,j) = 0;}}
  for(int i=0; i<k; i++){
    e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/(n(i)-1);
    ve(i) = vy(i)*0.5;
    vb(i,i) = ve(i)/MSx(i);
    rho(i,i) = 1;}
  iG = solve(vb);
  // Convergence control
  NumericMatrix beta0(p,k);
  int numit = 0.0;
  double cnv = 1.0;
  double logtol = log(tol);
  // Genetic correlation
  NumericMatrix GC(k,k);
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
    GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  // Convergence analysis
  NumericVector AddToDiag0 = AddToDiag*(vy/MSx);
  NumericVector StoreConv(maxit), StoreConvGC(maxit), ve0(k);
  NumericVector StoreConvBV(maxit), StoreConvVC(maxit);
  NumericMatrix StoreH2(maxit,k), GC0(k,k), fit0(n0,k), vb0(k,k);
  // Loop
  while(numit<maxit){
    // Store pre-iteration
    beta0 = b+0.0; fit0 = fit+0.0;
    vb0 = vb+0.0; ve0 = ve+0.0;
    GC0 = GC+0.0;
    // Gauss-Seidel loop
    for(int j=0; j<p; j++){
      b0 = b(j,_);
      LHS = iG+0;
      for(int i=0; i<k; i++){
        LHS(i,i) = iG(i,i)+(xx(j,i)/ve(i));
        RHS(i) = (SOR*sum(e(_,i)*X(_,j))+(2.0-SOR)*xx(j,i)*b0(i))/ve(i);}
      // Update effects
      b1 = solve(LHS, RHS);
      b(j,_) = b1;
      // Update residuals
      for(int i=0; i<k; i++){
        e(_,i) = (e(_,i)-X(_,j)*(b1(i)-b0(i)))*o(_,i);}
    }
    // Fitting model
    for(int i=0; i<n0; i++){
      for(int j=0; j<k; j++){
        fit(i,j) = sum(X(i,_)*b(_,j));
      }
    }
    // Residual variance components update
    for(int i=0; i<k; i++){
      ve(i) = MultiplyDiag * (sum(e(_,i)*y(_,i)))/n(i);
      }
    // Genetic covariance components update
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
      // Diag VC
      if(i==j){
           //vb(i,j) = (1.01*vy(i)-ve(i))/MSx(i);
           tmp = (sum(fit(_,i)*y(_,j)) )/( n(i) *MSx(i) );
           vb(i,j) = tmp*MultiplyDiag + AddToDiag0(i);
        }else{
          if(i>j){
            tmp = (sum(fit(_,i)*y(_,j))+sum(fit(_,j)*y(_,i))) / ((n(i)*MSx(i))+(n(j)*MSx(j)));
            vb(i,j) = tmp*MultiplyOffDiag;
            vb(j,i) = vb(i,j);
          }
        }}}
    iG = solve(vb);
    // Genetic correlations
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
    // Decay on ridging & Successive Over Relaxation
    if (numit%5==0){if(SOR>1){SOR=SOR-0.01;}; if(SOR<1) SOR=SOR+0.01; } 
    if( DecayAddToDiag & (numit%20==0) & (AddToDiag>0) )  AddToDiag0 = (AddToDiag*sqrt(20/numit)) * (vy/MSx);
    // Convergence
    cnv = log(sum((beta0-b)*(beta0-b)));
    StoreConv[numit] = cnv;
    cnv = log(sum((vb0-vb)*(vb0-vb))+sum((ve0-ve)*(ve0-ve)));
    StoreConvVC[numit] = cnv;
    cnv = log(sum((GC0-GC)*(GC0-GC)));
    StoreConvGC[numit] = cnv;
    cnv = log(sum((fit0-fit)*(fit0-fit)));
    StoreConvBV[numit] = cnv;
    StoreH2(numit,_) = 1-ve/vy;
    ++numit;
    if(numit % 50 == 0){ cout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){break;}
  }
  // Fitting the model
  NumericVector h2(k); h2 = 1-ve/vy;
  for(int i=0; i<k; i++){ vb(i,i) = sum(fit(_,i)*y(_,i))/(n(i)*MSx(i)); }
  for(int i=0; i<n0; i++){for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j))+mu(j);}}
  // Output lists
  List convergence = List::create(Named("ConvCoef")=StoreConv,
                                  Named("ConvGC")=StoreConvGC,
                                  Named("ConvBV")=StoreConvBV,
                                  Named("ConvVC")=StoreConvBV,
                                  Named("AddToDiag")=AddToDiag0,
                                  Named("StoredH2")=StoreH2);
  List covariances = List::create(Named("h2")=h2,
                                  Named("Vb")=vb,
                                  Named("GC")=GC,
                                  Named("Ve")=ve,
                                  Named("Vy")=vy,
                                  Named("MSx")=MSx);
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("var")=covariances,
                      Named("cnv")=convergence);}
