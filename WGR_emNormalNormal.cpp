#include <Rcpp.h>
using namespace Rcpp;
//' @name Regularized normal-normal model
//' @title Conjugated whole-genome regression
//' @description Machinery for hierarchical hyper-dimensional modeling.
//' @param y Numeric vector. Response variable, NA not allowed.
//' @param gen Numeric matrix. Prediction variables, NA not allowed.
//' @param h2 Numeric (0>X>1). Prior heritability. Set 0 to estimate variance components.
//' @param B Numeric vector. Prior regression coefficients. 
//' @param w Numeric (0>X>1). Weight of prior regression coefficients (BDA3, pg 42, eq 2.21).
//' @return List containing the intercept (mu), regression coefficients (b), fitted values (hat) and heritability (h2).
//' @details
//' Fits a ridge regression, if prior are provided it fits the model through a normal-normal model (BDA3, pg 40-43). If h2 is provided, lambda is estimated as: lambda=S*(1-h2)/h2, where scale is the sum of the variance of the SNPs, as S=2*sum(p*(1-p)).
//' Regression coefficients are updated via expectation-maximization (EM), which in this case is a variant of coordinate descent method to minimize the loss function based on the full-conditional optimization through Gauss-Seidel Residual Update (GSRU). GSRU is the consensus algorithm for whole-genome regression methods, being widely utilized in most the Bayesian and machine learning packages.
//' @export nnrr
//' @aliases nnrr
// [[Rcpp::export]]
SEXP nnrr(NumericVector y, NumericMatrix gen, double h2 = 0.5,
          Nullable<NumericVector> B = R_NilValue, double w = 0.5){
  // Convergence criteria
  int maxit = 200;
  double tol = 10e-7;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and residuals
  NumericVector b(p);
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Prior beta
  NumericVector A(p);
  if(B.isNotNull()){A = B;}
  NumericVector tilde(n);
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx);
  // Regulation coefficients
  double Lmb = cxx*((1-h2)/h2);
  if(h2==0){Lmb = cxx;};
  double b1;
  // Variance components
  double va = 1;
  double ve = 1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Full conditional y
      b0 = b[j];
      tilde = e + gen(_,j)*b0;
      // Update coefficient
      if(B.isNotNull()){
        b1 = (sum(gen(_,j)*tilde))/(Lmb+xx(j))*(1-w)+A[j]*w;
        }else{
        b1 = (sum(gen(_,j)*tilde))/(Lmb+xx(j));}
        b[j] = b1;
      // Update residuals
      e = tilde-gen(_,j)*b1;
    }
    // Alternative framework
    if(h2==0){
      // Variance components
      ve = sum(e*y)/(n-1);
      va = sum(b*b+(ve/(xx+Lmb)))/p;
      Lmb = sqrt(cxx*ve/va);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  if(h2==0){h2 = cxx*va/(ve+cxx*va);}
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu;}
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("h2")=h2);}
