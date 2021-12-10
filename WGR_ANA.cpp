#include <Rcpp.h>
#include <iostream>
#include <random>
using namespace Rcpp;

//' @title
//' Average Net algorithm
//' @description
//' New method for hyper-dimensional modeling 
//' @param y Numeric vector. Response variable, NA not allowed.
//' @param gen Numeric matrix. Prediction variables, NA not allowed.
//' @param h2 Numeric 0>X>1. Prior heritability.
//' @param alpha Numeric vector. Alpha parameters for multiple elastic-nets.
//' @param maxit Integer. Maximum number of interactions.
//' @param tol Numeric. Tolerance of absolute convergence.
//' @details Fits an average elastic-net, using 3 alpha values and lambda estimated from h2.
//' @export ANA

// [[Rcpp::export]]
SEXP ANA(NumericVector y,
         NumericMatrix gen,
         double h2 = 0.5,
         NumericVector alpha = NumericVector::create(0.1,0.01,0.001),
         int maxit = 1000, double tol = 10e-6){
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  int a = alpha.length();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  double ssy = sum(e*e);
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx);
  
  // Regulation coefficients
  NumericVector Lmb1(a);
  NumericVector Lmb2(a);
  for(int k=0; k<a; k++){
    Lmb1[k] = cxx*((1-h2)/h2)*alpha[k]*0.5;
    Lmb2[k] = cxx*((1-h2)/h2)*(1-alpha[k]);
  }
  NumericVector Gs(a);
  double OLS;
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j, J=0; J<p; J++){
      j = RGSvec[J];
      // Ordinary Least Square
      b0 = b[j];
      OLS = (sum(gen(_,j)*e)+xx[j]*b0);
      // Regularization of positive OLS
      if(OLS>0){
        // Regularization of positive OLS
        for(int k=0; k<a; k++){Gs[k] = (OLS-Lmb1[k])/(Lmb2[k]+xx(j));if(Gs[k]<0){Gs[k]=0;};}
      }else{
        // Regularization of negative OLS
        for(int k=0; k<a; k++){Gs[k] = (OLS+Lmb1[k])/(Lmb2[k]+xx(j));if(Gs[k]>0){Gs[k]=0;};}
      }
      b[j] = mean(Gs);
      // Residuals update
      e = e-gen(_,j)*(b[j]-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  double H2 = 1 - sum(y*e)/ssy;
  // Output
return List::create(Named("mu")=mu,Named("b")=b,
                    Named("hat")=fit,Named("h2")=H2,
                    Named("numit")=numit,Named("cnv")=cnv);}
