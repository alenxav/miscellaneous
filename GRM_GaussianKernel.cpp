// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd GAU(Eigen::MatrixXd X, double phi = 1.0, int cores = 2){
  Eigen::setNbThreads(cores);
  int n = X.rows(); double tmp;
  Eigen::MatrixXd XXp = X*X.transpose();
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){ if(i>j){
    tmp = sqrt(XXp(i,i) + XXp(j,j) - 2*XXp(i,j));
    XXp(i,j) = tmp*1.0; XXp(j,i) = tmp*1.0;}}};
  for(int i=0; i<n; i++){XXp(i,i) = 0.0;}
  tmp = phi * (-n*(n-1)) / (XXp.colwise().sum()).sum();
  XXp *= tmp; return exp(XXp.array());}
