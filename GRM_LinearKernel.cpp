// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd GRM(Eigen::MatrixXd X, bool centralizeZ = true, int cores = 2){
  Eigen::setNbThreads(cores); int p = X.cols(); double tmp;
  if(centralizeZ){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXd XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean());
  XXp *= tmp; return XXp;}
