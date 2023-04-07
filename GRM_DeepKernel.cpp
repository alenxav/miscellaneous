// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>

// [[Rcpp::export]]
Eigen::MatrixXd ArcCosK(Eigen::MatrixXd X, bool centralizeX = true, int cores = 2){
  // cseweb.ucsd.edu/~saul/papers/nips09_kernel.pdf
  Eigen::setNbThreads(cores); int p = X.cols(), n = X.rows(); 
  double tmp, Npi=3.1416, theta, J1, Kij, Norm;
  if(centralizeX){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXd XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean()); XXp *= tmp;
  Eigen::VectorXd DiagXXp = XXp.diagonal().array();
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){ if(i>=j){ 
    Norm = sqrt(DiagXXp(i)*DiagXXp(j));
    theta = acos( XXp(i,j)/Norm);
    J1 = sin(theta) + (Npi-theta)*cos(theta);
    Kij = Norm/Npi*J1;
    XXp(i,j) = Kij*1.0; XXp(j,i) = Kij*1.0;}}}
  return XXp;}
