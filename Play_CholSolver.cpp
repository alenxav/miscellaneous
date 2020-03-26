// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;

using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;

// [[Rcpp::export]]
SEXP CHOL(MatrixXd Y, // Phenotype
          MatrixXd X, // Genotype
          double h2 = 0.5){ 
  
  // Get dimensions of inputs
  int n = X.rows();
  int p = X.cols();
  int k = Y.cols();
  
  cout << "Observations: " << n << "\n";
  cout << "Paramaters: " << p << "\n";
  cout << "Responses: " << k << "\n";
  
  // Cross-products for Gauss-Seidel
  VectorXd xx = X.colwise().squaredNorm();
  
  // Shrinkage (L2 penalization)
  double shrk = (1-h2)/h2;
  double lmb1 = xx.mean()*shrk;
  
  // Factorize X
  cout << "Factorizing X for pre-conditioning\n";
  MatrixXd XpX(MatrixXd(p, p).setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
  for(int j=0; j<p; ++j){ XpX(j,j) =  XpX(j,j)+lmb1; }
  LLT<MatrixXd> CholX;
  CholX.compute(XpX);
  
  // Memory placeholder for weights and intercept
  MatrixXd W1 = MatrixXd(p,k);
  VectorXd I1 = Y.colwise().mean();
  for(int j=0; j<k; ++j){ Y.col(j).array()-=I1[j]; }
  
  // Solve
  cout << "Solve linear systems\n";
  VectorXd rhs(p);
  for(int j=0; j<k; ++j){
    rhs = X.transpose() * Y.col(j);
    W1.col(j) = CholX.solve(rhs);
    }

  // Fit
  cout << "Fit linear model\n";
  MatrixXd Hat(n,k);
  Hat = X * W1;
  for(int j=0; j<k; ++j){ Hat.col(j).array()+=I1[j]; }
    
  // Return list with weights and intercepts
  return Rcpp::List::create(Rcpp::Named("b")=W1,
                            Rcpp::Named("mu")=I1,
                            Rcpp::Named("hat")=Hat);
  
}


