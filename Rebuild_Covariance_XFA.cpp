// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd XFA(const Eigen::MatrixXd& G, int pcs=3) {
  Eigen::VectorXd Var = G.diagonal();
  Eigen::VectorXd std_dev = Var.array().sqrt();
  Eigen::VectorXd inv_std_dev = std_dev.array().unaryExpr([](double v) {return (v > 0) ? 1.0 / v : 0.0;});
  Eigen::MatrixXd C = inv_std_dev.asDiagonal() * G * inv_std_dev.asDiagonal();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(C);
  if (eigen_solver.info() != Eigen::Success) {throw std::runtime_error("Eigendecomposition of the correlation matrix failed.");}
  Eigen::VectorXd eigenvalues = eigen_solver.eigenvalues();
  Eigen::MatrixXd eigenvectors = eigen_solver.eigenvectors();
  int n = eigenvalues.size();
  if (pcs > n) { pcs = n; }
  Eigen::MatrixXd V_reduced = eigenvectors.rightCols(pcs);
  Eigen::VectorXd D_reduced_diag = eigenvalues.tail(pcs);
  Eigen::MatrixXd NewC = V_reduced * D_reduced_diag.asDiagonal() * V_reduced.transpose();
  NewC.diagonal().setOnes();
  Eigen::MatrixXd NewG = std_dev.asDiagonal() * NewC * std_dev.asDiagonal();
  return NewG;
}

/*** R
G = bWGR::SimGC();
sigma = rchisq(50,10)
G = G * tcrossprod(sigma)
G_new = XFA(G,3);
plot(G,G_new);
*/