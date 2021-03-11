// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::isnan;

// [[Rcpp::export]]
SEXP EigenMV(MatrixXd Y, MatrixXd X, int maxit = 500, 
             double tol = 10e-10, double deflate = 0.95){
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();

  // Incidence matrix Z
  MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  VectorXd n = Z.colwise().sum();
  VectorXd iN = n.array().inverse();
  
  // Centralize y
  VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X
  MatrixXd XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  VectorXd MSx = XSX.colwise().sum();
  VectorXd TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  VectorXd ve = vy * 0.5;
  VectorXd iVe = ve.array().inverse();
  MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  MatrixXd iG = vb.inverse();
  VectorXd h2 = 1 - ve.array()/vy.array();
  
  // Beta tilde;
  MatrixXd tilde = X.transpose() * y;
  
  // Initialize coefficient matrices
  MatrixXd LHS(k,k);
  VectorXd RHS(k);
  MatrixXd b = MatrixXd::Zero(p,k);
  VectorXd b0(k), b1(k);
  MatrixXd e(n0,k); e = y*1.0;
  
  // Convergence control
  MatrixXd beta0(p,k);
  VectorXd CNV(maxit);
  double cnv = 10.0;
  int numit = 0;
  double logtol = log10(tol);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Gauss-Seidel loop
    for(int j=0; j<p; j++){

      // Update coefficient
      b0 = b.row(j)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(j).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(j).transpose()*e).array() + XX.row(j).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(j) = b1;
      
      // Update residuals
      e = (e-(X.col(j)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();

    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ vb(i,i) = TildeHat(i,i)/TrXSX(i); }else{
          vb(i,j) = deflate*(TildeHat(i,j)+TildeHat(j,i))/(TrXSX(i)+TrXSX(j));}}}
    iG = vb.inverse();
    
    // Print status
    cnv = log10((beta0.array()-b.array()).square().sum());
    CNV(numit) = cnv;
    ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){break;}
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations
  MatrixXd GC(k,k);
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
    GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("cnv")=cnv);
  
}

