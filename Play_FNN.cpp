// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;

// Activation function
MatrixXd ReLU(MatrixXd A){
  for(int i=0; i<A.rows(); ++i){
    for(int j=0; j<A.cols(); ++j){
      //if(A(i,j)<0.0){A(i,j)=A(i,j)/100.0;} // Leaky ReLU
      if(A(i,j)<0.0){A(i,j)=0.0;} // Regular ReLU
      }}
  return A;}

// Fit node function
MatrixXd Node(MatrixXd H, MatrixXd W, VectorXd I){
  MatrixXd Hat(H.rows(),W.cols());
  Hat = H * W;
  for(int j=0; j<W.cols(); ++j){ Hat.col(j).array()+=I[j]; }
  return Hat;}

// Fit DNN function
MatrixXd FitDNN(MatrixXd X,
             MatrixXd W1, MatrixXd W2, MatrixXd W3,
             VectorXd I1, VectorXd I2, VectorXd I3){
  MatrixXd Hat(X.rows(),W3.cols());
  Hat = Node(ReLU(Node(ReLU(Node(X,W1,I1)),W2,I2)),W3,I3);
  return Hat;}

// Gauss-Seidel
// [[Rcpp::export]]
VectorXd CoordDesc(VectorXd y, MatrixXd X, VectorXd b,
                   VectorXd xx, double lmb, int GSit = 20){
  double b0, b1;
  VectorXd B = b;
  VectorXd e = y;
  for(int i=0; i<GSit; ++i){
    for(int j=0; j<X.cols(); ++j){
      b0 = B(j)+0.0;
      b1 = (e.transpose()*X.col(j)+xx(j)*b(j))/(xx(j)+lmb);
      e = e - X.col(j)*(b1-b0);
      B(j) = b1;}}
  return B;}

// Dropout function
MatrixXd Dropout(MatrixXd X, double doRate){
  VectorXd unif = VectorXd::Random(X.cols());
  for(int i=0; i<X.cols(); ++i){ unif[i] = (unif[i]+1.0)/2.0;}
  MatrixXd DO(X.rows(),X.cols());
  for(int i=0; i<X.cols(); ++i){
    if(unif[i]<doRate){
      DO.col(i) = VectorXd::Zero(X.rows());
    }else{DO.col(i) = X.col(i);}}
  return DO;}

// [[Rcpp::export]]
SEXP FNN(MatrixXd Y, // Phenotype
         MatrixXd X, // Genotype
         double rate = 0.1, // Learning rate
         int epochs = 10, // Full-data iterations
         int HL1=200, // Size of hidden layers 1
         int HL2=100, // Size of hidden layers 2
         double do1 = 0.50, // Droupout rate layer 1
         double do2 = 0.25, // Droupout rate layer 2
         double h2 = 0.8, // Heritability / Shrinkage
         int GSit = 20){ // Inner Gauss-Seidel iterations
  
  // Get dimensions of inputs
  int n = X.rows();
  int p = X.cols();
  int k = Y.cols();
  
  cout << "\n Network setup: \n\n";
  cout << "Observations: " << n << "\n";
  cout << "Paramaters: " << p << "\n";
  cout << "Responses: " << k << "\n";
  cout << "Nodes 1st layer: " << HL1 << "\n";
  cout << "Nodes 2nd layer: " << HL2 << "\n";
  cout << "Learning rate: " << rate << "\n";
  cout << "Dropout rate (1): " << do1 << "\n";
  cout << "Dropout rate (2): " << do2 << "\n";
  cout << "Number of epochs: " << epochs << "\n\n";
  
  // Memory placeholder for weights and intercept
  MatrixXd W1 = MatrixXd::Random(p,HL1); W1 *= (2/sqrt(p));
  VectorXd I1 = VectorXd::Random(HL1);
  MatrixXd W2 = MatrixXd::Random(HL1,HL2); W2 *= (2/sqrt(HL1));
  VectorXd I2 = VectorXd::Random(HL2);
  MatrixXd W3 = MatrixXd::Random(HL2,k); W3 *= (2/sqrt(HL2));
  VectorXd I3 = Y.colwise().mean();
  
  // Memory placeholder for fitted hidden layers
  MatrixXd H1(n,HL1); H1 = ReLU(Node(X,W1,I1));
  MatrixXd H2(n,HL2); H2 = ReLU(Node(H1,W2,I2));
  MatrixXd H3(n,k);   H3 = Node(H2,W3,I3);
  
  // Cross-products for Gauss-Seidel
  VectorXd xx = X.colwise().squaredNorm();
  VectorXd hh1 = H1.colwise().squaredNorm();
  VectorXd hh2 = H2.colwise().squaredNorm();
  
  // Shrinkage (L2 penalization)
  double shrk = (1-h2)/h2;
  double lmb1 = xx.mean()*shrk;
  double lmb2 = hh1.mean()*shrk;
  double lmb3 = hh2.mean()*shrk;
  
  // Memory placeholder for derivatives
  MatrixXd MSE(epochs,k);
  MatrixXd dH1(n,HL1);
  MatrixXd dH2(n,HL2);
  MatrixXd dH3(n,k);
  MatrixXd dW1(p,HL1);
  MatrixXd dW2(HL1,HL2);
  MatrixXd dW3(HL2,k);
  VectorXd dI1(HL1);
  VectorXd dI2(HL2);
  VectorXd dI3(k);
  
  // Loop to fit the DNN
  for(int i=0; i<epochs; ++i){
    
    cout << "Iteration " << i+1;
    
    /////////////
    // Layer 1 //
    /////////////
    
    cout << " (1) ";
    // Back-propagate residuals
    dH3 = Y-H3;
    dH2 = ReLU(dH3*W3.transpose()); 
    dH2 = Dropout(dH2,do2);
    dH1 = ReLU(dH2*W2.transpose());
    dH1 = Dropout(dH1,do1);
    // Compute gradient
    for(int j=0; j<HL1; ++j){
      dW1.col(j) = CoordDesc(dH1.col(j),X,W1.col(j),xx,lmb1,GSit);}
    dI1 = H1.colwise().mean();
    // Update 
    W1 = W1 + dW1*rate*(1-do1);
    I1 = I1 - dI1;
    
    /////////////
    // Layer 2 //
    /////////////
    
    cout << " (2) ";
    // Fit model
    H1 = ReLU(Node(X,W1,I1));
    H2 = ReLU(Node(H1,W2,I2));
    H3 = Node(H2,W3,I3);
    // Update cross products and regularizer
    hh1 = H1.colwise().squaredNorm();
    lmb2 = hh1.mean()*shrk;
    // Back-propagate residuals
    dH3 = Y-H3;
    dH2 = ReLU(dH3*W3.transpose());
    dH2 = Dropout(dH2,do2);
    // Compute gradient
    for(int j=0; j<HL2; ++j){
      dW2.col(j) = CoordDesc(dH2.col(j),H1,W2.col(j),hh1,lmb2,GSit);}
    dI2 = H2.colwise().mean();
    // Update 
    W2 = W2 + dW2*rate*(1-do2);
    I2 = I2 - dI2;
    
    /////////////
    // Layer 3 //
    /////////////
    
    cout << " (3) ";
    // Fit model
    H2 = ReLU(Node(H1,W2,I2));
    H3 = Node(H2,W3,I3);
    // Update cross products and regularizer
    hh2 = H2.colwise().squaredNorm();
    lmb3 = hh2.mean()*shrk;
    // Back-propagate residuals
    dH3 = Y-H3;
    // Compute gradient
    for(int j=0; j<k; ++j){
      dW3.col(j) = CoordDesc(dH3.col(j),H2,W3.col(j),hh2,lmb3,GSit);}
    dI3 = H3.colwise().mean();
    // Update 
    W3 = W3 + dW3*rate;
    I3 = I3 - dI3;
    
    //////////////////
    // Report Error //
    //////////////////
    
    cout << "\n";
    H3 = Node(H2,W3,I3);
    dH3 = Y-H3;
    MSE.row(i) = dH3.colwise().squaredNorm();
    cout << "MSE " << MSE.row(i) << "\n";
    
    
  }
  
  
  // Return list with weights and intercepts
  return Rcpp::List::create(Rcpp::Named("H1")=H1,
                            Rcpp::Named("H2")=H2,
                            Rcpp::Named("H3")=H3,
                            Rcpp::Named("W1")=W1,
                            Rcpp::Named("I1")=I1,
                            Rcpp::Named("W2")=W2,
                            Rcpp::Named("I2")=I2,
                            Rcpp::Named("W3")=W3,
                            Rcpp::Named("I3")=I3,
                            Rcpp::Named("MSE")=MSE);
  
}


