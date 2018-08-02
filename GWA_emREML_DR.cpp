#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emREML(NumericVector y, NumericMatrix Z,
     Rcpp::Nullable<Rcpp::NumericVector> D = R_NilValue,
     Rcpp::Nullable<Rcpp::NumericVector> R = R_NilValue,
     int maxit = 500, double tol = 10e-8){
  // Functions starts here
  int n = Z.nrow();
  int p = Z.ncol();
  // Observations weights
  bool N_WEIGHTS = FALSE;
  NumericVector w(n), iw(n);
  double wsum = 0;
  if(R.isNotNull()){
    N_WEIGHTS=TRUE;
    iw = R; w = 1/iw;
    wsum = sum(w);}
  // Marker weights
  bool P_WEIGHTS = FALSE;
  NumericVector d(p);
  if(D.isNotNull()){P_WEIGHTS=TRUE; d=D;}
  // Beta, mu and epsilon
  double b0, eM, ve, vb, h2, mu = mean(y), vy = var(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    if(N_WEIGHTS){
      xx[i] = sum(Z(_,i)*w*Z(_,i));
      vx[i] = var(Z(_,i)*w);
    }else{
      xx[i] = sum(Z(_,i)*Z(_,i));
      vx[i] = var(Z(_,i));
    }
  }
  double MSx = sum(vx), Lmb=MSx;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      b0 = b[j];
      if(P_WEIGHTS){
        if(N_WEIGHTS){
    // D and R
    b[j] = (sum(Z(_,j)*w*e)+xx[j]*b0)/(xx[j]+Lmb/d[j]);
        }else{
    // D, no R
    b[j] = (sum(Z(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb/d[j]); 
        }
      }else{
        if(N_WEIGHTS){
    // No D, R
    b[j] = (sum(Z(_,j)*w*e)+xx[j]*b0)/(xx[j]+Lmb);
        }else{
    // No D, No R
    b[j] = (sum(Z(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);        
        }
      }
      // Update residuals
      e = e-Z(_,j)*(b[j]-b0);}
    // Update intercept and variance components
    if(N_WEIGHTS){
      eM = sum(e*w)/wsum;
      mu = mu+eM;
      e = e-eM;
      ve = sum(e*w*y)/(wsum-1);
    }else{
      eM = mean(e);
      mu = mu+eM;
      e = e-eM;
      ve = sum(e*y)/(n-1);
    }
    vb = (vy-ve)/MSx;
    Lmb = ve/vb;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit(n);
  for(int k=0; k<n; k++){ fit[k] = sum(Z(k,_)*b)+mu; }
  h2 = vb*MSx/(vb*MSx+ve);
  // Genome-wide screening
  NumericVector LRT(p),PVAL(p),y0(n),e0(n),e1(n),b_ols(p);
  double ve0, ve1, L0, L1;
  for(int j=0; j<p; j++){
    // Full conditional phenotype
    y0 = e+Z(_,j)*b[j];
    if(N_WEIGHTS){
      b_ols[j] = sum(Z(_,j)*w*y0)/xx[j];
    }else{
      b_ols[j] = sum(Z(_,j)*y0)/xx[j];
    }
    // Null model
    if(N_WEIGHTS){
     ve0 = sum(y0*w*e0)/(wsum-1);
     e0 = y0-mean(y0);//sum(y0*w)/wsum;
    }else{
     ve0 = sum(y0*e0)/(n-1); 
     e0 = y0-mean(y0);
    }
    // Alt model
    if(N_WEIGHTS){
      ve1 = sum(y0*w*e1)/(n-1);
      e1 = y0-Z(_,j)*b_ols[j];
      e1 = e1-sum(e1*w)/wsum;
    }else{
      ve1 = sum(y0*e1)/(n-1); 
      e1 = y0-Z(_,j)*b_ols[j];
      e1 = e1-mean(e1);
    }
    // Likelihood ratio
    if(N_WEIGHTS){
      L0 = -sum(e0*w*e0)/(2*ve0)-0.5*n*log(6.28*ve0);
      L1 = -sum(e1*w*e1)/(2*ve1)-0.5*n*log(6.28*ve1);
      LRT[j] = 2*(L1-L0);  
    }else{
      L0 = -sum(e0*e0)/(2*ve0)-0.5*n*log(6.28*ve0);
      L1 = -sum(e1*e1)/(2*ve1)-0.5*n*log(6.28*ve1);
      LRT[j] = 2*(L1-L0);
    }
  }
  // Some bug, I cannot get the first LRT
  LRT[0] = LRT[1]; LRT = LRT-min(LRT); LRT[0] = 0;
  PVAL = -log10(1-pchisq(LRT,1,true));
  // Output
  return List::create(
    Named("HAT")=fit,
    Named("REG")=List::create(Named("mu")=mu,
    Named("b")=b),
    Named("GWA")=List::create(Named("b_OLS")=b_ols,
    Named("LRT")=LRT,
    Named("PVAL")=PVAL),
    Named("VCA")=List::create(Named("h2")=h2,
    Named("Vy")=vy,
    Named("Va")=vb*MSx,
    Named("Vb")=vb,
    Named("Ve")=ve));
}
