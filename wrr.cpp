#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP wrr(NumericVector y, NumericMatrix gen, double model = 1){
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
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx);
  double Lmb = cxx;
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
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+b0*xx(j))/(Lmb*vx(j)+xx(j));
      e = e-gen(_,j)*(b[j]-b0);
    }
    if(model==1){
      ve = sum(e*y)/(n-1);
      va = sum(b*b+ve/(Lmb+xx))/p;
      Lmb = sqrt(cxx*ve/va);
      vx = sqrt(va/(b*b+ve/(xx+Lmb)));
    }
    if(model==2){
      ve = sum(e*y)/(n-1);
      va = sum(b*b+ve/(Lmb+xx))/p;
      Lmb = sqrt(cxx*ve/va);
      vx = va/(b*b+ve/(xx+Lmb));
    }
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  double h2 = cxx*va/(ve+cxx*va);
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu;}
  // Output
  return List::create(Named("mu")=mu,Named("b")=b,
                      Named("vb")=va,Named("ve")=ve,
                      Named("cxx")=cxx,Named("w")=vx,
                      Named("hat")=fit,Named("h2")=h2);}