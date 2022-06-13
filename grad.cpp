// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <cmath>
using namespace Rcpp;
using namespace arma;
NumericVector llog_C(NumericVector z) {
  int n = z.size();
  double eps=1/50;
  
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    if (z[i]<eps){out[i] = log(eps) - 1.5 + 2 * z[i]/eps - 0.5 * (z[i]/eps)*(z[i]/eps); } else {
      out[i] = log(z[i]);}
    
  }
  return(out);
}
Rcpp::NumericVector arma2vec(arma::vec x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}
double obj_fun_rcpp (mat X, vec Y, vec lambda,vec beta){
  mat  U = X.each_col() % (Y-X* beta);
  
  int n =  U.n_rows;
  
  arma::vec temp1 = (ones(n)/(ones(n)+ (U) * (lambda)));
  NumericVector temp2 = llog_C(arma2vec(temp1));
  double out = sum(temp2);
  return(out);
}
// [[Rcpp::export]]
arma::vec def_grad(mat X, vec Y, vec lambda,vec beta) {
  int n = beta.size();

  vec out(n);
  
  for (int i = 0; i < n; ++i) {
    vec beta_temp = beta;
    beta_temp(i) = beta_temp(i)+pow (10.0, -8.0);
    out(i)=(obj_fun_rcpp(X,Y,lambda,beta_temp)-obj_fun_rcpp(X,Y,lambda,beta))/pow (10.0, -8.0);
    
  }
  return out;
}


