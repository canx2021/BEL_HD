//--------------------------------------------------------------
// Header (header)
//--------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <cmath>
using namespace Rcpp;
using namespace arma;


//--------------------------------------------------------------
// Functions (Functions_cpp)
//--------------------------------------------------------------
arma::vec colSums(arma::mat mat) {
  int ncols = mat.n_cols;
  vec res = randn(ncols);
  for (int i=0; i<ncols; i++) {
    res(i) = sum(mat.col(i));
  }
  return(res);
}

arma::vec rowSums(arma::mat mat) {
  int nrows = mat.n_rows;
  vec res = randn(nrows);
  for (int i=0; i<nrows; i++) {
    res(i) = sum(mat.row(i));
  }
  return(res);
}
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
NumericVector llogp_C(NumericVector z) {
  int n = z.size();
  double eps=1/50;
  
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    if (z[i]<eps){out[i] = 2/eps - z[i]/eps/eps; } else {
      out[i] = 1/(z[i]);}
    
  }
  return(out);
}


Rcpp::NumericVector arma2vec(arma::vec x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

double penlik_new(arma::vec beta,
                  arma::vec tau,
                  arma::mat U,
                  arma::vec lambda){
  double pi = - sum(beta % beta/tau/2);
  arma::vec a = (1/(1+U * lambda));
  NumericVector aa = llog_C(arma2vec(a));
  double res = pi+sum(aa);
  return(res);
}

double obj_fun_rcpp (NumericVector q,NumericMatrix U, NumericVector v){
  int n =  as<arma::mat>(U).n_rows;
  arma::vec temp1 = ((ones(n)+ as<arma::mat>(U) * as<arma::vec>(q)));
  NumericVector temp2 = llog_C(arma2vec(temp1));
  double out = -sum(temp2)+sum( as<arma::vec>(q)% as<arma::vec>(q)/2/ as<arma::vec>(v));
  return(out);
}

NumericVector obj_grad_fun_rcpp (NumericVector q,NumericMatrix U, NumericVector v){
  int n =  as<arma::mat>(U).n_rows;
  arma::vec temp1 = as<arma::mat>(U) * as<arma::vec>(q) +ones(n);
  NumericVector temp2 = llogp_C(arma2vec(temp1));
  arma::vec out = -colSums(as<arma::mat>(U).each_col() % as<arma::vec>(temp2))+as<arma::vec>(q)/as<arma::vec>(v);
  return(arma2vec(out));
}

NumericMatrix arma2mat (arma::mat x) {
  NumericMatrix y = wrap(x) ;
  return(y) ;
}

// [[Rcpp::export]]
Rcpp::List Act_grad_cpp(
    arma::mat U_old,
    arma::mat X,
    arma::mat Sigma_lp,
    
    arma::vec lambda_old,
    arma::vec beta_old,
    arma::vec tau,
    arma::vec v,
    arma::vec Y,
    arma::vec mu_lp,
    
    int p
){
  Rcpp::Environment MASS("package:MASS");
  Rcpp::Function mvrnorm = MASS["mvrnorm"];
  
  Rcpp::Environment mvtnorm("package:mvtnorm");
  Rcpp::Function dmvnorm = mvtnorm["dmvnorm"];
  
  
  Rcpp::Environment stats("package:stats");
  Rcpp::Function rnorm = stats["rnorm"];
  Rcpp::Function optim = stats["optim"];
  Rcpp::Function runif = stats["runif"];
  
  NumericVector mu_lp_vec = arma2vec(mu_lp);
  NumericMatrix Sigma_lp_vec = arma2mat(Sigma_lp);
  
  NumericVector beta_new = mvrnorm(1,mu_lp_vec,Sigma_lp_vec);
  
  mat  U = X.each_col() % (Y-X* as<arma::vec>(beta_new));
  
  vec lambda_int = 0.0001*ones(p);
  NumericVector lambda_int_vec = arma2vec(lambda_int);
  
  NumericMatrix U_mat = arma2mat(U);
  NumericVector v_vec = arma2vec(v);
  
  NumericVector lambda= arma2vec(lambda_old) - obj_grad_fun_rcpp(arma2vec(lambda_old),U_mat,v_vec)*pow (10.0, -20.0);
  
  NumericVector gnew = dmvnorm(beta_new,mu_lp_vec,Sigma_lp_vec);
  NumericVector gold = dmvnorm(arma2vec(beta_old),mu_lp_vec,Sigma_lp_vec); 
  double fold = penlik_new(beta_old,tau,U_old,lambda_old);
  double fnew = penlik_new(beta_new,tau,U,lambda);
  
  double r = fnew + (log(gold))[0] - fold - (log(gnew))[0];
  
  NumericVector nun = (runif(1));
  double un = log(nun)[0];
  if(r>un){
    return List::create(Named("accept") = 1 , 
                        Named("beta") = beta_new ,
                        Named("lambda") = lambda ,
                        Named("U") = U);
  }else{
    return List::create(Named("accept") = 0 , 
                        Named("beta") = beta_old ,
                        Named("lambda") = lambda_old ,
                        Named("U") = U_old);
  }
  
  
}




