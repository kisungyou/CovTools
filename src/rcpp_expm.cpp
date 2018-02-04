#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cx_vec rcpp_geigen(arma::mat& A, arma::mat& B){
  arma::cx_vec out(A.n_cols,fill::zeros);
  eig_pair(out,A,B);
  return(out);
}
