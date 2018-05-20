#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * (1) [2013.Cai] Optimal hypothesis testing
 */
// [[Rcpp::export]]
double rcpptest1_cai11(arma::mat& X){
  // 1-1. set params
  const int n = X.n_rows;
  const int p = X.n_cols;
  arma::rowvec Xi(p,fill::zeros);
  arma::rowvec Xj(p,fill::zeros);
  // 1-2. iteration
  double hXiXj = 0;
  for (int j=1;j<n;j++){
    Xj = X.row(j);
    for (int i=0;i<j;i++){
      Xi = X.row(i);
      hXiXj += std::pow(dot(Xi,Xj),2)-(dot(Xi,Xi)+dot(Xj,Xj)) + p;
    }
  }
  // 1-3. return hXiXj
  return(hXiXj);
}




