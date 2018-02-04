#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double rcpp_objective(arma::mat& S, arma::mat& X,arma::mat& Z, double lambda){
  double term1 = arma::trace(S*X);
  double term2 = log(arma::det(X));
  double term3 = lambda*arma::norm(vectorise(Z),1);
  double output = term1-term2+term3;
  return(output);
}
double rcpp_shrinkage(double a, double kappa){
  double term1=0;
  double term2=0;
  if (a>kappa){    term1 = a-kappa;  }
  if (a<-kappa){   term2 = -a-kappa; }
  double output = term1-term2;
  return(output);
}

/*
* 1. ADMM for Sparse Precision Matrix Estimation
*    Core Part from Sparse Precision Matrix, in
*
*    minimize  trace(S*X) - log det X + lambda*||X||_1
*
*    with input S      : sample covariance
*               lambda : parameter
*               rho and alpha set as 1.0
*/
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat rcpp_ADMMprecision(arma::mat& S, const double lambda){
  // 1. parameters and set up
  const int max_iter  = 1000;
  const double abstol = 1e-6;
  const double reltol = 1e-3;
  const int    n      = S.n_cols;

  arma::mat X(n,n,fill::zeros);
  arma::mat X_hat(n,n,fill::zeros);
  arma::mat Z(n,n,fill::zeros);
  arma::mat Zold(n,n,fill::zeros);
  arma::mat U(n,n,fill::zeros);

  double rho   = 1.0;
  double alpha = 1.0;

  arma::colvec es(n,fill::zeros);
  arma::colvec xi(n,fill::zeros);
  arma::mat Q(n,n,fill::zeros);

  arma::vec objval(max_iter,fill::zeros);
  arma::vec r_norm(max_iter,fill::zeros);
  arma::vec s_norm(max_iter,fill::zeros);
  arma::vec eps_pri(max_iter,fill::zeros);
  arma::vec eps_dual(max_iter,fill::zeros);

  for (int k=0;k<max_iter;k++){
    // update X
    eig_sym(es,Q,rho*(Z-U)-S);
    for (int i=0;i<n;i++){
      xi(i) = (es(i)+sqrt(pow(es(i),2)+4*rho))/(2*rho);
    }
    X = Q*arma::diagmat(xi)*Q.t();

    // update Z with relazation
    Zold = Z;
    X_hat = alpha*X + (1-alpha)*Zold;
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        Z(i,j) = rcpp_shrinkage(X_hat(i,j)+U(i,j), lambda/rho);
      }
    }

    // update U
    U = U + (X_hat-Z);

    // diagnostics
    objval(k) = rcpp_objective(S, X, Z, lambda);
    r_norm(k) = arma::norm(X-Z,"fro");
    s_norm(k) = arma::norm(-rho*(Z-Zold),"fro");
    if (norm(X,"fro")>norm(Z,"fro")){
      eps_pri(k) = static_cast<double>(n)*abstol + reltol*norm(X,"fro");
    } else {
      eps_pri(k) = static_cast<double>(n)*abstol + reltol*norm(Z,"fro");
    }
    eps_dual(k) = static_cast<double>(n)*abstol + reltol*norm(rho*U, "fro");

    if ((r_norm(k)<eps_pri(k))&&(s_norm(k)<eps_dual(k))){
      break;
    }
  }
  return(X);
}
