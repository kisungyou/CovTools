#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat covest_2003LW_computeF(arma::mat S, double rbar){
  // parameter and setup
  int N = S.n_rows;
  arma::mat matF(N,N,fill::zeros);

  // diagonal
  for (int i=0;i<N;i++){
    matF(i,i) = S(i,i);
  }
  // off-diagonal
  double odval = 0.0;
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      odval = rbar*std::sqrt(S(i,i)*S(j,j));
      matF(i,j) = odval;
      matF(j,i) = odval;
    }
  }
  // return the output
  return(matF);
}

// Y : (N,T), Ybar : (N), matS : (N,N)
// [[Rcpp::export]]
arma::mat covest_2003LW_computePi(arma::mat Y, arma::vec Ybar, arma::mat matS){
  // parameters
  int dimN = Y.n_rows;
  int dimT = Y.n_cols;

  // compute with many iterations..
  double tmpval = 0.0;
  double inners = 0.0;
  arma::mat matPi(dimN,dimN,fill::zeros);
  for (int i=0;i<dimN;i++){
    for (int j=0;j<dimN;j++){
      tmpval = 0.0;
      for (int t=0;t<dimT;t++){
        inners = (Y(i,t)-Ybar(i))*(Y(j,t)-Ybar(j)) - matS(i,j);
        tmpval += inners*inners;
      }
      matPi(i,j) = tmpval/static_cast<double>(dimT);
    }
  }

  // output
  return(matPi);
}

// Y : (N,T), Ybar : (N), matS : (N,N); only compute the second term for rho.hat object
// [[Rcpp::export]]
double covest_2003LW_computeRho(arma::mat Y, arma::vec Ybar, arma::mat matS, double rbar){
  // parameters
  int NN = Y.n_rows;
  int TT = Y.n_cols;
  double dTT = static_cast<double>(TT);

  // iterate
  double output = 0.0;
  double term1 = 0.0;
  double term2 = 0.0;
  for (int i=0;i<NN;i++){
    for (int j=0;j<NN;j++){
      if (i!=j){
        // reset the terms
        term1 = 0.0;
        term2 = 0.0;
        // iterate for theta.hat
        for (int t=0;t<TT;t++){
          term1 += ((Y(i,t)-Ybar(i))*(Y(i,t)-Ybar(i)) - matS(i,i))*((Y(i,t)-Ybar(i))*(Y(j,t)-Ybar(j)) - matS(i,j));
        }
        for (int t=0;t<TT;t++){
          term2 += ((Y(j,t)-Ybar(j))*(Y(j,t)-Ybar(j)) - matS(j,j))*((Y(i,t)-Ybar(i))*(Y(j,t)-Ybar(j)) - matS(i,j));
        }
        // adjust the terms
        term1 *= std::sqrt(matS(j,j)/matS(i,i))/dTT;
        term2 *= std::sqrt(matS(i,i)/matS(j,j))/dTT;

        // input
        output += (rbar/2.0)*(term1+term2);
      }
    }
  }

  // return output
  return(output);
}
