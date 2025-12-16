
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;

// [[Rcpp::export]]
Rcpp::List pvp(Eigen::MatrixXd P, int test = 0) {
  MatrixXd Q = P.transpose();
  int n2 = Q.cols(), n = Q.rows();
  MatrixXd Va = Q * P;
  MatrixXd Pa(n, n2);
  MatrixXd Vainv = Va.inverse();

  if (test == 0) {
    Pa = Vainv * Q;
  }

  Rcpp::List output = Rcpp::List::create(Rcpp::Named("Va") = Va,
                                         Rcpp::Named("Pa") = Pa);
  return(output);
}
