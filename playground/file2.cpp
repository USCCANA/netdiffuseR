// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "adjmat.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cumulative_adopt(const arma::mat & adopt) {
  int n = adopt.n_rows;
  int T = adopt.n_cols;

  arma::mat cumadopt(3,T);

  // Computing cumulative adoptes
  cumadopt.row(0) = cumsum(sum(adopt,0));
  cumadopt.row(1) = cumadopt.row(0)/n;

  // Calculating rate
  cumadopt(2,0) = 0.0;
  for(int t=1;t<T;t++)
    cumadopt(2,t) = (cumadopt(1,t) - cumadopt(1,t-1))/cumadopt(1,t-1);

  return cumadopt;
}

/***R
set.seed(123)
times <- sample(1:5, 10, TRUE)
adoptmat <- adopt_mat_cpp(times)
adoptmat$adoptmat
cumulative_adopt(adoptmat$adoptmat1)
rowMeans(cumulative_adopt(adoptmat$adoptmat1))

diffusiontest::cumulativeAdopters(adoptmat$adoptmat1)
*/
