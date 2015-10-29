// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::cube cubeinout(arma::cube x) {
  return x;
}
/***R
cubeinout(array(1, c(2,2,2)))
*/
