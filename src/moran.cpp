// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double moran_cpp(const arma::colvec & x, const arma::sp_mat & w) {
  double xmean = mean(x);
  int n = x.size();

  // Checking dims
  if (w.n_cols != w.n_rows) stop("-w- is not a square matrix");
  if (n != w.n_cols) stop("-x- and -w- dimensions differ");

  // Weights sum
  double wsum = accu(w);

  arma::colvec xcent = x - xmean;
  double numer = accu((xcent * xcent.t()) % w);
  double denom = accu(pow(xcent, 2.0));

  return (n/wsum)*(numer/(denom + 1e-15));
}

