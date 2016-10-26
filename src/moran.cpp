// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> moran_cpp(const arma::colvec & x, const arma::sp_mat & w) {
  double xmean = mean(x);

  int N = x.n_rows;

  // Checking dims
  if (w.n_cols != w.n_rows) stop("-w- is not a square matrix");
  if (N != (int) w.n_cols) stop("-x- and -w- dimensions differ");

  // Weights sum
  double wsum = accu(w);

  arma::colvec xcent = x - xmean;
  double numer = accu((xcent * xcent.t()) % w);
  double denom = accu(pow(xcent, 2.0));

  /*/ Computing standard error
  double s1, s2, s3, s4, s5;
  arma::sp_mat::const_iterator b = w.begin();
  arma::sp_mat::const_iterator e = w.end();

  int i,j;
  for (arma::sp_mat::const_iterator iter=b; iter!=e; ++iter) {
    i = iter.col();
    j = iter.row();

    // Computing each count
    if (i >= j) s1+= pow(w.at(i,j) + w.at(j,i),2.0);

  }
  s1 = sqrt(s1);

  arma::colvec tmp(w.n_cols);
  for (int i=0; i < (int) w.n_rows; i++)
    tmp.at(i) = pow(accu(w.row(i)) + accu(w.col(i)), 2.0);

  s2 = sum(tmp);

  s3 = (1.0/N * accu(pow(xcent,4.0))) / pow(1.0/N * accu(pow(xcent,2.0)), 2.0);
  s4 = (N*N - 3.0*N + 3.0)*s1 - N*s2 + 3.0 * wsum * wsum;
  s5 = (N*N - N)*s1 - 2.0*N*s2 + 6.0*wsum*wsum;
  */
  std::vector<double> res(2);
  res[0] = (x.size()/wsum)*(numer/(denom + 1e-15));

  /*res[1] = (N*s4 - s3*s5)/((N - 1.0)*(N - 2.0)*(N - 3.0)*wsum*wsum) -
    pow(-1.0/(N - 1.0), 2.0);*/
  res[1] = 0;

  return res;
}

