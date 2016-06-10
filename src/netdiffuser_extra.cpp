#include <RcppArmadillo.h>
using namespace Rcpp;

// Returns a M x 2 matrix (M: # of non-zero elements) with the set of coordinates
// of the non-zero elements of an sparse matrix. Indices are from 0 to (n-1)
// [[Rcpp::export]]
arma::umat sparse_indexes(const arma::sp_mat & mat) {

  int n  = mat.n_nonzero;
  arma::umat indices(n,2);

  // If the matrix is empty (which makes no sense)
  if (!n) return indices;

  // More efficient implementation
  arma::sp_mat::const_iterator begin = mat.begin();
  arma::sp_mat::const_iterator end   = mat.end();

  int i = 0;
  for (arma::sp_mat::const_iterator it = begin; it != end; ++it) {
    indices.at(i, 0) = it.row();
    indices.at(i++, 1) = it.col();
  }

  // return indices;
  return indices;
}

//   I  ) (x(i) < x(j)) & (y(i) < y(j)) = + Works fine iff
//   II ) (x(i) > x(j)) & (y(i) < y(j)) = - must add pi
//   III) (x(i) > x(j)) & (y(i) > y(j)) = +
//   IV ) (x(i) < x(j)) & (y(i) > y(j)) = - (works ok)
// [[Rcpp::export]]
double angle(double x0, double y0, double x1, double y1) {
  // Computing distances and angles
  double xdist = x1 - x0;
  double ydist = y1 - y0;
  double alpha = atan(ydist/(xdist)+1e-15);

  // Setting cases
  if      ((xdist < 0) && (ydist > 0)) return(alpha + PI);
  else if ((xdist < 0) && (ydist < 0)) return(alpha + PI);

  return(alpha);
}


// [[Rcpp::export]]
arma::sp_mat sp_trimatl(const arma::sp_mat & x) {
  // Getting start-end
  arma::sp_mat::const_iterator start = x.begin();
  arma::sp_mat::const_iterator end   = x.end();

  // Empty mat
  int n = x.n_cols;
  arma::sp_mat newx(n,n);

  for (arma::sp_mat::const_iterator it = start; it != end; ++it) {
    int i = it.row();
    int j = it.col();
    if (i >= j)
      newx.at(i,j) = *it;
  }

  return newx;
}

// [[Rcpp::export]]
arma::sp_mat sp_diag(const arma::sp_mat & x, const arma::vec & v) {
  // Checking dimensions
  if (x.n_cols != x.n_rows) Rcpp::stop("-x- must be square.");
  if (x.n_cols != v.n_elem) Rcpp::stop("length(v) must be equal to ncol(x)");

  arma::sp_mat out(x);
  out.diag() = v;

  return out;
}


// [[Rcpp::export]]
int unif_rand_w_exclusion(int n, int e) {
  double r = unif_rand();
  int num = floor(r*(n-1));
  if (e != n) {
    if (num >= e) ++num;
  }
  return(num);
}
