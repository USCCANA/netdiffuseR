#include <RcppArmadillo.h>
using namespace Rcpp;

// Returns a M x 2 matrix (M: # of non-zero elements) with the set of coordinates
// of the non-zero elements of an sparse matrix. Indices are from 0 to (n-1)
// [[Rcpp::export]]
arma::umat sparse_indexes(const arma::sp_mat & mat) {

  int n  = mat.n_nonzero;
  int nc = mat.n_cols;
  arma::umat indices(n,2);

  // If the matrix is empty (which makes no sense)
  if (!n) return indices;

  int curcol = 1;
  for (int i=0;i<n;i++) {

    // Incrementing column while there's no new element in the col
    while (mat.col_ptrs[curcol] <= i) ++curcol;
    // while (mat.col_ptrs[curcol-1] == mat.col_ptrs[curcol]) ++curcol;

    // Asigning indexes
    indices.at(i, 0) = mat.row_indices[i];
    // indices.at(i,1) = j;
    indices.at(i, 1) = curcol - 1;
  }
  // return indices;
  return indices;
}

