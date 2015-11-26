// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include "stats.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat scale_free_cpp(int m0 = 1, int m = 1, int t = 10) {

  // Creating the empty graph
  int n = m0 + t;
  arma::sp_mat graph(n,n);
  arma::colvec pattach(n, arma::fill::zeros);
  arma::colvec dgr(n, arma::fill::zeros);
  dgr.at(0) = 1;
  
  // Start the process
  int count = 1;
  for(int i=0;i<t;i++) {

    for(int j=0;j<count;j++) {
      // Calculating probabilies of attachment
      // dgr.head_rows(count) = degree_cpp(graph.submat(0, 0, count, count), 2, false, true);
      continue;
    }

    ++count;
  }

  return graph;
}
