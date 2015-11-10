// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat infection_cpp(
    NumericVector graph,
    const arma::colvec & times) {

  // Coersing a NumericVector into a cube for ease of use
  IntegerVector dims=graph.attr("dim");
  const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
  const int n = graph_cube.n_rows;
  const int T = graph_cube.n_slices;
  arma::mat infect(n,n, arma::fill::zeros);

  return infect;
}


/***R
library(netdiffuseR)
set.seed(123)
graph <- rand_graph()
graph <- array(unlist(lapply(1:3, function(x,...) graph)), dim=c(10,10,3))
times <- sample(1:3, 10, TRUE)
*/
