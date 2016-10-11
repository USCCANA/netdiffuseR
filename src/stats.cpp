/*******************************************************************************
 stats.cpp: FUNCTIONS FOR CALCULATING NETWORK STATISTICS

 created: nov 2, 2015
 copyright: MIT + LICENCE file

 This set of functions helps handling network data, in particular, diffussion
 network data. Importing edgelist as adjacency matrices and generating auxiliary
 matrices used for the model.

 This file contains the following functions:
  - toa_mat_cpp: Creates an adoption matrix of size nxT
  - edgelist_to_adjmat_cpp: Creates an adjacency matrix from an edgelist
  - adjmat_to_edgelist_cpp: The converse of the previous function
  - toa_mat_cpp: Creates a Time of Adoption Matrix of size nxT
  - isolated_cpp: Identifies the isolated nodes in a network
  - drop_isolated_cpp: Removes isolated networks from an adjmat
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"

using namespace Rcpp;



//' @export
//' @rdname vertex_covariate_dist
// [[Rcpp::export]]
arma::sp_mat vertex_covariate_dist(const arma::sp_mat & graph,
                                   const arma::mat & X,
                                   double p = 2.0) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  // Iterating over elements of graph
  int i,j;
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++) {

    i = iter.row();
    j = iter.col();

    ans.at(i,j) = pow(
      arma::as_scalar((X.row(i) - X.row(j)) * (X.row(i) - X.row(j)).t()),
      1/p
    );
  }

  return ans;
}

// [[Rcpp::export]]
arma::sp_mat vertex_mahalanobis_dist_cpp(const arma::sp_mat & graph,
                                   const arma::mat & X,
                                   const arma::mat & S) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  arma::mat Sinv = S.i();

  // Iterating over elements of graph
  int i,j;
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++) {

    i = iter.row();
    j = iter.col();

    ans.at(i,j) = pow(
      arma::as_scalar((X.row(i) - X.row(j)) * Sinv * (X.row(i) - X.row(j)).t()),
      0.5
    );
  }

  return ans;
}

//' Comparisons at dyadic level
//' @param graph A matrix of size \eqn{n\times n}{n*n} of class \code{dgCMatrix}.
//' @param X A numeric vector of length \eqn{n}.
//' @param funname Character scalar. Comparison to make (see details).
//' @details
//'
//' This auxiliary function takes advantage of the sparcity of \code{graph} and
//' applies a function in the form of \eqn{funname(x_i,x_j)}{funname(X[i],X[j])}
//' only to \eqn{(i,j)} that have no empty entry. In other words, applies a compares
//' elements of \code{X} only between vertices that have a link; making
//' \code{nlinks(graph)} comparisons instead of looping through \eqn{n\times n}{n*n},
//' which is much faster.
//'
//' \code{funname} can take any of the following values:
//' \code{"distance"}, \code{"^2"} or \code{"quaddistance"}, \code{">"} or \code{"greater"},
//' \code{"<"} or \code{"smaller"}, \code{">="} or \code{"greaterequal"},
//' \code{"<="} or \code{"smallerequal"}, \code{"=="} or \code{"equal"}.
//' @return A matrix \code{dgCMatrix} of size \eqn{n\times n}{n*n} with values in
//' the form of \eqn{funname(x_i,x_j)}{funname(X[i],X[j])}.
//' @examples
//'
//' # Basic example ------------------------------------------------------------
//' set.seed(1313)
//' G <- rgraph_ws(10, 4, .2)
//' x <- rnorm(10)
//'
//' vertex_covariate_compare(G, x, "distance")
//' vertex_covariate_compare(G, x, "^2")
//' vertex_covariate_compare(G, x, ">=")
//' vertex_covariate_compare(G, x, "<=")
//' @export
// [[Rcpp::export]]
arma::sp_mat vertex_covariate_compare(const arma::sp_mat & graph, const NumericVector & X,
                                      std::string funname) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  // Iterating over elements of graph
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++)
    ans.at(iter.row(),iter.col()) = fun(X[iter.row()], X[iter.col()]);

  return ans;
}
