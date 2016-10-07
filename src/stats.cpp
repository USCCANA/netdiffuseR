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

// [[Rcpp::export]]
arma::mat cumulative_adopt_count_cpp(const arma::mat & cumadopt) {
  int n = cumadopt.n_rows;
  int T = cumadopt.n_cols;

  arma::mat adoptcount(3,T);

  // Computing cumulative adoptes
  adoptcount.row(0) = cumsum(sum(cumadopt,0));
  adoptcount.row(1) = adoptcount.row(0)/n;

  // Calculating rate
  adoptcount(2,0) = 0.0;
  for(int t=1;t<T;t++)
    adoptcount(2,t) = (adoptcount(1,t) - adoptcount(1,t-1))/(adoptcount(1,t-1) + 1e-10);

  return adoptcount;
}


// [[Rcpp::export]]
arma::rowvec hazard_rate_cpp(const arma::mat & cumadopt) {
  int n = cumadopt.n_rows;
  int T = cumadopt.n_cols;

  arma::rowvec cumadoptcount(T);
  arma::rowvec hazard(T);

  // Computing cumulative adoptes
  cumadoptcount = cumsum(sum(cumadopt,0));

  hazard(0) = 0.0;
  for(int t=1;t<T;t++)
    hazard(t) = (cumadoptcount(t) - cumadoptcount(t-1)) /
        (n - cumadoptcount(t-1) + 1e-10);

  return hazard;
}


// [[Rcpp::export]]
NumericVector threshold_cpp(
    const arma::mat & exposure,
    const arma::vec & toa,
    bool include_censored = false
    ) {

  int n = exposure.n_rows;
  // int T = exposure.n_cols;

  NumericVector threshold(n);

  for(int i=0;i<n;i++) {
    // If NA (aka nan in Armadillo), then NA.
    if (!arma::is_finite(toa(i))) {
      threshold(i) = arma::datum::nan;
      continue;
    }

    // If left censored and specified, then don't compute
    if ((toa(i)==1) & !include_censored) {
      threshold(i) = NA_REAL;
      continue;
    }

    threshold(i) = exposure(i,toa(i)-1);
  }

  return threshold;
}

/** *R

set.seed(123)
graph <- rand_dyn_graph_cpp(n=10,t=10)
tadopt <- sample(1:dim(graph)[3], dim(graph)[2], TRUE)

adopt <- toa_mat_cpp(tadopt)
exp_mat <- exposure(graph, adopt$cumadopt)$unweight

threshold_cpp(exp_mat, tadopt)

*/
/*
// [[Rcpp::export]]
double graph_density(arma::sp_mat graph, bool undirected=false) {
  int n = graph.n_cols;
  double dens = 0.0;

  return graph.n_nonzero / (n * (n-1));

}
*/

//' Computes p-norm between connected vertices
//' @param graph A square matrix of size \eqn{n} of class dgCMatrix.
//' @param X A numeric matrix of size \eqn{n\times K}{n * K}. Vertices attributes
//' @param p Numeric scalar. Norm to compute
//' @return A matrix of size \eqn{n\times n}{n*n} of class \code{dgCMatrix}. Will
//' be symmetric only if \code{graph} is symmetric.
//' @details For each par of vertices, the function computes the following
//' \deqn{%
//' D_{ij} = \left(\sum_{k=1}^K (X_{ik} - X_{jk})^{p} \right)^{1/p}\mbox{ if }graph_{i,j}\neq 0
//' }{%
//' D(i,j) = [\sum_k (X(i,k) - X(j,k))^p]^(1/p)  if graph(i,j) != 0
//' }
//' @export
//' @examples
//' set.seed(123)
//' G <- rgraph_ws(20, 4, .1)
//' X <- matrix(runif(40), ncol=2)
//'
//' vertex_covariate_dist(G, X)
// [[Rcpp::export]]
arma::sp_mat vertex_covariate_dist(const arma::sp_mat & graph, const NumericMatrix & X, double p=2.0) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  // Iterating over elements of graph
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++) {
    for (int k=0;k<X.ncol();k++)
      ans.at(iter.row(), iter.col()) += pow(fabs(X(iter.col(),k) - X(iter.row(),k)), p);

    ans.at(iter.row(), iter.col()) = pow(ans.at(iter.row(), iter.col()), 1/p);
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
