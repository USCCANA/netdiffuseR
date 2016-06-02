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

using namespace Rcpp;

/* cmode:
 *  0: Indegree
 *  1: Outdegree
 *  2: Degree
 */
// [[Rcpp::export]]
arma::colvec degree_cpp(
    const arma::sp_mat & adjmat0, const int & cmode=2,
    bool undirected=true, bool self=false, bool valued=false) {

  if (cmode < 0 || cmode > 2) stop("Invalid degree");

  // Checking if it is valued or not
  int n = adjmat0.n_cols;
  arma::sp_mat adjmat(n,n);
  if (!valued) adjmat = arma::spones(adjmat0);
  else adjmat = adjmat0;

  // Calculating row/col sums
  arma::sp_mat degree(adjmat.n_cols, 1);
  if (cmode == 2) {
    if (undirected) degree = sum(adjmat,0).t();
    else degree = sum(adjmat,0).t() + sum(adjmat, 1);
  }
  else if (cmode == 0) degree = sum(adjmat, 0).t();
  else if (cmode == 1) degree = sum(adjmat, 1);

  // Checking if self or not
  arma::mat loop = arma::conv_to<arma::mat>::from(adjmat.diag());
  if (!self) {
    if (cmode == 2) {
      if (!undirected) degree = degree - loop*2;
      else             degree = degree - loop;
    }
    else degree = degree - loop;
  }

  arma::mat output = arma::conv_to< arma::mat >::from(degree);
  return output.col(0);
}

/* **R
 edgelist <- rbind(c(2,1),c(3,1),c(3,2))
 adjmat <- edgelist_to_adjmat_cpp(edgelist)
 degree_cpp(adjmat,0,FALSE)
 degree_cpp(adjmat,1,FALSE)
 degree_cpp(adjmat,2,FALSE)
*/

/* which:
 * 0: Unweighted
 * 1: Structure equiv weighted
 * 2: Indegree
 * 3: Outdegree
 * 4: Degree
 */

// [[Rcpp::export]]
arma::colvec exposure_cpp(
  const arma::sp_mat & graph,
  const arma::colvec & cumadopt,
  const arma::colvec & attrs,
  bool outgoing = true,
  bool valued = true,
  bool normalized = true
) {

  // Getting parameters
  unsigned n = cumadopt.n_rows;
  if (n != graph.n_cols) stop("-graph- is not squared.");

  // Creating output and auxiliar objects
  arma::colvec exposure(n);

  // Checking if need to fill with ones
  arma::sp_mat graph0(graph);
  if (!valued) graph0 = arma::spones(graph0);

  // Checking if incomming
  if (!outgoing) graph0 = graph0.t();

  // Computing numerator
  // NUMERATOR = (attrs.col(t) % cumadopt.col(t));
  arma::colvec NUMERATOR   = graph0 * (attrs % cumadopt);
  arma::colvec DENOMINATOR(n);

  // Exposure for the time period
  if (normalized) {
    DENOMINATOR = graph0 * attrs + 1e-15;
    exposure = NUMERATOR / DENOMINATOR;
  } else {
    exposure = NUMERATOR;
  }

  // Returning
  return exposure;
}

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
