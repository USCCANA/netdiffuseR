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
    if (cmode == 2) degree = degree - loop*2;
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
arma::mat exposure_cpp(
  List graph,
  arma::mat cumadopt,
  arma::mat attrs,
  bool outgoing = true,
  bool valued = true,
  bool normalized = true
) {

  // Getting parameters
  int T = graph.size();
  int n = cumadopt.n_rows;

  // Creating output and auxiliar objects
  arma::mat exposure(n,T);
  arma::colvec NUMERATOR(n);
  arma::colvec DENOMINATOR(n);

  for (int t=0;t<T;t++) {
    arma::sp_mat graph_t = graph[t];

    // Checking if need to fill with ones
    if (!valued) {
      graph_t = arma::spones(graph_t);
      graph_t.diag().zeros();
    }

    // Checking if incomming
    if (!outgoing) graph_t = graph_t.t();

    // Computing numerator
    NUMERATOR = graph_t * (attrs.col(t) % cumadopt.col(t));

    // Exposure for the time period
    if (normalized) {
      DENOMINATOR = graph_t * attrs.col(t) + 1e-15;
      exposure.col(t) = NUMERATOR / DENOMINATOR;
    } else {
      exposure.col(t) = NUMERATOR;
    }
  }

  // Returning
  return exposure;
}

/*
// [[Rcpp::export]]
arma::mat exposure_cpp(
    List graph, const arma::mat & cumadopt, int wtype = 0,
    double v = 1.0, bool undirected=true, bool normalized=true,
    int n=0, int T=0) {

  // Checking dimensions
  if (graph.size() != cumadopt.n_cols)
    stop("The length of -graph- and number of columns in -cumadopt- do not coincide.");

  // Variables initialization
  arma::mat exposure(n,T);
  arma::colvec NUMERATOR(n);
  arma::colvec DENOMINATOR(n);

  for(int t=0;t<T;t++) {

    arma::sp_mat graph_t = graph[t];

    // Computing weights
    if (wtype==0 || wtype == 1) { // Unweighted or SE
      NUMERATOR = graph_t*cumadopt.col(t);
      if (normalized) DENOMINATOR = arma::conv_to<arma::mat>::from(sum(graph_t,1)) + 1e-15;
    }
    else if (wtype > 1 && wtype <= 4) { // Degree
      arma::colvec degree = degree_cpp(graph_t, wtype - 2, undirected);

      NUMERATOR = graph_t*(cumadopt.col(t) % degree);
      if (normalized) DENOMINATOR = arma::conv_to<arma::mat>::from(sum(graph_t,1)) + 1e-15;
    }
    else {
      stop("Invalid weight code.");
    }

    // Filling the output
    if (normalized) exposure.col(t) = NUMERATOR / DENOMINATOR;
    else exposure.col(t) = NUMERATOR;

  }

  return exposure;
}
*/

/** *R
library(sna)
library(network)
library(netdiffuseR)
set.seed(123)
graph <- rand_dyn_graph_cpp(n=10,t=10)
adopt <- toa_mat_cpp(sample(1:dim(graph)[3], dim(graph)[2], TRUE))

exposure <- function(dynmat, adopt) {
  list(
    indegree=exposure_cpp(dynmat, adopt,2),
    structeq=exposure_cpp(dynmat, adopt,1),
    unweight=exposure_cpp(dynmat, adopt,0)
  )
}

# library(microbenchmark)
#
# microbenchmark(
#   old=ExposureCalc(graph, adopt$cumadopt),
#   new=exposure(graph, adopt$cumadopt),
#   times=100
# )
# Unit: microseconds
# expr      min        lq      mean    median        uq       max neval cld
#  old 4474.098 4621.5280 5190.7354 4756.6200 4926.9885 93214.670  1000   b
#  new   50.578   55.1375  101.2266   86.3295   97.3335  3711.223  1000  a
# 4756.6200/86.3295 = 55.09843
*/

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


/** *R
set.seed(123)
times <- sample(1:5, 10, TRUE)
adoptmat <- toa_mat_cpp(times)
adoptmat$adoptmat
new = cumulative_adopt_count_cpp(adoptmat$cumadopt)
old = netdifusseR::cumulativeAdopters(adoptmat$cumadopt)
new;old
microbenchmark::microbenchmark(
  old = netdiffusseR::cumulativeAdopters(adoptmat$cumadopt),
  new = cumulative_adopt_count_cpp(adoptmat$cumadopt)
)

*/

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

/** *R
old = hazardrate(adoptmat$cumadopt)
new = hazard_rate_cpp(adoptmat$cumadopt)
*/


// [[Rcpp::export]]
arma::colvec threshold_cpp(
    const arma::mat & exposure,
    const arma::vec & toa
    ) {

  int n = exposure.n_rows;
  // int T = exposure.n_cols;

  arma::colvec threshold(n);

  for(int i=0;i<n;i++) {
    // If NA (aka nan in Armadillo), then NA.
    if (!arma::is_finite(toa(i))) {
      threshold(i) = arma::datum::nan;
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
