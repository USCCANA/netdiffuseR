#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "adjmat.hpp"
#include "struct_equiv.hpp"

using namespace Rcpp;

/* which:
 * 0: Unweighted
 * 1: Structure equiv weighted
 * 2: Indegree
 * 3: Outdegree
 * 4: Degree
 */

// [[Rcpp::export]]
arma::mat exposure_cpp(
    NumericVector graph, const arma::mat & cumadopt, int wtype = 0,
    double v = 1.0, bool undirected=true) {

  // Coersing a NumericVector into a cube for ease of use
  IntegerVector dims=graph.attr("dim");
  const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
  const int n = graph_cube.n_rows;
  const int T = graph_cube.n_slices;
  arma::mat exposure(n,T);

  // Initializing containers
  List se(3);
  arma::mat semat(n,n);

  arma::colvec NUMERATOR(n);
  arma::colvec DENOMINATOR(n);
  arma::mat graph_t(n,n);

  // Filling the first and last with zeros
  // may skip this in the future.
  exposure.col(0).fill(0.0);
  exposure.col(T).fill(0.0);

  for(int t=0;t<(T-1);t++) {
    graph_t = graph_cube.slice(t);
    // Computing weights

    if (wtype==0) { // Unweighted
      NUMERATOR = graph_t*cumadopt.col(t);
      DENOMINATOR = sum(graph_t,1) + 1e-15;
    }
    else if (wtype == 1) {// SE

      // Calculating the inverse of the SE distance
      se = struct_equiv_cpp(graph_t, v, false, true);
      semat       = (as< arma::mat >(se["SE"]));

      NUMERATOR   = semat * cumadopt.col(t);
      DENOMINATOR = sum(semat, 1) + 1e-15;
    }
    else if (wtype > 1 && wtype <= 4) { // Degree
      arma::colvec degree = degree_cpp(graph_t, wtype - 2, undirected);

      NUMERATOR = graph_t*(cumadopt.col(t) % degree);
      DENOMINATOR = sum(graph_t,1) + 1e-15;
    }
    else {
      stop("Invalid weight code.");
    }

    // Filling the output
    exposure.col(t+1) = NUMERATOR / DENOMINATOR;

  }

  return exposure;
}


/***R
library(sna)
library(network)
library(diffusiontest)
set.seed(123)
graph <- rand_dyn_graph_cpp(n=10,t=10)
adopt <- adopt_mat_cpp(sample(1:dim(graph)[3], dim(graph)[2], TRUE))

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

/** *R
set.seed(123)
times <- sample(1:5, 10, TRUE)
adoptmat <- adopt_mat_cpp(times)
adoptmat$adoptmat
new = cumulative_adopt_count_cpp(adoptmat$cumadopt)
old = diffusiontest::cumulativeAdopters(adoptmat$cumadopt)
new;old
microbenchmark::microbenchmark(
  old = diffusiontest::cumulativeAdopters(adoptmat$cumadopt),
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
    const arma::vec & toe
    ) {

  int n = exposure.n_rows;
  int T = exposure.n_cols;

  arma::colvec threshold(n);

  for(int i=0;i<n;i++) {
    threshold(i) = exposure(i,toe(i)-1);
  }

  return threshold;
}

/** *R

set.seed(123)
graph <- rand_dyn_graph_cpp(n=10,t=10)
tadopt <- sample(1:dim(graph)[3], dim(graph)[2], TRUE)

adopt <- adopt_mat_cpp(tadopt)
exp_mat <- exposure(graph, adopt$cumadopt)$unweight

threshold_cpp(exp_mat, tadopt)

*/
