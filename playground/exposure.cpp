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
arma::mat exposure_cpp(NumericVector graph, const arma::mat & adopt,
                  int wtype = 0) {

  // Coersing a NumericVector into a cube for ease of use
  IntegerVector dims=graph.attr("dim");
  const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
  const int n = graph_cube.n_rows;
  const int T = graph_cube.n_slices;
  arma::mat exposure(n,T, arma::fill::zeros);

  // Initializing containers
  List se(3);
  arma::colvec degree(n);
  arma::mat semat(n,n);

  arma::colvec NUMERATOR(n);
  arma::colvec DENOMINATOR(n);
  arma::mat graph_t(n,n);

  for(int t=0;t<(T-1);t++) {
    graph_t = graph_cube.slice(t);
    // Computing weights

    if (wtype==0) { // Unweighted
      NUMERATOR = graph_t*adopt.col(t);
      DENOMINATOR = sum(graph_t,1);
    }
    else if (wtype == 1) {// SE

      // Calculating the inverse of the SE distance
      se = struct_equiv_cpp(graph_t, 1.0, true, true);
      semat       = (as< arma::mat >(se["d"]));
      NUMERATOR   = semat * adopt.col(t);
      DENOMINATOR = sum(semat, 1);
    }
    else if (wtype > 1 && wtype <= 4) { // Degree
      degree = degree_cpp(graph_t, wtype - 2);
      NUMERATOR = graph_t*(adopt.col(t) % degree);
      DENOMINATOR = sum(graph_t,1);
    }
    else {
      stop("Invalid weight code.");
    }

    // Filling the output
    exposure.col(t+1) = NUMERATOR / (DENOMINATOR + 1e-10);

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

library(microbenchmark)

microbenchmark(
  old=ExposureCalc(graph, adopt$adoptmat),
  new=exposure(graph, adopt$adoptmat),
  times=100
)
# Unit: microseconds
# expr      min        lq      mean    median        uq       max neval cld
#  old 4474.098 4621.5280 5190.7354 4756.6200 4926.9885 93214.670  1000   b
#  new   50.578   55.1375  101.2266   86.3295   97.3335  3711.223  1000  a
# 4756.6200/86.3295 = 55.09843
*/
