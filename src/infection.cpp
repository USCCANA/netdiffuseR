// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat infection_cpp(
    NumericVector graph,
    const arma::colvec & times,
    bool normalize = true) {

  // Coersing a NumericVector into a cube for ease of use
  IntegerVector dims=graph.attr("dim");
  const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
  const int n = graph_cube.n_rows;
  const int T = graph_cube.n_slices;
  arma::colvec infect(n, arma::fill::zeros);

  // Variables to use within loop
  arma::mat graph_t(n,n);
  int ti, tj;

  for(int i=0;i<n;i++) {
    // Capturing variables
    ti = times(i);
    if (ti == T) continue;

    double numerator = 0.0;
    double denominator = 0.0;

    // For the adjusted verion, see the mathematical supplement on Valente et al. (2015)
    int nadopt_t = 0;

    graph_t = graph_cube.slice(ti);

    for(int j=0;j<n;j++) {
      if (i==j) continue;

      tj = times(j);
      // Adding up for t+1
      if (graph_t(j,i) != 0) {
        if (tj ==(ti+1)) {
          numerator += 1.0;
        }
        // Adding up for t+1 <= t <= T
        if (tj >= ti) denominator += 1.0;
      }

      // Has adopted so far?
      if (tj <= (ti+1)) ++nadopt_t;
    }

    // Putting all together
    infect(i) = (numerator / (denominator + 1e-10)) / (nadopt_t + 1e-10);
    if (normalize) infect(i) = infect(i) / (nadopt_t + 1e-10);
  }

  return infect;
}

// [[Rcpp::export]]
arma::colvec susceptibility_cpp(
    NumericVector graph,
    const arma::colvec & times,
    bool normalize = true) {

  // Coersing a NumericVector into a cube for ease of use
  IntegerVector dims=graph.attr("dim");
  const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
  const int n = graph_cube.n_rows;
  const int T = graph_cube.n_slices;
  arma::colvec suscep(n, arma::fill::zeros);

  // Variables to use within loop
  arma::mat graph_t(n,n);
  int ti, tj;

  for(int i=0;i<n;i++) {
    // Capturing variables
    ti = times(i);
    if (ti == 1) continue;

    double numerator = 0.0;
    double denominator = 0.0;

    // For the adjusted verion, see the mathematical supplement on Valente et al. (2015)
    int nadopt_t = 0;

    graph_t = graph_cube.slice(ti-1);

    for(int j=0;j<n;j++) {
      if (i==j) continue;

      tj = times(j);
      // Adding up for t+1
      if (graph_t(i,j) != 0) {
        if (tj ==(ti-1)) {
          numerator += 1.0;
        }
        // Adding up for t+1 <= t <= T
        if (tj < ti) denominator += 1.0;
      }

      // Has adopted so far?
      if (tj <= (ti+1)) ++nadopt_t;
    }

    // Putting all together
    suscep(i) = (numerator / (denominator + 1e-10));
    if (normalize) suscep(i) = suscep(i) / (nadopt_t + 1e-10);
  }

  return suscep;
}

/***R
library(netdiffuseR)
set.seed(123)
graph <- rand_graph(undirected = FALSE)
graph <- array(unlist(lapply(1:3, function(x,...) graph)), dim=c(10,10,3))
times <- sample(1:3, 10, TRUE)
infection(graph, times)
plot_diffnet(graph, toa_mat(times)$cumadopt)

# Paper data
graph <- edgelist_to_adjmat(matrix(c(5,1,5,2,5,3,5,4,6,5,7,5,8,5), byrow = TRUE, ncol=2), simplify = FALSE)
graph <- array(unlist(lapply(1:5, function(x,...) graph)), dim=c(8,8,5))
times <- c(1,2,5,5,3,4,5,5)
plot_diffnet(graph, toa_mat(times)$cumadopt)
infection(graph, times, TRUE)
susceptibility(graph, times, TRUE)
*/
