// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rgraph_ba_cpp(
    arma::sp_mat graph,
    arma::colvec dgr, int m = 1, int t = 10) {

  // Creating the empty graph
  int m0 = graph.n_cols;
  int n = m0 + t;

  // Creating new graph and vector of degree
  arma::sp_mat graph_new(n,n);
  graph_new.submat(0,0,m0-1, m0-1) = graph;

  arma::colvec dgr_new(n, arma::fill::zeros);
  dgr_new.subvec(0, m0-1) = dgr;

  // Start the process, K is sum(dgr)
  int K = sum(dgr);
  // std::cout << dgr ;

  for(int i=0;i<t;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // The number of conections is trucated by the number of vertices in the graph
    // so -m1- is actually -m-, but if currently there are less vertices than m,
    // then its truncated.
    int m1 = m;
    if (m > m0) m1 = m0;

    for (int j=0;j<m1;j++) {
      // Incrementing the degree of the first
      dgr_new.at(m0) += 1.0;

      // std::cout << j << " Iter, "<< m << " m\n";
      // Random selection
      double randdraw = unif_rand();

      // Calculating probabilities of been drawn. -cump- is the cumsum of the
      // probabilities
      double cump = 0.0;
      for (int k=0; k<(m0+1); k++) {

        cump += dgr_new.at(k)/(K + 1);
        // if (cump > 1) printf("over1 (%04d/%04d): %9.4g\n", k+1, m0+1,cump );

        // Links to the set of previous vertices
        if (randdraw <= cump) {
          graph_new.at(m0, k) += 1.0; // , graph_new.at(k, m0) += 1.0
          dgr_new.at(k) += 1.0;

          // Sumation of degrees
          K += 2;
          break;
        }
      }
    }
    ++m0;
  }

  return graph_new;
}

// [[Rcpp::export]]
arma::sp_mat rgraph_ba_new_cpp(int m0 = 1, int m = 1, int t = 10) {
  int n = m0;
  arma::sp_mat graph(n, n);
  graph.diag() = arma::ones(n);

  arma::colvec dgr(n, arma::fill::ones);
  dgr = dgr*2.0;
  return rgraph_ba_cpp(graph, dgr, m, t);
}


/* **R
library(Matrix)
set.seed(123)
graph <- rgraph_ba_cpp(t = 500, m=1)
# graph
library(sna)
library(SparseM)
gplot(methods::as(graph, "matrix.csc"))

library(netdiffuseR)
deg <- dgr(graph)
deg2 <- degree(methods::as(graph, "matrix.csc"))

all(table(deg2) == table(deg))

library(microbenchmark)
ig <- methods::as(graph, "matrix.csc")
microbenchmark(
  new= netdiffuseR::dgr(graph),
  old= sna::degree(ig), times=1000
)

*/
