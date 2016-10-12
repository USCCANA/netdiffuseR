// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rgraph_ba_cpp(
    const arma::sp_mat & graph,
    const arma::colvec & dgr, int m = 1, int t = 10, bool self=true) {

  // Creating the empty graph
  // m0: Size of the graph (changes throughout time).
  // n: Final size of the graph.
  int m0 = graph.n_cols;
  int n = m0 + t;

  // Creating new graph and vector of degree
  arma::sp_mat graph_new(n,n);
  graph_new.submat(0,0,m0-1, m0-1) = graph;

  arma::colvec dgr_new(n, arma::fill::zeros);
  dgr_new.subvec(0, m0-1) = dgr;

  // Start the process, K is sum(dgr)
  int K = sum(dgr);

  int m0_trunc;
  double randdraw, cump;

  // If self=true, then the prob are computed over m0+1, otherwise only over m0
  int extra = self? 1 : 0;
  for(int i=0;i<t;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // The number of conections is trucated by the number of vertices in the graph
    // so -m0_trunc- is actually -m-, but if currently there are less vertices than m,
    // then its truncated.
    m0_trunc = (m > m0)? m0 : m;
    // if (m > m0) m0_trunc = m0;

    int K0 = K;
    for (int j=0;j<m0_trunc;j++) {
      // Incrementing the degree of the one that is been added
      // one by one until having degree m0_trunk.
      dgr_new.at(m0) += 1.0;

      // std::cout << j << " Iter, "<< m << " m\n";
      // Random selection
      randdraw = unif_rand();

      // Calculating probabilities of been drawn. -cump- is the cumsum of the
      // probabilities
      cump = 0.0;
      for (int k=0; k<m0+extra; k++) {

        // In the case that in iter i the total degree is zero (no links)
        // then all individuals are equally likely to receive a new link.
        if (K0 != 0) cump += dgr_new.at(k)/(K + extra);
        else cump += 1/m0;

        // if (cump > 1) printf("over1 (%04d/%04d): %9.4g\n", k+1, m0+1,cump );

        // Links to the set of previous vertices
        if (randdraw <= cump) {
          graph_new.at(m0, k) += 1.0; // , graph_new.at(k, m0) += 1.0
          dgr_new.at(k) += 1.0;

          // Sumation of degrees
          K += 2; // outgoing + incoming (2)
          break;
        }
      }
    }
    ++m0;
  }

  return graph_new;
}

// [[Rcpp::export]]
arma::sp_mat rgraph_ba_new_cpp(int m0 = 1, int m = 1, int t = 10, bool self=true) {
  int n = m0;

  arma::sp_mat graph(n, n);
  arma::colvec dgr(n);

  if (!self) {
    graph.diag().fill(0.0);
    dgr.fill(0.0);
  }
  else {
    graph.diag().fill(1.0);
    dgr.fill(2.0);
  }

  return rgraph_ba_cpp(graph, dgr, m, t, self);
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
