// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include "stats.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat scale_free_cpp(int m0 = 1, int m = 1, int t = 10) {

  // Creating the empty graph
  int n = m0 + t;
  arma::sp_mat graph(n,n);
  arma::colvec dgr(n, arma::fill::zeros);
  dgr.at(0) = 2.0;
  graph.at(0,0) = 2.0;

  // Start the process, K is sum(dgr)
  int K = 2;
  for(int i=0;i<t;i++) {
    // The number of conections is trucated by the number of vertices in the graph
    int m1 = m;
    if (m > m0) m1 = m0;

    for (int j=0;j<m1;j++) {
      // Random selection
      double randdraw = unif_rand();

      // Calculating probabilities of been drawn. -cump- is the cumsum of
      double cump = 0.0;
      for (int k=0; k<m0; k++) {
        cump += dgr.at(k)/(K + 1);

        // Links to the set of previous vertices
        if (randdraw <= cump) {
          graph.at(m0, k) += 1.0, graph.at(k, m0) += 1.0;
          dgr.at(k) += 1.0, dgr.at(m0) += 1.0;

          // Sumation of degrees
          K += 2;
          break;
        }

        // Otherwise, it can link itself
        if ((k+1) == m0) {
          // printf("yes\n");
          graph.at(m0,m0) += 2.0;
          dgr.at(m0) += 2.0;

          // Sumation of degrees
          K += 2;
        }
      }
//       std::cout << "New iter\n";
    }
    // Increase the number of existing edges
    // std::cout << "Degree:" << dgr.t();
    ++m0;
  }

  return graph;
}

/***R
library(Matrix)
set.seed(123)
graph <- scale_free_cpp(t = 5000, m=1)
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
  new= netdiffuseR:::dgr.dgCMatrix(graph),
  old= sna:::degree(ig)
)

*/
