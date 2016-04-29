// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rewire_swap(
    const arma::sp_mat & graph, int nsteps=100,
    bool self=false, bool multiple=false,
    bool undirected=false) {

  // Clonning graph
  arma::sp_mat newgraph(graph);

  // Getting the indexes
  arma::umat indexes(graph.n_nonzero, 2);
  if (undirected) indexes = sparse_indexes(sp_trimatl(newgraph));
  else            indexes = sparse_indexes(newgraph);

  int s = 0;
  while (s < nsteps) {

    // Checking user interrupt
    if (s % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // Indexes
    int ij = floor(unif_rand()*indexes.n_rows);
    int i = indexes.at(ij, 0);
    int j = indexes.at(ij, 1);

    // If self, then remove from the indexes and go for the next
    if (!self && i==j) continue;

    // New end(s)
    int newij, newi, newj;

    // Choosing the other rand
    newij = floor( unif_rand()*indexes.n_rows );

    // The newij shouldn't be the same as ij!
    if (ij == newij) {
      s++;
      continue;
    }

    newi = indexes.at(newij, 0);
    newj = indexes.at(newij, 1);

    // If any of the edges coincides and no self edge is allowed,
    // then we must check for self pointing
    if (!self && ((i == newj) | (newi == j))) {
      s++;
      continue;
    }

    // Checking multiple
    if (!multiple && (newgraph.at(i, newj) != 0) | (newgraph.at(newi, j) != 0)) {
      s++;
      continue;
    }

    // Cleaning 2: Changing the graph patterns
    double v0 = newgraph.at(i, j);
    double v1 = newgraph.at(newi, newj);

    newgraph.at(i, j) = 0.0;
    if (undirected) newgraph.at(j, i) = 0.0;

    newgraph.at(newi, newj) = 0.0;
    if (undirected) newgraph.at(newj, newi) = 0.0;

    // Setting new values
    newgraph.at(i, newj) += v0;
    if (undirected) newgraph.at(newj, i) += v0;

    newgraph.at(newi, j) += v1;
    if (undirected) newgraph.at(j, newi) += v1;

    // And changing indexes
    indexes.at(ij,1) = newj;
    indexes.at(newij,1) = j;

    s++;
  }

  return newgraph;
}


/***R
x <- netdiffuseR::ring_lattice(10, 2)
set.seed(19237)

t(cbind(netdiffuseR::dgr(x, "indegree"), netdiffuseR::dgr(x, "outdegree")))
x2 <- rewire_swap(x, 1)
t(cbind(netdiffuseR::dgr(x2, "indegree"), netdiffuseR::dgr(x2, "outdegree")))

# Benchmark
library(igraph)
library(microbenchmark)

x <- barabasi.game(1e4)
y <- as_adj(x)

microbenchmark(
  ig = rewire(x, keeping_degseq(niter = 100)),
  nd = netdiffuseR:::rewire_swap(y, 1)
)
*/
