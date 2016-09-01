// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rewire_swap(
    const arma::sp_mat & graph, int nsteps=100,
    bool self=false, bool multiple=false,
    bool undirected=false, double pr_rewire=0.5,
    bool althexagons=false) {

  // Clonning graph
  arma::sp_mat newgraph(graph);

  // Getting the indexes
  arma::umat indexes(graph.n_nonzero, 2);
  if (undirected) indexes = sparse_indexes(sp_trimatl(newgraph));
  else            indexes = sparse_indexes(newgraph);

  double dens = graph.n_nonzero/(graph.n_cols*graph.n_cols);
  int s = 0;
  while (s++ < nsteps) {

    // Does it need to be rewired?
    if (unif_rand() > pr_rewire) continue;

    // Checking user interrupt
    if (s % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // Indexes
    int ij = floor(unif_rand()*indexes.n_rows);
    int i = indexes.at(ij, 0);
    int j = indexes.at(ij, 1);

    // If self, then remove from the indexes and go for the next
    if (!self && i==j) continue;

    // New end(s) (and values)
    int newij, newi, newj;
    int i1,i2,i3;
    double v0,v1,v2;

    // Choosing the other rand
    newij = floor( unif_rand()*indexes.n_rows );

    // The newij shouldn't be the same as ij!
    if (ij == newij) continue;

    newi = indexes.at(newij, 0);
    newj = indexes.at(newij, 1);

    bool ismultiple = !multiple &&
      (newgraph.at(i, newj) != 0) | (newgraph.at(newi, j) != 0);

    // Alternating Hexagons
    // Ramachandra Rao, et al, The Indian Journal of Statistics
    if (ismultiple && althexagons && (unif_rand() < 0.5)) {
      // Case 1: i1i2 and i2i3 (j == newi)
      if ((j==newi) && (newgraph.at(newj, i))) {

        i1 = i, i2 = j, i3 = newj;

        // All alt are zero
        if ((newgraph.at(i1,i3) != 0) ||
            (newgraph.at(i2,i1) != 0) ||
            (newgraph.at(i3,i2) != 0)) continue;

        v0 = newgraph.at(i1, i2);
        v1 = newgraph.at(i2, i3);
        v2 = newgraph.at(i3, i1);

        newgraph.at(i1, i2) = 0;
        newgraph.at(i2, i3) = 0;
        newgraph.at(i3, i1) = 0;

        newgraph.at(i1, i3) = v0;
        newgraph.at(i2, i1) = v1;
        newgraph.at(i3, i2) = v2;

      } else if ((i==newj) && (newgraph.at(j, newi))) {

        i1 = i, i2 = j, i3 = newi;

        // All alt are zero
        if ((newgraph.at(i1,i3) != 0) ||
            (newgraph.at(i2,i1) != 0) ||
            (newgraph.at(i3,i2) != 0)) continue;

        v0 = newgraph.at(i1, i2);
        v1 = newgraph.at(i2, i3);
        v2 = newgraph.at(i3, i1);

        newgraph.at(i1, i2) = 0;
        newgraph.at(i2, i3) = 0;
        newgraph.at(i3, i1) = 0;

        newgraph.at(i1, i3) = v0;
        newgraph.at(i2, i1) = v1;
        newgraph.at(i3, i2) = v2;

      }

      continue;
    }

    // Checking multiple
    if (ismultiple)
      continue;

    // If any of the edges coincides and no self edge is allowed,
    // then we must check for self pointing
    if (!self && ((i == newj) | (newi == j))) continue;

    // Cleaning 2: Changing the graph patterns
    v0 = newgraph.at(i, j);
    v1 = newgraph.at(newi, newj);

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
library(netdiffuseR)

set.seed(1133)
x <- barabasi.game(1e4)
y <- as_adj(x)

ind <- netdiffuseR:::sparse_indexes(y)

microbenchmark(
  ig      = rewire(x, keeping_degseq(niter = 100)),
  nd      = netdiffuseR:::rewire_swap(y, 100),
  # nd_fast = netdiffuseR:::rewire_swap_fast(y, ind, 100),
  unit="relative"
)

*/
