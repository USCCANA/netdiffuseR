/*******************************************************************************
 rand_graph.cpp: FUNCTIONS FOR RANDOM GRAPH GENERATION

 created: nov 2, 2015
 copyright: MIT + LICENCE file

 This file contains the following functions:
 - adopt_mat_cpp: Creates an adoption matrix of size nxT
 - edgelist_to_adjmat_cpp: Creates an adjacency matrix from an edgelist
 - adjmat_to_edgelist_cpp: The converse of the previous function
 - toa_mat_cpp: Creates a Time of Adoption Matrix of size nxT
 - isolated_cpp: Identifies the isolated nodes in a network
 - drop_isolated_cpp: Removes isolated networks from an adjmat
*******************************************************************************/


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rand_graph_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {
  arma::mat graph(n, n, arma::fill::zeros);

  // Using Rcpp (R's) RNG since it uses R's seed
  NumericVector datasource = runif(n*n);

  double w = 0.0;
  for(int i=0;i<n;i++) {

    /* Setting the length of the subloop acordingly to type of graph */
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;j++) {

      /* Assessing if include self */
      if (!self && (i==j)) continue;

      /* Setting the value of the tie */
      double val = datasource[i*n+j];
      w = val;

      if (val > (1-p)) {
        if (!weighted) w=1.0;
        graph(i,j) = w;
        if (undirected) graph(j,i) = w;
      }
    }
  }

  return graph;
}

// [[Rcpp::export]]
arma::cube rand_dyn_graph_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {

  arma::cube graphs(n,n,t);
  for(int i=0;i<t;i++)
    graphs.slice(i) = rand_graph_cpp(n, p, undirected, weighted, self);

  return graphs;

}

/* **R
 set.seed(123)
 rand_graph_cpp()

 rand_graph(undirected=FALSE)

 x <- rand_graph(1000, p=.1)
 sum(x)/length(x)
 */
