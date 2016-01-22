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
arma::sp_mat rgraph_er_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {

  arma::sp_mat graph(n, n);

  GetRNGstate();
  for(int i=0;i<n;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    /* Setting the length of the subloop acordingly to type of graph */
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;j++) {

      /* Assessing if include self */
      if (!self && (i==j)) continue;

      /* Setting the value of the tie */
      double w = unif_rand();
      if (w > (1-p)) {
        if (!weighted) w=1.0;
        graph.at(i,j) = w;
        if (undirected) graph.at(j,i) = w;
      }
    }
  }
  PutRNGstate();

  return graph;
}

// [[Rcpp::export]]
List rgraph_er_dyn_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {

  List graphs(t);
  for(int i=0;i<t;i++)
    graphs[i] = rgraph_er_cpp(n, p, undirected, weighted, self);

  return graphs;

}

/* **R
 set.seed(123)
 rand_graph_cpp()

 rand_graph(undirected=FALSE)

 x <- rand_graph(1000, p=.1)
 sum(x)/length(x)
 */
