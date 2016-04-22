/*******************************************************************************
 rand_graph.cpp: FUNCTIONS FOR RANDOM GRAPH GENERATION

 created: nov 2, 2015
 copyright: MIT + LICENCE file
*******************************************************************************/


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rgraph_er_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {

  arma::sp_mat graph(n, n);

  for(int i=0;i<n;i++) {
    // Checling user interrup
    if (i % 200 == 0)
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
