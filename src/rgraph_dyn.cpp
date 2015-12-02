// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "rgraph_ba.h"
// #include "stats.h"

using namespace Rcpp;

// // [[Rcpp::export]]
// arma::sp_mat rgraph_dyn_next_cpp(arma::sp_mat & graph, arma::colvec & dgr, double a=.3, double p=.3, double r=.3) {
//
// }

// [[Rcpp::export]]
double persistant(const List & graph, int i, int j, int n) {
  int T = graph.size();

  double p = 0.0;
  for(int t=(T-1);t>0;t--) {

    arma::sp_mat graph_t1 = graph[t];
    arma::sp_mat graph_t0 = graph[t-1];
    bool prod = ((graph_t0.at(i,j) - graph_t1.at(i,j))==0)  && graph_t0.at(i,j);
    if (prod) p += 1.0;
  }

  return p;
}

/***R
library(netdiffuseR)
set.seed(123)
graph <- rgraph_ba()
graph <- list(graph, graph)
graph[[2]][1,] <- 1
graph[[3]] <- rgraph_ba()
graph
sapply(0:10, function(x) persistant(graph, 0,x ,11))
*/
