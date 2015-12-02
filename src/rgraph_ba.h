// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_RAND_GRAPH_
#define NETDIFFUSER_RAND_GRAPH_
using namespace Rcpp;

arma::sp_mat rgraph_ba_cpp(
    arma::sp_mat graph,
    arma::colvec dgr, int m = 1, int t = 10);

arma::sp_mat rgraph_ba_new_cpp(int m0 = 1, int m = 1, int t = 10);
#endif
