/*******************************************************************************
* rgraph_ws.h header function for rgraph_ws,cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_RGRAPH_WS_
#define NETDIFFUSER_RGRAPH_WS_
using namespace Rcpp;

arma::umat sparse_indexes(const arma::sp_mat & mat);

arma::sp_mat ring_lattice(int n, int k);

arma::sp_mat rewire_graph_cpp(
    const arma::sp_mat & graph, double p,
    bool both_ends=false, bool self=false, bool multiple=false,
    bool undirected=false);

#endif
