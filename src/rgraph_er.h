/*******************************************************************************
* adjmat.h header function for adjmat.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef DIFFTEST_RAND_GRAPH_
#define DIFFTEST_RAND_GRAPH_
using namespace Rcpp;

arma::sp_mat rgraph_er_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);

List rgraph_er_dyn_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);
#endif
