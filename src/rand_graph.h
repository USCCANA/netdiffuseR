/*******************************************************************************
* adjmat.h header function for adjmat.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef DIFFTEST_RAND_GRAPH_
#define DIFFTEST_RAND_GRAPH_
using namespace Rcpp;

arma::colvec degree_cpp(
    const arma::mat & adjmat, const int & cmode=2,
    bool undirected=true, bool self=false);

arma::mat rand_graph_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);

arma::cube rand_dyn_graph_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);
#endif
