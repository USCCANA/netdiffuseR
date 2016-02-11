/*******************************************************************************
* rgraph_er.h header function for rgraph_er.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_RGRAPH_ER_
#define NETDIFFUSER_RGRAPH_ER_
using namespace Rcpp;

arma::sp_mat rgraph_er_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);

List rgraph_er_dyn_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);
#endif
