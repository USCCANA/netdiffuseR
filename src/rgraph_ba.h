/*******************************************************************************
* rgraph_ba.h header function for rgraph_ba.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_RGRAPH_BA_
#define NETDIFFUSER_RGRAPH_BA_
using namespace Rcpp;

arma::sp_mat rgraph_ba_cpp(
    arma::sp_mat graph,
    arma::colvec dgr, int m = 1, int t = 10, bool self=true);

arma::sp_mat rgraph_ba_new_cpp(int m0 = 1, int m = 1, int t = 10, bool self=true);

arma::sp_mat rgraph_sf_homo(
    const arma::colvec & eta, const arma::sp_mat & graph,
    const arma::colvec & dgr, int m = 1, int t = 10, bool self=true);

arma::sp_mat rgraph_sf_homo_new(
    const arma::colvec & eta, int m0 = 1, int m = 1, int t = 10, bool self=true);

#endif
