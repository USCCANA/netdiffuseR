/*******************************************************************************
* rewire_swap.h header function for rewire_swap.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_REWIRE_SWAP_
#define NETDIFFUSER_REWIRE_SWAP_

using namespace Rcpp;

arma::sp_mat rewire_swap(
    const arma::sp_mat & graph, int nsteps=100,
    bool self=false, bool multiple=false,
    bool undirected=false);

#endif
