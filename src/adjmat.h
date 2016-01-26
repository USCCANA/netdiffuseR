/*******************************************************************************
* adjmat.h header function for adjmat.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef DIFFTEST_ADJMAT_
#define DIFFTEST_ADJMAT_

using namespace Rcpp;

List adopt_mat_cpp(const IntegerVector & year);

arma::sp_mat edgelist_to_adjmat_cpp(
    const arma::mat & data, NumericVector weights = NumericVector::create(),
    int n = 0,bool undirected = false, bool self = false, bool multiple = false);

arma::mat adjmat_to_edgelist_cpp(
    const arma::sp_mat & adjmat, bool undirected = true);

// arma::mat adjmat_to_dyn_edgelist_cpp(NumericVector adjmat, bool undirected=true);

IntegerMatrix toa_mat_cpp(const IntegerVector & year, int t0, int t1);

arma::icolvec isolated_cpp(const arma::sp_mat & adjmat, bool undirected=true);

arma::sp_mat drop_isolated_cpp(const arma::sp_mat & adjmat, arma::icolvec isolated, bool undirected=true);

#endif
