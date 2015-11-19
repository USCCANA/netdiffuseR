/*******************************************************************************
* adjmat.h header function for adjmat.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef DIFFTEST_ADJMAT_
#define DIFFTEST_ADJMAT_

using namespace Rcpp;

List adopt_mat_cpp(const IntegerVector & year);

arma::mat edgelist_to_adjmat_cpp(
    const arma::mat & data, NumericVector weights = NumericVector::create(),
    int n = 0,bool undirected = false, bool self = false, bool multiple = false);

arma::mat adjmat_to_edgelist_cpp(
    const arma::mat & adjmat, bool undirected = true);

IntegerMatrix toa_mat_cpp(const IntegerVector & year);

IntegerVector isolated_cpp(const arma::mat & adjmat, bool undirected=true);

arma::mat drop_isolated_cpp(const arma::mat & adjmat, arma::colvec isolated, bool undirected=true);

#endif
