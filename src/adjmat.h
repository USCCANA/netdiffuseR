/*******************************************************************************
* adjmat.h header function for adjmat.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_ADJMAT_
#define NETDIFFUSER_ADJMAT_

using namespace Rcpp;

void is_square_sp_mat(const arma::sp_mat & mat);

double min_int_na_cpp(const IntegerVector & x, const LogicalVector & isna);
double max_int_na_cpp(const IntegerVector & x, const LogicalVector & isna);

List toa_mat_cpp(const IntegerVector & year, int t0, int t1);

arma::sp_mat edgelist_to_adjmat_cpp(
    const arma::mat & edgelist, NumericVector weights = NumericVector::create(),
    int n = 0,bool undirected = false, bool self = false, bool multiple = false);

arma::mat adjmat_to_edgelist_cpp(
    const arma::sp_mat & adjmat, bool undirected = true);

IntegerMatrix toa_diff_cpp(const IntegerVector & year);

arma::icolvec isolated_cpp(const arma::sp_mat & adjmat, bool undirected=true);

arma::sp_mat drop_isolated_cpp(const arma::sp_mat & adjmat, arma::icolvec isolated, bool undirected=true);

#endif
