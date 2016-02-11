/*******************************************************************************
* selectionfunction.h header function for selectionfunction.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_SELECTIONFUNCTION_
#define NETDIFFUSER_SELECTIONFUNCTION_
using namespace Rcpp;

DataFrame select_egoalter_cpp(
    const arma::sp_mat & adjmat_t0,
    const arma::sp_mat & adjmat_t1,
    const NumericVector & adopt_t0,
    const NumericVector & adopt_t1);

#endif
