/*******************************************************************************
* egonets.h header function for egonets.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_EGONETS_
#define NETDIFFUSER_EGONETS_

using namespace Rcpp;

double egonet_attrs_cpp(
    const arma::sp_mat & graph, const arma::uvec V,
    NumericMatrix attrs, bool outer=true, bool self=true, bool valued=true);
#endif
