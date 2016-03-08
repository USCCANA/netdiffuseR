/*******************************************************************************
* egonets.h header function for egonets.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_EXTRA_
#define NETDIFFUSER_EXTRA_

using namespace Rcpp;

arma::umat sparse_indexes(const arma::sp_mat & mat);
#endif
