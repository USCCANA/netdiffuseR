/*******************************************************************************
* moran.h header function for moran.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_MORAN_
#define NETDIFFUSER_MORAN_

using namespace Rcpp;

double moran_cpp(const arma::colvec & x, const arma::sp_mat & w);
#endif
