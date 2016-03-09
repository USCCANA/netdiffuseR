/*******************************************************************************
* struct_equiv.h
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_STRUCT_EQUIV_
#define NETDIFFUSER_STRUCT_EQUIV_

using namespace Rcpp;

List struct_equiv_cpp(const arma::sp_mat & gdist, double v = 1.0);

#endif
