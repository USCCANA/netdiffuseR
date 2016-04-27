/*******************************************************************************
* egonets.h header function for egonets.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_EXTRA_
#define NETDIFFUSER_EXTRA_

using namespace Rcpp;

arma::umat sparse_indexes(const arma::sp_mat & mat);

double angle(double x0, double y0, double x1, double y1);

arma::sp_mat sp_trimatl(const arma::sp_mat & x);

int unif_rand_w_exclusion(int n, int e);

#endif
