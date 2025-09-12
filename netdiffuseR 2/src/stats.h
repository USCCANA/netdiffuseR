/*******************************************************************************
* stats.hpp header function for stats.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_STATS_
#define NETDIFFUSER_STATS_
using namespace Rcpp;

arma::sp_mat vertex_covariate_dist(
    const arma::sp_mat & graph,
    const NumericMatrix & X,
    double p=2.0
);

arma::sp_mat vertex_mahalanobis_dist_cpp(
    const arma::sp_mat & graph,
    const arma::mat & X,
    const arma::mat & W
  );

arma::sp_mat vertex_covariate_compare(
    const arma::sp_mat & graph,
    const NumericVector & X,
    std::string symbol
  );

List moran_cpp(
    const arma::colvec & x,
    const arma::sp_mat & w
  );

List struct_equiv_cpp(
    const arma::sp_mat & gdist,
    double v = 1.0
);
#endif
