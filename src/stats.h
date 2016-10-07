/*******************************************************************************
* stats.hpp header function for stats.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_STATS_
#define NETDIFFUSER_STATS_
using namespace Rcpp;

arma::mat cumulative_adopt_count_cpp(const arma::mat & cumadopt);

arma::rowvec hazard_rate_cpp(const arma::mat & cumadopt);

arma::colvec threshold_cpp(const arma::mat & exposure,
                           const arma::vec & toa, bool include_censored=false);

arma::sp_mat vertex_covariate_dist(const arma::sp_mat & graph,
                                   const NumericMatrix & X, double p=2.0);

arma::sp_mat vertex_covariate_compare(const arma::sp_mat & graph, const NumericVector & X,
                                      std::string symbol);

#endif
