/*******************************************************************************
* stats.hpp header function for stats.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_STATS_
#define NETDIFFUSER_STATS_
using namespace Rcpp;

arma::colvec degree_cpp(
    const arma::sp_mat & adjmat0, const int & cmode=2,
    bool undirected=true, bool self=false, bool valued=false);

arma::mat exposure_cpp(
    List graph, arma::mat cumadopt, arma::mat attrs,
    bool outgoing =true, bool valued=true, bool normalized=true);

arma::mat cumulative_adopt_count_cpp(const arma::mat & cumadopt);

arma::rowvec hazard_rate_cpp(const arma::mat & cumadopt);

arma::colvec threshold_cpp(const arma::mat & exposure,
                           const arma::vec & toa);

#endif
