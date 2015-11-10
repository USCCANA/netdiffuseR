/*******************************************************************************
* stats.hpp header function for stats.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifndef DIFFTEST_STATS_
#define DIFFTEST_STATS_
using namespace Rcpp;

arma::colvec degree_cpp(
    const arma::mat & adjmat, const int & cmode=2,
    bool undirected=true, bool self=false);

arma::mat exposure_cpp(
    NumericVector graph, const arma::mat & cumadopt, int wtype = 0,
    double v = 1.0, bool undirected=true, bool normalized=true);

arma::mat cumulative_adopt_count_cpp(const arma::mat & cumadopt);

arma::rowvec hazard_rate_cpp(const arma::mat & cumadopt);

arma::colvec threshold_cpp(const arma::mat & exposure,
                           const arma::vec & toe);

#endif
