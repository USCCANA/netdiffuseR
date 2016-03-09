/*******************************************************************************
* infection.h header function for infection.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_INFECTSUSCEPT_
#define NETDIFFUSER_INFECTSUSCEPT_

using namespace Rcpp;

arma::mat infection_cpp(
    List graph, const arma::colvec & times, bool normalize = true,
    int K = 1, double r = 0.5, bool expdiscount = false, int n=0, int T=0,
    bool valued=false, bool outgoing=true);

arma::mat susceptibility_cpp(
    List graph, const arma::colvec & times, bool normalize = true,
    int K = 1, double r = 0.5, bool expdiscount = false, int n=0, int T=0,
    bool valued=false, bool outgoing=true);

#endif
