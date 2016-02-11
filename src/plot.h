/*******************************************************************************
* plot.h header function for plot.cpp
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef NETDIFFUSER_PLOT_
#define NETDIFFUSER_PLOT_

using namespace Rcpp;

arma::vec seq_cpp(double from, double to, int lengthout);

List grid_distribution(const arma::vec & x, const arma::vec & y, int nlevels=100);

NumericMatrix edges_coords(
    const arma::sp_mat & graph,
    const arma::colvec & toa,
    const arma::colvec & x,
    const arma::colvec & y,
    const arma::colvec & vertex_cex,
    bool undirected=true,
    bool no_contemporary=true);
#endif
