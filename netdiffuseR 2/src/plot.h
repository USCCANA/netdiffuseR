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

arma::mat edges_arrow(
    const double & x0,
    const double & y0,
    const double & x1,
    const double & y1,
    const double & height,
    const double & width,
    const double beta = 1.5707963267949, // PI/2
    NumericVector dev = NumericVector::create(),
    NumericVector ran = NumericVector::create()
);

NumericMatrix edges_coords(
    const arma::sp_mat & graph,
    const arma::colvec & toa,
    const arma::colvec & x,
    const arma::colvec & y,
    const arma::colvec & vertex_cex,
    bool undirected=true,
    bool no_contemporary=true,
    NumericVector dev = NumericVector::create(),
    NumericVector ran = NumericVector::create()
);

List vertices_coords(
    const arma::colvec & x,
    const arma::colvec & y,
    const arma::colvec & size,
    const arma::colvec & nsides,
    const arma::colvec & rot,
    NumericVector dev = NumericVector::create(),
    NumericVector ran = NumericVector::create()
);
#endif
