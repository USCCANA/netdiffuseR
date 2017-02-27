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

arma::sp_mat sp_as_undirected(const arma::sp_mat & x);

typedef double (*funcPtr)(double y0, double y1);

double st_dist(double y0, double y1);
double st_quaddist(double y0, double y1);
double st_greater(double y0, double y1);
double st_greaterequal(double y0, double y1);
double st_smaller(double y0, double y1);
double st_smallerequal(double y0, double y1);
double st_equal(double y0, double y1);
double st_min(double y0, double y1);
double st_max(double y0, double y1);
double st_mean(double y0, double y1);

void st_getfun(std::string funname, funcPtr & fun);

arma::sp_mat bootnet_fillself(
    arma::sp_mat & graph, const IntegerVector & index, const NumericVector & E
  );

#endif
