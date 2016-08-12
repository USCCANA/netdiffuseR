/*******************************************************************************
 * struct_test.h header function for struct_test.cpp
 *******************************************************************************/
#include <RcppArmadillo.h>
using namespace Rcpp;


#ifndef NETDIFFUSER_STRUCT_TEST_
#define NETDIFFUSER_STRUCT_TEST_

typedef double (*funcPtr)(double y0, double y1);

double st_dist(double y0, double y1);
double st_greater(double y0, double y1);
double st_greaterequal(double y0, double y1);
double st_smaller(double y0, double y1);
double st_smallerequal(double y0, double y1);
double st_equal(double y0, double y1);

void st_getfun(std::string funname, funcPtr & fun);

double struct_test_mean(NumericVector & y, std::string funname, bool self=false);
double struct_test_var(NumericVector & y, std::string funname, bool self=false);

#endif
