/*******************************************************************************
 * struct_test.h header function for struct_test.cpp
 *******************************************************************************/
#include <RcppArmadillo.h>
using namespace Rcpp;


#ifndef NETDIFFUSER_STRUCT_TEST_
#define NETDIFFUSER_STRUCT_TEST_

double struct_test_mean(NumericVector & y, std::string funname, bool self=false);
double struct_test_var(NumericVector & y, std::string funname, bool self=false);

#endif
