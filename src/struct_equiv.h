/*******************************************************************************
* adjmat.h
*
* List of functions for creating and manipulating adjacency matrices. Here are
* included:
*  - vec_comb: combining vectors and extracting unique elements
*  - adopt_mat_cpp: Based on moment of adoption, creates a matrix indicanting
*    when it was adopted (n x T)
*  - edgelist_to_adjmat_cpp: Converts an edgelist to a adjmat
*  - toa_mat_cpp: creates the time of adoption matrix
*
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef DIFFTEST_STRUCT_EQUIV_
#define DIFFTEST_STRUCT_EQUIV_

using namespace Rcpp;

List struct_equiv_cpp(const arma::sp_mat & gdist,
                      double v = 1.0,
                      bool unscaled = false,
                      bool inv = false, double invrep = 0.0);

#endif
