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
#include <RcppArmadillo.h>

#ifndef STRUCT_EQUIV_H
#define STRUCT_EQUIV_H

using namespace Rcpp;

List struct_equiv_cpp(const arma::mat & gdist,double v = 1.0,
                      bool donly = false, bool inv = false, double invrep = 0.0);

#endif
