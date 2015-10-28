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

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef ADJMAT_H
#define ADJMAT_H

using namespace Rcpp;

NumericVector vec_comb(const NumericVector & a, const NumericVector & b, bool uni = false);

List adopt_mat_cpp(const IntegerVector & year);

arma::matrix edgelist_to_adjmat_cpp(
    const NumericMatrix & data, NumericVector weights = NumericVector::create(),
    int n = 0,bool undirected = false);

IntegerMatrix toa_mat_cpp(const IntegerVector & year);

IntegerVector isolated_cpp(const NumericMatrix & adjmat);

arma::colvec degree_cpp(const NumericMatrix & adjmat, const int & cmode=2,
                           bool undirected=true, bool self=false);

arma::mat rand_graph_cpp(int n=10, bool weighted=false, bool undirected=true);

#endif
