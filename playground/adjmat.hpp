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

#ifndef ADJMAT_H
#define ADJMAT_H

using namespace Rcpp;

List adopt_mat_cpp(const IntegerVector & year);

arma::mat edgelist_to_adjmat_cpp(
    const arma::mat & data, NumericVector weights = NumericVector::create(),
    int n = 0,bool undirected = false);

arma::mat adjmat_to_edgelist_cpp(
    const arma::mat & adjmat, bool undirected = true);

IntegerMatrix toa_mat_cpp(const IntegerVector & year);

IntegerVector isolated_cpp(const arma::mat & adjmat);

arma::mat drop_isolated_cpp(const arma::mat & adjmat, bool undirected=true) {
  int n = adjmat.n_cols;
  arma::colvec isolated = isolated_cpp(adjmat, undirected);

arma::colvec degree_cpp(const arma::mat & adjmat, const int & cmode=2,
                           bool undirected=true, bool self=false);

arma::mat rand_graph_cpp(int n=10, double p = 0.3, bool undirected=true,
                         bool weighted=false, bool self=false);

arma::cube rand_dyn_graph_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false);

#endif
