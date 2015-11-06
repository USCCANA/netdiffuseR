/*******************************************************************************
 adjmat.cpp: FUNCTIONS FOR NETWORK DATA MANAGEMENT

 created: nov 2, 2015
 copyright: MIT + LICENCE file

 This set of functions helps handling network data, in particular, diffussion
 network data. Importing edgelist as adjacency matrices and generating auxiliary
 matrices used for the model.

 This file contains the following functions:
   - adopt_mat_cpp: Creates an adoption matrix of size nxT
   - edgelist_to_adjmat_cpp: Creates an adjacency matrix from an edgelist
   - adjmat_to_edgelist_cpp: The converse of the previous function
   - toa_mat_cpp: Creates a Time of Adoption Matrix of size nxT
   - isolated_cpp: Identifies the isolated nodes in a network
   - drop_isolated_cpp: Removes isolated networks from an adjmat
*******************************************************************************/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "rand_graph.h"

using namespace Rcpp;

// [[Rcpp::export]]
List adopt_mat_cpp(const IntegerVector & year) {

  // Measuring time
  int T0 = min(year);
  int T = max(year) - T0 + 1;
  int n = year.size();

  // Creating output
  List out(2);
  arma::mat adopt(n,T,arma::fill::zeros);

  for(int i=0;i<n;i++)
    adopt(i,year[i]-T0) = 1.0;

  arma::mat cumadopt = cumsum(adopt, 1);

  /* Adopt_mat -> cumadopt
    Adopt_mat1 -> adopt
  */
  return List::create(_["adopt"]=adopt, _["cumadopt"]=cumadopt);
}

// [[Rcpp::export]]
arma::mat edgelist_to_adjmat_cpp(
    const arma::mat & edgelist,
    NumericVector weights = NumericVector::create(),
    int n = 0,
    bool undirected = false) {

  // Getting the dimensions and creating objects
  int m = edgelist.n_rows;

  // Checking out weights
  NumericVector w(m,1.0);
  if (weights.size() == m) w = clone(weights);

  // Identifying the unique nodes and creating the adjmat (base)
  if (n == 0) {
    n = (int) max(max(edgelist));
  }
  arma::mat adjmat(n,n, arma::fill::zeros);

  for(int i=0;i<m;i++) {
    adjmat(edgelist(i,0)-1,edgelist(i,1)-1) += w[i];
    if (undirected) adjmat(edgelist(i,1)-1,edgelist(i,0)-1) += w[i];
  }

  return adjmat;
}


// [[Rcpp::export]]
arma::mat adjmat_to_edgelist_cpp(const arma::mat & adjmat, bool undirected = true) {
  std::vector< double > ego;
  std::vector< double > alter;

  int n = adjmat.n_cols;

  for(int i=0;i<n;i++) {
    /* Setting the length of the subloop acordingly to type of graph */
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;j++) {
      if (adjmat(i,j))
        ego.push_back(i+1.0), alter.push_back(j+1.0);
    }
  }

  // Creating colvectors to be used with join_rows.
  arma::mat egom(ego);
  arma::mat alterm(alter);

  arma::mat edgelist = join_rows(egom, alterm);

  return edgelist;
}

// [[Rcpp::export]]
arma::mat adjmat_to_dyn_edgelist_cpp(NumericVector adjmat, bool undirected=true) {

  // Coersing a NumericVector into a cube for ease of use
  IntegerVector dims=adjmat.attr("dim");
  const arma::cube adjmat_cube(adjmat.begin(), dims[0], dims[1], dims[2], false);

  int T = adjmat_cube.n_slices;
  int n = adjmat_cube.n_cols;

  std::vector< double > ego;
  std::vector< double > alter;
  std::vector< double > time;

  for(int t=0;t<T;t++)
    for(int i=0;i<n;i++) {
      /* Setting the length of the subloop acordingly to type of graph */
      int m = n;
      if (undirected) m=i;
      for(int j=0;j<m;j++)
        if (adjmat_cube(i,j,t))
          ego.push_back(i+1.0), alter.push_back(j+1.0), time.push_back(t+1.0);
    }

  // Creating colvectors to be used with join_rows.
  arma::mat egom(ego);
  arma::mat alterm(alter);
  arma::mat timem(time);

  arma::mat edgelist = join_rows(join_rows(egom, alterm), timem);

  return edgelist;
}


// [[Rcpp::export]]
IntegerMatrix toa_mat_cpp(const IntegerVector & year) {
  int n = year.size();

  IntegerMatrix toa(n,n);

  for(int i=0;i<n;i++)
    for(int j=0;j<i;j++)
      toa(i,j) = year[j]-year[i], toa(j,i)=year[i]-year[j];

  return toa;
}

// [[Rcpp::export]]
arma::colvec isolated_cpp(
    const arma::mat & adjmat,
    bool undirected=true) {

  int n = adjmat.n_cols;
  arma::colvec isolated(n,arma::fill::ones);

  // Looping through (all) the matrix. Setting the value to 0
  // whenever there's a link (for both individuals).

  for(int i=0;i<n;i++) {
    // First, need to get the size of the loop
    int m=n;
    if (undirected) m = i;
    for(int j=0;j<m;j++)
      if (adjmat(i,j))
        isolated(i) = 0, isolated(j) = 0;
  }

  return isolated;
}

// [[Rcpp::export]]
arma::mat drop_isolated_cpp(const arma::mat & adjmat, arma::colvec isolated, bool undirected=true) {
  int n = adjmat.n_cols;
  if (isolated.n_rows==0) isolated = isolated_cpp(adjmat, undirected);

  int m = sum(isolated);
  arma::mat newadjmat(n-m,n-m);

  // Rcpp::warning()

  // Indexes of the new adjacency matrix
  int ii=0;
  int ji=0;
  for(int i=0;i<n;i++) {

    // If an isolated was found, continue next
    if (isolated(i)) continue;
    for(int j=0;j<n;j++) {
      if (isolated(j)) continue;
      newadjmat(ii,ji++)=adjmat(i,j);
    }
    // Continue next
    ii++,ji=0;
  }

  return newadjmat;
}

/** *R
set.seed(123)
lonenet <- rand_graph_cpp(8, p=.8,undirected=TRUE)
lonenet[c(1,4),] <- 0
lonenet[,c(1,4)] <- 0
isolated_cpp(lonenet)
*/
