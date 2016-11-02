/*******************************************************************************
 adjmat.cpp: FUNCTIONS FOR NETWORK DATA MANAGEMENT

 created: nov 2, 2015
 copyright: MIT + LICENCE file

 This set of functions helps handling network data, in particular, diffussion
 network data. Importing edgelist as adjacency matrices and generating auxiliary
 matrices used for the model.

*******************************************************************************/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"

using namespace Rcpp;

// [[Rcpp::export]]
List toa_mat_cpp(const IntegerVector & year, int t0, int t1) {

  // Pin down NAs
  LogicalVector isna = is_na(year);

  // Measuring time
  int n = year.size();

  // Creating output
  List out(2);
  arma::mat adopt(n,t1 - t0 + 1, arma::fill::zeros);

  for(int i=0;i<n;i++) {
    if (isna[i]) continue;
    adopt(i,year[i]-t0) = 1.0;
  }

  arma::mat cumadopt = cumsum(adopt, 1);

  /* Adopt_mat -> cumadopt
    Adopt_mat1 -> adopt
  */
  return List::create(_["adopt"]=adopt, _["cumadopt"]=cumadopt);
}

// [[Rcpp::export]]
arma::sp_mat edgelist_to_adjmat_cpp(
    const arma::mat & edgelist,
    NumericVector weights = NumericVector::create(),
    int n = 0,
    bool undirected = false,
    bool self = false,
    bool multiple = false) {

  // Getting the dimensions and creating objects
  int m = edgelist.n_rows;

  // Checking out weights
  NumericVector w(m,1.0);
  if (weights.size() == m) w = clone(weights);

  // Identifying the unique nodes and creating the adjmat (base)
  // If n is not provided (0), then we need to calculate the number of
  // nodes in the network.
  if (n == 0) {
    n = (int) max(max(edgelist));
  }

  // Creating output mat and aux mat. The counts mat is used in logical expressions
  // since the arma::mat class can lead to error when compared to 0 (because of
  // the precision)
  arma::sp_mat adjmat(n,n);

  for(int i=0;i<m;i++) {
    // Ids of the vertices
    int ego   = ((int) edgelist(i,0)) - 1;
    int alter = ((int) edgelist(i,1)) - 1;

    // If self edges are not allowed
    if (!self && (ego == alter)) continue;

    // If undirected, the order does not matters
    if (undirected) {
      int tmp=ego;
      if (ego > alter) ego = alter, alter=tmp;
    }

    // If multiple edges are not allowed.
    if (multiple | (adjmat(ego, alter) == 0))
      adjmat.at(ego, alter) += w[i];

    // If undirected, must include in the switch
    if (undirected)
      if (multiple | (adjmat(alter, ego) == 0))
        adjmat.at(alter, ego) += w[i];
  }

  return adjmat;
}


// [[Rcpp::export]]
arma::mat adjmat_to_edgelist_cpp(
    const arma::sp_mat & adjmat,
    bool undirected = true) {

  arma::umat coords = sparse_indexes(adjmat);
  int m = coords.n_rows;
  arma::mat edgelist(m, 3);

  for (int i=0;i<m;i++) {
    edgelist.at(i,0) = coords.at(i, 0) + 1;
    edgelist.at(i,1) = coords.at(i, 1) + 1;
    edgelist.at(i,2) = adjmat.at(coords.at(i,0), coords.at(i,1));
  }

  return edgelist;
}

// [[Rcpp::export]]
IntegerMatrix toa_diff_cpp(const IntegerVector & year) {
  int n = year.size();

  IntegerMatrix diff(n,n);
  LogicalVector isna = is_na(year);

  for(int i=0;i<n;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // If na, then fill the diff with NA
    if (isna[i]) {
      for(int j=0;j<n;j++)
        diff(i,j) = NA_INTEGER, diff(j,i) = NA_INTEGER;
      continue;
    }
    for(int j=0;j<i;j++) {
      if (isna[j]) continue;
      diff(i,j) = year[j]-year[i], diff(j,i)=year[i]-year[j];
    }
  }

  return diff;
}

// [[Rcpp::export]]
arma::icolvec isolated_cpp(
    const arma::sp_mat & adjmat,
    bool undirected=true) {

  int n = adjmat.n_cols;

  // Checking size
  if ((int) adjmat.n_rows != n)
    stop("It must be an square matrix.");

  arma::icolvec isolated(n, arma::fill::ones);

  // Looping through (all) the matrix. Setting the value to 0
  // whenever there's a link (for both individuals).

  for(int i=0;i<n;i++) {
    // First, need to get the size of the loop
    int m=n;
    if (undirected) m = i;
    for(int j=0;j<m;j++)
      if ((adjmat.at(i,j)!=0))
        isolated[i]=0, isolated[j]=0;
  }

  return isolated;
}

// [[Rcpp::export]]
arma::sp_mat drop_isolated_cpp(
    const arma::sp_mat & adjmat,
    arma::icolvec isolated, bool undirected=true) {

  // Checking size
  if (adjmat.n_rows != adjmat.n_cols)
    stop("It must be an square matrix.");

  if (isolated.n_rows==0u) isolated = isolated_cpp(adjmat, undirected);
  else if (isolated.n_rows != adjmat.n_cols) stop("-isolated- must have the same length as nrow(adjmat)");

  int m = sum(isolated);
  arma::sp_mat newadjmat(adjmat.n_cols-m,adjmat.n_cols-m);

  // Rcpp::warning()

  // Indexes of the new adjacency matrix
  int ii=0;
  int ji=0;
  for(unsigned i=0;i<adjmat.n_cols;i++) {

    // If an isolated was found, continue next
    if (isolated[i]) continue;
    for(unsigned j=0;j<adjmat.n_cols;j++) {
      if (isolated[j]) continue;
      newadjmat.at(ii,ji++)=adjmat.at(i,j);
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

// [[Rcpp::export]]
List egonet_attrs_cpp(
    const arma::sp_mat & graph, const arma::uvec V,
    bool outer=true, bool self=true, bool self_attrs=false, bool valued=true) {

  // General variables
  int N = V.n_elem;

  // Column names
  CharacterVector cnames(2);

  cnames[0] = "value";
  cnames[1] = "id";

  // Depending on inner or outer edges:
  //  - outer: Accounts for the rows
  //  - inner: Accounts for the columns
  arma::sp_mat tgraph;
  if (outer) tgraph = graph.t();
  else tgraph = graph;

  List data(N);

  for (int v=0;v<N;v++) {
    // Index
    int e = V.at(v);

    // Analyzing the case when is undirected
    int rm = 0;
    if (graph.at(e,e) && !self) rm = 1;

    // Individual specific variables
    arma::sp_mat g = tgraph.col(e);
    NumericMatrix out(g.n_nonzero-rm + (self_attrs ? 1 : 0),2);
    // We add two so:
    //  - we can include the value of the edge
    //  - we cna include the id number of the vertex (wich goes from 1 to n)

    // Retrieving the desired set of attributes
    if (self_attrs) {
      out(0,0) = tgraph.at(v,v);
      out(0,1) = v + 1;
    }
    int nloop = (self_attrs ? 1 : 0);
    for (unsigned i = 0u;i < g.n_nonzero;i++) {

      // Edge index
      int index = g.row_indices[i];
      if (!self && index == e) continue;
      out(nloop,1) = index + 1;

      // Edge value
      if (valued) out(nloop,0) = g.values[i];
      else out(nloop,0) = 1.0;

      // Increasing after success
      nloop++;
    }

    colnames(out) = cnames;

    data[v] = out;
  }


  return data;
}

/***R
library(netdiffuseR)
set.seed(123)
n <- 1e3
graph <- rgraph_ba(t=n-1)
mat <- cbind(id = 1:n, rand=runif(n), d=dgr(graph))

diffnet_attrs <- function(diffnet) {
  attrs
}



onestep_attrs <- egonet_attrs_cpp(graph,1:n-1, self = FALSE)
# twostep_attrs <- nsteps_attrs_cpp(graph %*% graph,1:n, mat, self = FALSE)
stop()
# Computing the degree weighted rand
w <- lapply(onestep_attrs, function(x) {
  (t(x[,4]) %*% x[,5])/sum(x[,5])
})

w <- do.call(rbind, w)
*/

/***R
library(sna)

g <- matrix(0, ncol=6, nrow=6)
g[2:3,1] <- g[1,2:3] <- 1
g[4,2] <- g[2,4] <-1
g[6,5] <- g[5,6] <- 1

oldpar <- par(no.readonly = TRUE)
coord <- gplot.layout.fruchtermanreingold(g, list())
par(mfrow=c(2,2))
gplot(g, displaylabels = TRUE, coord = coord, edge.lwd = 1)
gplot(g %*% g, displaylabels = TRUE, coord = coord, edge.lwd = 1)
gplot(g %*% g %*% g, displaylabels = TRUE, coord = coord, edge.lwd = 2)
gplot(g %*% g %*% g %*% g, displaylabels = TRUE, coord = coord, edge.lwd = 1)
par(oldpar)

*/
