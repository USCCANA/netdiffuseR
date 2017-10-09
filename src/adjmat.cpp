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
    unsigned int ego   = ((unsigned int) edgelist(i,0)) - 1u;
    unsigned int alter = ((unsigned int) edgelist(i,1)) - 1u;

    // If self edges are not allowed
    if (!self && (ego == alter))
      continue;

    // If undirected, the order does not matters
    if (undirected && (ego > alter))
      std::swap(alter, ego);

    // If multiple edges are allowed.
    if (multiple | (adjmat(ego, alter) == 0.0))
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

  unsigned int m = adjmat.n_nonzero, i = 0u;
  arma::mat edgelist(m, 3);

  for (arma::sp_mat::const_iterator it = adjmat.begin(); it != adjmat.end(); it++) {
    edgelist.at(i,0) = it.row() + 1;
    edgelist.at(i,1) = it.col() + 1;
    edgelist.at(i++,2) = (*it);
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
    bool outer=true, bool self=true, bool valued=true) {

  // Creating containers
  std::vector< std::vector< unsigned int > > id(graph.n_rows);
  std::vector< std::vector< double > > value(graph.n_rows);

  typedef arma::sp_mat::const_iterator spiter;

  // We will fill it with zeros
  arma::uvec useit(graph.n_rows, arma::fill::zeros);
  for (unsigned int i = 0u; i < V.size(); i++)
    useit.at(V.at(i)) = 1u;

  // Finding values
  unsigned int i, j;
  for (spiter it = graph.begin(); it != graph.end(); it ++) {

    // Depending on outer or not
    if (outer) i = it.row(), j = it.col();
    else j = it.row(), i = it.col();

    if (useit.at(i) == 0u)
      continue;

    // Adding coordinates
    id.at(i).push_back(j + 1u);
    value.at(i).push_back( (valued)? (*it) : 1.0 );
  }

  // Coercing output
  List data(V.size());

  for (i = 0u; i < V.size(); i++) {

    // If self
    if (self)
      id.at(V.at(i)).push_back( V.at(i) + 1u),
      value.at(V.at(i)).push_back( graph.at(V.at(i), V.at(i)));

    // Creating the dataframe
    data.at(i) = DataFrame::create(
      _["id"]    = Rcpp::wrap(id.at(V.at(i))),
      _["value"] = Rcpp::wrap(value.at(V.at(i)))
      );
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


//[[Rcpp::export]]
arma::sp_mat approx_geodesicCpp(
    const arma::sp_mat & G,
    unsigned int n = 6,
    bool warn = false
) {

  int N = (int) G.n_cols;
  arma::sp_mat ans(N,N);

  typedef arma::sp_mat::const_iterator spiter;

  // Going through the steps
  arma::sp_mat pG = G;
  arma::sp_mat G0 = G;
  int change_count = 0;
  n++;
  unsigned int nsteps;
  for (unsigned int i=1u; i<n; i++) {

    nsteps = 0u;

    // Computing nsteps
    std::vector< unsigned int > irow;
    std::vector< unsigned int > icol;
    std::vector< unsigned int > ival;

    // Iterating throught the power graph's elements
    for (spiter it = pG.begin(); it != pG.end(); it ++) {

      if (it.row() == it.col())
        continue;

      // Checking user interrupt
      if (++nsteps % 1000)
        Rcpp::checkUserInterrupt();

      if (ans.at(it.row(), it.col()) == 0u) {
        // Storing coordinates
        irow.push_back(it.row());
        icol.push_back(it.col());

        ival.push_back(i);
        ++change_count;
      }
    }

    // Refilling the ans using batch insertion
    arma::sp_mat tmp(
        arma::join_cols(
          arma::conv_to< arma::urowvec >::from(irow),
          arma::conv_to< arma::urowvec >::from(icol)
        ),
        arma::conv_to< arma::colvec >::from(ival),
        N, N,
        true,
        false
      );

    // Adding the delta
    ans = ans + tmp;

    // Was there any change?
    if (!change_count) {
      if (warn)
        warning("The algorithm stopped at %i iterations.", i);
      break;
    } else change_count = 0;

    // Graph power
    pG *= G0;
  }

  return ans;
}


/***R

library(sna)
  library(netdiffuseR)
  library(igraph)
  library(microbenchmark)

  set.seed(123)
  g_sp  <- rgraph_ws(n=500, k = 3, p = .2, self = FALSE)
  g_mat <- as.matrix(g_sp)
  g_ig  <- graph_from_adjacency_matrix(g_sp)

  microbenchmark(
    sna = geodist(g_mat),
    nd  = approx_geodesic(g_sp, 6),
    ig  = distances(g_ig),
    times = 500,
    unit = "ms"
  )

  # Unit: milliseconds
  # expr       min        lq     mean    median       uq      max neval
  # sna 32.578531 37.298168 43.86050 40.834380 46.49778 196.5955   500
  # nd  7.959418  8.679596 10.24444  9.251891 10.61416  22.9517   500
  # ig 11.516345 13.012492 15.67670 14.049763 16.31242 171.5637   500

  ans0 <- geodist(g_mat)[[2]]
ans1 <- as.matrix(approx_geodesic(g_sp, 10))

  are0 <- which(ans1[] != 0, arr.ind = TRUE)
  prop.table(table(ans0 - ans1))

  g_sp  <- kfamilyDiffNet$graph$`1`

  microbenchmark(
    sna = geodist(as.matrix(g_sp)),
    nd  = approx_geodesic(g_sp, 10), times = 30,
    unit = "ms"
  )

  ans0 <- geodist(as.matrix(g_sp))[[2]]
ans1 <- as.matrix(approx_geodesic(g_sp, 20))

  are0 <- which(ans1[] != 0, arr.ind = TRUE)
  prop.table(table(ans0 - ans1))




  */
