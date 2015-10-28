#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Combine two vectors and (if) get the unique vect
NumericVector vec_comb(const NumericVector & a, const NumericVector & b, bool uni = false) {

  // Parameters
  int n_a = a.size(), n_b = b.size();

  NumericVector c(n_a+n_b);

  for(int i=0;i<n_a;i++)
    c[i] = a[i];

  for(int i=0;i<n_b;i++)
    c[i+n_a] = b[i];

  // If returning a unique vector
  if (uni) return(unique(c));
  return c;
}

// [[Rcpp::export]]
List adopt_mat_cpp(const IntegerVector & year) {

  // Measuring time
  int T0 = min(year);
  int T = max(year) - T0 + 1;
  int n = year.size();

  // Creating output
  List out(2);
  IntegerMatrix adoptmat1(n,T);

  for(int i=0;i<n;i++)
    adoptmat1(i,year[i]-T0) = 1;

  IntegerMatrix adoptmat = clone(adoptmat1);
  for(int i=0;i<n;i++)
    for(int j=1;j<T;j++)
      adoptmat(i,j) += adoptmat(i,j-1);

  return List::create(_["adoptmat"]=adoptmat,_["adoptmat1"]=adoptmat1);
}

// [[Rcpp::export]]
arma::mat edgelist_to_adjmat_cpp(
    const NumericMatrix & edgelist,
    NumericVector weights = NumericVector::create(),
    int n = 0,
    bool undirected = false) {

  // Getting the dimensions and creating objects
  int m = edgelist.nrow();

  // Checking out weights
  NumericVector w(m,1.0);
  if (weights.size() == m) w = clone(weights);

  // Identifying the unique nodes and creating the adjmat (base)
  if (n == 0) {
    NumericVector nodes = vec_comb(edgelist(_,0), edgelist(_,1), true);
    n = max(nodes);
  }
  arma::mat mat(n,n, arma::fill::zeros);

  for(int i=0;i<m;i++) {
    mat(edgelist(i,0)-1,edgelist(i,1)-1) += w[i];
    if (undirected) mat(edgelist(i,1)-1,edgelist(i,0)-1) += w[i];
  }

  return mat;
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
IntegerVector isolated_cpp(const NumericMatrix & adjmat) {

  int n = adjmat.ncol();
  IntegerVector isolated(n,1);

  // Looping through (all) the matrix. Setting the value to 0
  // whenever there's a link (for both individuals).
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      if (adjmat(i,j))
        isolated[i] = 0, isolated[j] = 0;

  return isolated;
}

/* cmode:
 *  0: Indegree
 *  1: Outdegree
 *  2: Degree
 */
// [[Rcpp::export]]
arma::colvec degree_cpp(
    const NumericMatrix & adjmat, const int & cmode=2,
    bool undirected=true, bool self=false) {

  int n = adjmat.ncol();
  NumericVector indegree(n);
  NumericVector oudegree(n);
  arma::colvec degree(n);

  for(int i=0;i<n;i++) {
    int m=n;
    if (undirected) m = i;
    for(int j=0;j<m;j++) {

      // Checking out whether compute self or not
      if (!self && i==j) continue;

      double val = adjmat(i,j);
      if (val!=0.0) {
        if (cmode!=1) indegree[j] += val;
        if (cmode!=0) oudegree[i] += val;
      }
    }
  }

  // Adding all up
  for(int i=0;i<n;i++)
    degree[i]+=indegree[i]+oudegree[i];

  return degree;
}

/* **R
edgelist <- rbind(c(2,1),c(3,1),c(3,2))
adjmat <- edgelist_to_adjmat_cpp(edgelist)
degree_cpp(adjmat,0,FALSE)
degree_cpp(adjmat,1,FALSE)
degree_cpp(adjmat,2,FALSE)
*/

// [[Rcpp::export]]
arma::mat rand_graph_cpp(
    int n=10, double p = 0.5,
    bool weighted=false, bool undirected=true, bool self=false) {
  arma::mat graph(n, n, arma::fill::zeros);

  NumericVector datasource = runif(n*n);

  for(int i=0;i<n;i++) {
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;j++) {
      if (!self && (i==j)) continue;
      double val = datasource[i*n+j];
      if (val > (1-p)) {
        graph(i,j) = 1.0;
        if (undirected) graph(j,i) = 1.0;
      }
    }
  }

  return graph;
}

// [[Rcpp::export]]
arma::cube mycube_cpp(int n, int t) {
  arma::cube x(n,n,t, arma::fill::zeros);
  return x;
}

/***R
set.seed(123)
rand_graph_cpp()

set.seed(123)
rand_graph_cpp()

set.seed(123)
rand_graph_cpp(undirected=FALSE)

library(igraph)
x <- rand_graph_cpp(10)
z <- graph_from_adjacency_matrix(x);plot(z)

x <- rand_graph_cpp(10000, p=.1)
sum(x)/length(x)
*/
