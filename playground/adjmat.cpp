#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

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
  arma::mat mat(n,n, arma::fill::zeros);

  for(int i=0;i<m;i++) {
    mat(edgelist(i,0)-1,edgelist(i,1)-1) += w[i];
    if (undirected) mat(edgelist(i,1)-1,edgelist(i,0)-1) += w[i];
  }

  return mat;
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

  arma::mat output = join_rows(egom, alterm);

  return output;
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
IntegerVector isolated_cpp(const arma::mat & adjmat) {

  int n = adjmat.n_cols;
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
    const arma::mat & adjmat, const int & cmode=2,
    bool undirected=true, bool self=false) {

  int n = adjmat.n_cols;
  arma::colvec indegree(n, arma::fill::zeros);
  arma::colvec oudegree(n, arma::fill::zeros);

  for(int i=0;i<n;i++) {
    int m=n;
    if (undirected) m = i;
    for(int j=0;j<m;j++) {

      // Checking out whether compute self or not
      if (!self && i==j) continue;

      double val = adjmat(i,j);
      if (val!=0.0) {
        if (cmode!=1) indegree(j) += val;
        if (cmode!=0) oudegree(i) += val;
      }
    }
  }

  arma::colvec degree = indegree+oudegree;

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
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {
  arma::mat graph(n, n, arma::fill::zeros);

  NumericVector datasource = runif(n*n);

  double w = 0.0;
  for(int i=0;i<n;i++) {

    /* Setting the length of the subloop acordingly to type of graph */
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;j++) {

      /* Assessing if include self */
      if (!self && (i==j)) continue;

      /* Setting the value of the tie */
      double val = datasource[i*n+j];
      w = val;

      if (val > (1-p)) {
        if (!weighted) w=1.0;
        graph(i,j) = w;
        if (undirected) graph(j,i) = w;
      }
    }
  }

  return graph;
}

// [[Rcpp::export]]
arma::cube rand_dyn_graph_cpp(
    int n=10, int t=3, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {

  arma::cube graphs(n,n,t);
  for(int i=0;i<t;i++)
    graphs.slice(i) = rand_graph_cpp(n, p, undirected, weighted, self);

  return graphs;

}

// [[Rcpp::export]]
arma::cube mycube_cpp(int n, int t) {
  arma::cube x(n,n,t, arma::fill::zeros);
  return x;
}

/* **R
set.seed(123)
rand_graph_cpp()

rand_graph(undirected=FALSE)

x <- rand_graph(1000, p=.1)
sum(x)/length(x)
*/
