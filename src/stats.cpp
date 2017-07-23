/*******************************************************************************
 stats.cpp: FUNCTIONS FOR CALCULATING NETWORK STATISTICS

 created: nov 2, 2015
 copyright: MIT + LICENCE file

 This set of functions helps handling network data, in particular, diffussion
 network data. Importing edgelist as adjacency matrices and generating auxiliary
 matrices used for the model.

 This file contains the following functions:
  - toa_mat_cpp: Creates an adoption matrix of size nxT
  - edgelist_to_adjmat_cpp: Creates an adjacency matrix from an edgelist
  - adjmat_to_edgelist_cpp: The converse of the previous function
  - toa_mat_cpp: Creates a Time of Adoption Matrix of size nxT
  - isolated_cpp: Identifies the isolated nodes in a network
  - drop_isolated_cpp: Removes isolated networks from an adjmat
*******************************************************************************/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"

using namespace Rcpp;



//' @export
//' @rdname vertex_covariate_dist
// [[Rcpp::export]]
arma::sp_mat vertex_covariate_dist(const arma::sp_mat & graph,
                                   const arma::mat & X,
                                   double p = 2.0) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  // Iterating over elements of graph
  int i,j;
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++) {

    i = iter.row();
    j = iter.col();

    ans.at(i,j) = pow(
      arma::as_scalar((X.row(i) - X.row(j)) * (X.row(i) - X.row(j)).t()),
      1/p
    );
  }

  return ans;
}

// [[Rcpp::export]]
arma::sp_mat vertex_mahalanobis_dist_cpp(const arma::sp_mat & graph,
                                   const arma::mat & X,
                                   const arma::mat & S) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  arma::mat Sinv = S.i();

  // Iterating over elements of graph
  int i,j;
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++) {

    i = iter.row();
    j = iter.col();

    ans.at(i,j) = pow(
      arma::as_scalar((X.row(i) - X.row(j)) * Sinv * (X.row(i) - X.row(j)).t()),
      0.5
    );
  }

  return ans;
}

//' Comparisons at dyadic level
//' @param graph A matrix of size \eqn{n\times n}{n*n} of class \code{dgCMatrix}.
//' @param X A numeric vector of length \eqn{n}.
//' @param funname Character scalar. Comparison to make (see details).
//' @details
//'
//' This auxiliary function takes advantage of the sparseness of \code{graph} and
//' applies a function in the form of \eqn{funname(x_i,x_j)}{funname(X[i],X[j])}
//' only to \eqn{(i,j)} that have no empty entry. In other words, applies a compares
//' elements of \code{X} only between vertices that have a link; making
//' \code{nlinks(graph)} comparisons instead of looping through \eqn{n\times n}{n*n},
//' which is much faster.
//'
//' \code{funname} can take any of the following values:
//' \code{"distance"}, \code{"^2"} or \code{"quaddistance"}, \code{">"} or \code{"greater"},
//' \code{"<"} or \code{"smaller"}, \code{">="} or \code{"greaterequal"},
//' \code{"<="} or \code{"smallerequal"}, \code{"=="} or \code{"equal"}.
//' @return A matrix \code{dgCMatrix} of size \eqn{n\times n}{n*n} with values in
//' the form of \eqn{funname(x_i,x_j)}{funname(X[i],X[j])}.
//' @examples
//'
//' # Basic example ------------------------------------------------------------
//' set.seed(1313)
//' G <- rgraph_ws(10, 4, .2)
//' x <- rnorm(10)
//'
//' vertex_covariate_compare(G, x, "distance")
//' vertex_covariate_compare(G, x, "^2")
//' vertex_covariate_compare(G, x, ">=")
//' vertex_covariate_compare(G, x, "<=")
//' @export
// [[Rcpp::export]]
arma::sp_mat vertex_covariate_compare(const arma::sp_mat & graph, const NumericVector & X,
                                      std::string funname) {

  // Creating objects
  arma::sp_mat ans(graph.n_rows,graph.n_cols);
  arma::sp_mat::const_iterator b = graph.begin();
  arma::sp_mat::const_iterator e = graph.end();

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  // Iterating over elements of graph
  for (arma::sp_mat::const_iterator iter = b; iter != e; iter++)
    ans.at(iter.row(),iter.col()) = fun(X[iter.row()], X[iter.col()]);

  return ans;
}

// [[Rcpp::export]]
List moran_cpp(const arma::colvec & x, const arma::sp_mat & w) {
  double xmean = mean(x);

  int N = x.n_rows;

  // Checking dims
  if (w.n_cols != w.n_rows) stop("-w- is not a square matrix");
  if (N != (int) w.n_cols) stop("-x- and -w- dimensions differ");

  // Weights sum
  double wsum = accu(w);

  arma::colvec xcent = x - xmean;
  double numer = accu((xcent * xcent.t()) % w);
  double denom = accu(pow(xcent, 2.0));

  // Computing standard error
  double s1 = 0.0, s2 = 0.0, s3, s4, s5;
  arma::sp_mat::const_iterator b;

  // S1
  arma::sp_mat wtw = w + w.t();
  for (b = wtw.begin(); b!= wtw.end(); b++)
    s1 += pow((*b), 2.0);

  s1/= 2.0;

  // S2

  // row sums
  arma::sp_colvec w_rowsums = sum(w, 1);
  arma::sp_rowvec w_colsums = sum(w, 0);

  for (unsigned int i = 0; i < N; i++)
    s2 += pow(w_rowsums.at(i) + w_colsums.at(i), 2.0);

  // S3
  s3 = accu(pow(xcent, 4.0))/N/
    pow(accu(pow(xcent, 2.0))/N, 2.0);

  // S4
  s4 = (N*N - 3.0 * N + 3.0)*s1 - N*s2 + 3*pow(wsum, 2.0);

  // S5
  s5 = N*(N - 1.0)*s1 - 2.0*N*s2 + 6.0*pow(wsum, 2.0);

  return Rcpp::List::create(
    _["observed"] = (N/wsum)*(numer/(denom + 1e-20)),
    _["expected"] = -1.0 / (N - 1.0),
    _["sd"] = sqrt(
      (N*s4 - s3*s5)/((N - 1.0)*(N - 2.0)*(N - 3.0) * pow(wsum, 2.0)) -
      pow(-1.0/(N - 1.0), 2.0)
    )
  );

}


// [[Rcpp::export]]
List struct_equiv_cpp(
    const arma::sp_mat & graph, // Must be a geodesic distances graph
    double v = 1.0,
    bool unscaled = false,
    bool inv = false, double invrep = 0.0) {

  int n = graph.n_cols;
  if (graph.n_cols != graph.n_rows) stop("-graph- is not square.");

  NumericMatrix d(n,n);

  // Calculating Z vector as Z_i - Z_j = {z_ik - z_jk}
  NumericVector dmax(n, -1e100);
  for(int i=0;i<n;i++) {
    for(int j=0;j<i;j++) {

      // Computing sum(z_ik - z_jk)
      double sumik = 0.0;
      double sumki = 0.0;
      for(int k=0;k<n;k++) {
        // Summation accross all but i and j
        if (k == i || k == j) continue;
        sumik += pow(graph.at(i,k)-graph.at(j,k), 2.0);
        sumki += pow(graph.at(k,i)-graph.at(k,j), 2.0);
      }

      // Adding up the results
      d.at(i,j) = pow(pow(graph.at(i,j) - graph.at(j,i), 2.0) + sumik + sumki, 0.5 );

      // // If only inverse required
      // if (inv && unscaled) d.at(i,j) = 1.0/(d.at(i,j) + 1e-15);

      d.at(j,i) = d.at(i,j);
    }
  }

  // // If only distance must be computed
  // if (unscaled) return List::create(_["SE"]=d, _["d"]=d, _["gdist"]=graph);

  // Computing distances
  NumericMatrix SE(n,n);

  for(int i=0;i<n;i++) {

    // Getting the max of the line
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      if (dmax[i] < d.at(i,j)) dmax[i] = d.at(i,j);
    }

    // Computing sum(dmax - dkj)
    double sumdmaxd = 0.0;
    for(int k=0;k<n;k++) {
      if (k==i) continue;
      sumdmaxd += pow(dmax[i] - d.at(k,i), v);
    }

    // Computing (dmax - d)/sum(dmax - d)
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      SE.at(i,j) = pow(dmax[i] - d.at(j,i), v)/(sumdmaxd + 1e-15);
    }

    // // If inverse required
    // if (inv) {
    //   for(int j=0;j<n;j++) {
    //     SE.at(i,j) = 1/(SE.at(i,j) + 1e-10);
    //   }
    // }
  }

  return List::create(_["SE"]=SE, _["d"]=d, _["gdist"]=graph);
}

// [[Rcpp::export]]
arma::sp_mat matrix_compareCpp(
    const arma::sp_mat & A,
    const arma::sp_mat & B,
    Function fun
) {

  // Checking dimmensions
  int n = A.n_cols;
  int m = A.n_rows;


  // Comparing
  typedef arma::sp_mat::const_iterator spiter;

  arma::sp_umat done(n, m);
  arma::sp_mat   ans(n, m);

  // Iterating through matrix A
  for (spiter iter=A.begin(); iter!= A.end(); iter++)
    // Filling the value
    ans.at(iter.row(), iter.col()) =
      as< double >(fun(*iter, B.at(iter.row(), iter.col()) )),
      done.at(iter.row(), iter.col()) = 1u;

  // Iterating throught matrix B
  for (spiter iter=B.begin(); iter!= B.end(); iter++)
    if (done.at(iter.row(), iter.col()) != 1u)
      // Filling the value
      ans.at(iter.row(), iter.col()) =
        as< double >(fun(A.at(iter.row(), iter.col()), *iter));


  return ans;

}

/** *R
 set.seed(1234)
 adjmat <- rand_graph_cpp()

 struct_equiv <- function(adjmat, v=1, ...) {
 geod <- sna::geodist(adjmat, inf.replace = 0, ...)
 geod[["gdist"]] <- geod[["gdist"]]/max(geod[["gdist"]])
 struct_equiv_cpp(geod[["gdist"]], v)
 }

 microbenchmark::microbenchmark(sna::sedist, struct_equiv, times=10000)

# Use this implementation
 new <- struct_equiv(adjmat)
 old <- sna::sedist(adjmat, method = 'euclidean')
 rowSums(new$SE)
 new
 old
 */
