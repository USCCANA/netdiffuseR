#include <Rcpp.h>
using namespace Rcpp;

// Combine two vectors and (if) get the unique vect
NumericVector vec_comb(NumericVector a, NumericVector b, bool uni = false) {

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
NumericMatrix edgelist_to_adjmat_cpp(
    const NumericMatrix data,
    bool undirected = false) {

  // Getting the dimensions and creating objects
  int m = data.nrow();
  int k = data.ncol();

  // Checking out weights
  NumericVector w(m,1.0);
  if (k>2) w = data(_,2);

  // Identifying the unique nodes and creating the adjmat (base)
  NumericVector nodes = vec_comb(data(_,0), data(_,1), true);
  int n = nodes.size();
  NumericMatrix mat(n,n);

  for(int i=0;i<m;i++) {
                    mat(data(i,0)-1,data(i,1)-1) += w[i];
    if (undirected) mat(data(i,1)-1,data(i,0)-1) += w[i];
  }

  return mat;
}

// [[Rcpp::export]]
NumericVector uniquecpp(NumericVector x) {return unique(x);}

/***R
x <- c(1,1,2,3,4,4,5,1)
uniquecpp(x)

# Adj mat
#' Recodes an edgelist such that ids go from 1 to n
#' @param data Edgelist as either a matrix or dataframe with ego and alter
#' @param ... Further arguments for the method (ignored)
#' @return A recoded edgelist
#' @details Recomended for ease of use
recode <- function(...) UseMethod("recode")

#' @describeIn recode Method for data.frame
recode.data.frame <- function(data, ...) {
  cn <- colnames(data)
  data <- as.data.frame(recode.matrix(as.matrix(data), ...))
  colnames(data) <- cn
  data
}

#' @describeIn recode method for matrix
recode.matrix <- function(data, ...) {
  data <- as.factor(as.vector(data))
  n <- length(data)
  cbind(data[1:(n/2)], data[(n/2+1):n])
}

edgelist <- cbind(c(1,1,3,6),c(4,3,200,1))
edgelist
edgelist <- recode(edgelist)
edgelist

edgelist_to_adjmat_cpp(edgelist)
edgelist_to_adjmat_cpp(edgelist, undirected = TRUE)

set.seed(123)
edgelist <- cbind(edgelist, abs(rnorm(nrow(edgelist))))

edgelist_to_adjmat_cpp(edgelist)
edgelist_to_adjmat_cpp(edgelist, undirected = TRUE)

edgelist
*/
