#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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
    const NumericMatrix & data,
    NumericVector weights = NumericVector::create(),
    int n = 0,
    bool undirected = false) {

  // Getting the dimensions and creating objects
  int m = data.nrow();

  // Checking out weights
  NumericVector w(m,1.0);
  if (weights.size() == m) w = clone(weights);

  // Identifying the unique nodes and creating the adjmat (base)
  if (n == 0) {
    NumericVector nodes = vec_comb(data(_,0), data(_,1), true);
    n = max(nodes);
  }
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
#' @example
#' edgelist <- cbind(c(1,1,3,6),c(4,3,200,1))
#' edgelist
#' recode(edgelist)
recode <- function(...) UseMethod("recode")

#' @describeIn recode Method for data.frame
recode.data.frame <- function(data, ...) {
  cn <- colnames(data)
  data <- as.data.frame(recode.matrix(as.matrix(data), ...))
  colnames(data) <- cn
  data
}

#' @describeIn recode Method for matrix
recode.matrix <- function(data, ...) {

  # Checking the size of the matrix

  data <- as.factor(as.vector(data))
  n <- length(data)
  cbind(data[1:(n/2)], data[(n/2+1):n])
}

# Important difference with the previous version, this one accounts for duplicate
# dyads and also for self edges.
edgelist_to_adjmat <- function(
  edgelist, weights=NULL,
  times=NULL,
  undirected=FALSE, skip.recode=FALSE, no.self=FALSE, no.multiple=FALSE) {

  # Checking dim of edgelist
  if (ncol(edgelist) !=2) stop("Edgelist must have 2 columns")

  # Checking out the weights
  m <- nrow(edgelist)
  if (!length(weights)) weights <- rep(1, m)

  # Recoding nodes ids
  if (!skip.recode) dat <- recode(edgelist)
  else {
    warning('Skipping -recode- may cause unexpected behavior.')
    dat <- edgelist
  }
  n <- max(dat)

  # Checking out duplicates and self
  if (no.self)     dat <- dat[dat[,1]!=dat[,2]]
  if (no.multiple) dat <- unique(dat)

  # Checking out times
  if (!length(times)) times <- rep(1, m)
  t <- max(times)

  # Computing the adjmat
  adjmat <- vector("list", t)

  for(i in 1:length(adjmat)) {
    index <- which(times == i)
    adjmat[[i]] <- edgelist_to_adjmat_cpp(
      dat[index,,drop=FALSE], weights[index], n, undirected)
  }

  n <- nrow(adjmat[[1]])
  return(array(unlist(adjmat), dim=c(n,n,t)))
}

# Base data
set.seed(123)
n <- 2000
edgelist <- matrix(sample(1:n, size = n*10, replace = TRUE), ncol=2)
times <- sample.int(10, nrow(edgelist), replace=TRUE)
w <- abs(rnorm(nrow(edgelist)))

# # Simple example
# edgelist_to_adjmat(edgelist)
# edgelist_to_adjmat(edgelist, undirected = TRUE)
#
# # Using weights
# edgelist_to_adjmat(edgelist, w)
# edgelist_to_adjmat(edgelist, w, undirected = TRUE)
#
# # Using times
# edgelist_to_adjmat(edgelist, times = times)
# edgelist_to_adjmat(edgelist, times = times, undirected = TRUE)
#
# # Using times and weights
# edgelist_to_adjmat(edgelist, times = times, weights = w)
# edgelist_to_adjmat(edgelist, times = times, undirected = TRUE, weights = w)

# Benchmark with the previous version
library(microbenchmark)
library(diffusiontest)

dat <- as.data.frame(cbind(edgelist, w))
colnames(dat) <- c('ego','alter','tie')
microbenchmark(
  adjmatbuild(dat,n,1:n),
  edgelist_to_adjmat(edgelist, w), times=10)

old <- adjmatbuild(dat[,-3],n,1:n)
new <- (edgelist_to_adjmat(unique(edgelist), undirected = FALSE))[,,1]
arrayInd(which(old!=new), dim(old), dimnames(old))

*/



