Rcpp::sourceCpp("/home/george/Documents/usc/software/diffusiontest/playground/adjmat.cpp")
library(microbenchmark)
library(diffusiontest)

#' Erdos-Renyi (bernoulli) random graph
#'
#' Follows the G(N,p) model
#'
#' @param n Number of vertices
#' @param p Probability of connection between ego and alter.
#' @param undirected whether the graph is undirected or not.
#' @param weighted Whether the graph is weighted or not.
#' @param self Wheter it includes self-edges.
#' @return A graph represented by an adjacency matrix
#' @note The resulting adjacency matrix is dense (hence, be careful with the size)
#' @example
#' \dontrun{
#' # Setting the seed
#' set.seed(123)
#'
#' # Generating an directed graph
#' rand_graph(undirected=FALSE)
#'
#' # Comparing P(tie)
#' x <- rand_graph(1000, p=.1)
#' sum(x)/length(x)
#' }
rand_graph <- function(n=10, p=0.3, undirected=TRUE, weighted=FALSE, self=FALSE) {
  rand_graph_cpp(n, p, undirected, weighted, self)
}

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
recode <- function(data, ...) UseMethod("recode")

#' @describeIn recode Method for data frame
#' @export
recode.data.frame <- function(data, ...) {
  cn <- colnames(data)
  data <- as.data.frame(recode.matrix(as.matrix(data), ...))
  colnames(data) <- cn
  data
}

#' @describeIn recode Method for matrix
#' @export
recode.matrix <- function(data, ...) {

  # Checking the size of the matrixRcppArmadilloForward.h

  data <- as.factor(as.vector(data))
  n <- length(data)
  cbind(data[1:(n/2)], data[(n/2+1):n])
}

# Important difference with the previous version, this one accounts for duplicate
# dyads and also for self edges.
#' @param edgelist
#' @param weights
#' @param times
#' @param simplify
#' @param undirected
#' @param skip.recode
#' @param no.self
#' @param no.multiple
#' @example
#' # Base data
#' set.seed(123)
#' n <- 10
#' edgelist <- matrix(sample(1:n, size = n*10, replace = TRUE), ncol=2)
#' times <- sample.int(10, nrow(edgelist), replace=TRUE)
#' w <- abs(rnorm(nrow(edgelist)))
#'
#' # Simple example
#' edgelist_to_adjmat(edgelist)
#' edgelist_to_adjmat(edgelist, undirected = TRUE)
#'
#' # Using weights
#' edgelist_to_adjmat(edgelist, w)
#' edgelist_to_adjmat(edgelist, w, undirected = TRUE)
#'
#' # Using times
#' edgelist_to_adjmat(edgelist, times = times)
#' edgelist_to_adjmat(edgelist, times = times, undirected = TRUE)
#'
#' # Using times and weights
#' edgelist_to_adjmat(edgelist, times = times, weights = w)
#' edgelist_to_adjmat(edgelist, times = times, undirected = TRUE, weights = w)

edgelist_to_adjmat <- function(edgelist, ...) UseMethod("edgelist_to_adjmat")

#' @describeIn edgelist_to_adjmat Method for data frame
#' @export
edgelist_to_adjmat.data.frame <- function(edgelist, ...)
  edgelist_to_adjmat(as.matrix(edgelist), ...)

#' @describeIn edgelist_to_adjmat Method for matrix object
#' @export
edgelist_to_adjmat.matrix <- function(
  edgelist, weights=NULL,
  times=NULL, simplify=TRUE,
  undirected=FALSE, skip.recode=FALSE, no.self=FALSE, no.multiple=FALSE) {

  # Checking dim of edgelist (and others)
  if (ncol(edgelist) !=2) stop("Edgelist must have 2 columns")
  if (nrow(edgelist) != length(times)) stop("-times- should have the same length as number of rows in -edgelist-")
  if (nrow(edgelist) != length(weights)) stop("-weights- should have the same length as number of rows in -edgelist-")

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
  if (simplify & !length(times)) return(adjmat[[1]])
  else return(array(unlist(adjmat), dim=c(n,n,t)))
}


# # Benchmark with the previous version
# library(microbenchmark)
# library(diffusiontest)
#
# dat <- as.data.frame(cbind(edgelist, w))
# colnames(dat) <- c('ego','alter','tie')
# microbenchmark(
#   adjmatbuild(dat,n,1:n),
#   edgelist_to_adjmat(edgelist, w), times=100)
#
# old <- adjmatbuild(dat[,-3],n,1:n)
# new <- (edgelist_to_adjmat(unique(edgelist), undirected = FALSE))[,,1]
# arrayInd(which(old!=new), dim(old), dimnames(old))
#
# ## Dynamic
# microbenchmark(
#   adjByTime(cbind(year=times,dat),n,max(times)),
#   edgelist_to_adjmat(edgelist, w, times), times=100)


#'
#' @param time Integer vector containing time of adoption
#' @returns A n x T matrix with times of adoption.
adopt_mat <- function(time, ...) UseMethod("adopt_mat")

#' @describeIn adopt_mat Method for integer vector
#' @export
adopt_mat.integer <- function(time, ...) adopt_mat_cpp(x)

#' @describeIn adopt_mat Method for numeric vector
#' @export
adopt_mat.numeric <- function(time, ...) {
  if (inherits(time, 'numeric')) warning('-x- numeric. will be coersed to integer.')
  time <- as.integer(time)
  time <- time + min(time) + 1
  adopt_mat_cpp(time)
}


# set.seed(123)
# x <- sample(2000:2005, 10, TRUE)
# y <- as.numeric(as.factor(x))
#
# new <- adopt_mat(x)
# old <- adoptMat(y)
#
# sum(new[[1]] - old[[1]])
# sum(new[[2]] - old[[2]])
#
# microbenchmark(adoptMat(y), adopt_mat_cpp(x), times=1000)
# Unit: microseconds
#             expr    min     lq      mean median      uq      max neval cld
# adoptMat(y)      43.876 51.010 61.133262 53.002 55.9400 4070.201 10000   b
# adopt_mat_cpp(x)  4.620  6.226  7.921307  7.374  8.2605  114.874 10000  a

#' @param time Integer vector with time of adoption
toa_mat <- function(x,...) UseMethod("toa_mat")

#' @describeIn toa_mat Method for numeric vector
toa_mat.numeric <- function(x,...) {
  toa_mat_cpp(x)
}

set.seed(123)
x <- sample(2000:2005, 10, TRUE)
toa_mat(x)

microbenchmark(toaMat(x), toa_mat_cpp(x), times=1000)
# Unit: microseconds
#           expr     min      lq       mean   median       uq      max neval cld
# toaMat(x)      227.279 247.679 291.940566 272.4290 283.6845 3667.118  1000   b
# toa_mat_cpp(x)   3.539   4.623   6.887954   6.3755   7.4645   54.817  1000  a

