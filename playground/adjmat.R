Rcpp::sourceCpp("/home/george/Documents/usc/software/diffusiontest/playground/adjmat.cpp")
library(microbenchmark)
library(diffusiontest)

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

  # Checking the size of the matrixRcppArmadilloForward.h

  data <- as.factor(as.vector(data))
  n <- length(data)
  cbind(data[1:(n/2)], data[(n/2+1):n])
}

# Important difference with the previous version, this one accounts for duplicate
# dyads and also for self edges.
edgelist_to_adjmat <- function(
  edgelist, weights=NULL,
  times=NULL, simplify=TRUE,
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
  if (simplify & !length(times)) return(adjmat[[1]])
  else return(array(unlist(adjmat), dim=c(n,n,t)))
}

# Base data
set.seed(123)
n <- 10
edgelist <- matrix(sample(1:n, size = n*10, replace = TRUE), ncol=2)
times <- sample.int(10, nrow(edgelist), replace=TRUE)
w <- abs(rnorm(nrow(edgelist)))

# Simple example
edgelist_to_adjmat(edgelist)
edgelist_to_adjmat(edgelist, undirected = TRUE)

# Using weights
edgelist_to_adjmat(edgelist, w)
edgelist_to_adjmat(edgelist, w, undirected = TRUE)

# Using times
edgelist_to_adjmat(edgelist, times = times)
edgelist_to_adjmat(edgelist, times = times, undirected = TRUE)

# Using times and weights
edgelist_to_adjmat(edgelist, times = times, weights = w)
edgelist_to_adjmat(edgelist, times = times, undirected = TRUE, weights = w)

# Benchmark with the previous version
library(microbenchmark)
library(diffusiontest)

dat <- as.data.frame(cbind(edgelist, w))
colnames(dat) <- c('ego','alter','tie')
microbenchmark(
  adjmatbuild(dat,n,1:n),
  edgelist_to_adjmat(edgelist, w), times=100)
#
# old <- adjmatbuild(dat[,-3],n,1:n)
# new <- (edgelist_to_adjmat(unique(edgelist), undirected = FALSE))[,,1]
# arrayInd(which(old!=new), dim(old), dimnames(old))
#
# ## Dynamic
# microbenchmark(
#   adjByTime(cbind(year=times,dat),n,max(times)),
#   edgelist_to_adjmat(edgelist, w, times), times=100)


set.seed(123)
x <- sample(2000:2005, 10, TRUE)
y <- as.numeric(as.factor(x))

adopt_mat <- function(x,...) UseMethod("adopt_mat")
adopt_mat.integer <- function(x, ...) adopt_mat_cpp(x)
adopt_mat.numeric <- function(x, ...) {
  if (inherits(x, 'numeric')) warning('-x- numeric. will be coersed to integer.')
  x <- as.integer(x)
  x <- x + min(x) + 1
  adopt_mat_cpp(x)
}

new <- adopt_mat(x)
old <- adoptMat(y)

sum(new[[1]] - old[[1]])
sum(new[[2]] - old[[2]])

microbenchmark(adoptMat(y), adopt_mat_cpp(x), times=1000)
# Unit: microseconds
#             expr    min     lq      mean median      uq      max neval cld
# adoptMat(y)      43.876 51.010 61.133262 53.002 55.9400 4070.201 10000   b
# adopt_mat_cpp(x)  4.620  6.226  7.921307  7.374  8.2605  114.874 10000  a

toa_mat <- function(x,...) UseMethod("toa_mat")
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

