# Rcpp::sourceCpp("/home/george/Documents/usc/software/diffusiontest/playground/adjmat.cpp")
# library(microbenchmark)
# library(diffusiontest)

# Important difference with the previous version, this one accounts for duplicate
# dyads and also for self edges.

#' Conversion between adjacency matrix and edgelist
#'
#' Generates adjacency adjacency matrix from an edgelist and viceversa.
#'
#' @param edgelist Two column matrix/data.frame in the form of ego -source- and
#' alter -target- (see details).
#' @param adjmat Square matrix. An adjacency matrix.
#' @param weights Numeric vector. Strength of ties (optional).
#' @param times Integer vector. Periodicity of the ties (optional).
#' @param simplify Logical. When TRUE and no times vector it will return an adjacency
#' matrix, otherwise an array of adjacency matrices.
#' @param undirected Logical. TRUE when the graph is undirected.
#' @param skip.recode Logical. FALSE when recode of nodes's ids is performed (see details).
#' @param no.self Logical. TRUE when self edges are excluded.
#' @param no.multiple Logical. TRUE when multiple edges should not be included
#' (see details).
#' @param ... Further arguments for the method.
#' @details The edgelist must be coded from 1:n (otherwise it may cause an error).
#' By default, the function will \code{\link{recode}} the edgelist before starting.
#'
#' When multiple edges are included, each vertex between \{i,j\} will be accounted
#' as many times it appears in the edgelist. So if a vertex \{i,j\} appears 2
#' times, the adjacency matrix element (i,j) will have a 2.
#' @return In the case of \code{edgelist_to_adjmat} either an adjacency matrix
#' (if times is NULL) or an array of these (if times is not null). For
#' \code{adjmat_to_edgelist} the output is an edgelist.
#' @export
#' @examples
#' # Base data
#' set.seed(123)
#' n <- 10
#' edgelist <- rand_graph(n, as.edgelist=TRUE)
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

#' @rdname edgelist_to_adjmat
#' @export
edgelist_to_adjmat.data.frame <- function(edgelist, ...) {
  edgelist <- as.matrix(edgelist)
  NextMethod("edgelist_to_adjmat")
}


#' @rdname edgelist_to_adjmat
#' @export
edgelist_to_adjmat.matrix <- function(
  edgelist, weights=NULL,
  times=NULL, simplify=TRUE,
  undirected=FALSE, skip.recode=FALSE, no.self=FALSE, no.multiple=FALSE, ...) {

  # Checking out full observations (droping incomplete)
  index <- complete.cases(edgelist)
  edgelist <- edgelist[index,,drop=FALSE]
  if (length(times)) times <- times[index]
  if (length(weights)) weights <- weights[index]

  # For the output
  incomplete <- which(!index)
  if (length(incomplete))
    warning("Some vertices had NA/NULL values. These will not be considered",
            " in the graph: ",paste0(incomplete, collapse = ", "))

    # Checking dim of edgelist (and others)
  if (ncol(edgelist) !=2) stop("Edgelist must have 2 columns")
  if (length(times) && nrow(edgelist) != length(times)) stop("-times- should have the same length as number of rows in -edgelist-")
  if (length(weights) && nrow(edgelist) != length(weights)) stop("-weights- should have the same length as number of rows in -edgelist-")

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
  else {
    times <- times - min(times) + 1L
  }
  t <- max(times)

  # Computing the adjmat
  adjmat <- vector("list", t)

  for(i in 1:length(adjmat)) {
    index <- which(times == i)
    adjmat[[i]] <- edgelist_to_adjmat_cpp(
      dat[index,,drop=FALSE], weights[index], n, undirected)
  }

  n <- nrow(adjmat[[1]])
  if (t==1 & simplify) adjmat <- adjmat[[1]]
  else adjmat <- array(unlist(adjmat), dim=c(n,n,t))

  attr(adjmat, "incomplete") <- incomplete
  return(adjmat)
}

#' @rdname edgelist_to_adjmat
#' @export
adjmat_to_edgelist <- function(adjmat, undirected=TRUE) UseMethod("adjmat_to_edgelist")

#' @rdname edgelist_to_adjmat
#' @export
adjmat_to_edgelist.matrix <- function(adjmat, undirected=TRUE) {
  adjmat_to_edgelist_cpp(adjmat, undirected)
}

#' @rdname edgelist_to_adjmat
#' @export
adjmat_to_edgelist.array <- function(adjmat, undirected=TRUE) {
  edgelist <- matrix(ncol=2,nrow=0)
  times <- vector('integer',0L)
  for (i in 1:dim(adjmat)[3]) {
    x <- adjmat_to_edgelist.matrix(adjmat[,,i], undirected)
    edgelist <- rbind(edgelist, x)
    times <- c(times, rep(i,nrow(edgelist)))
  }

  return(list(edgelist, times))
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

#' Creates an adoption matrices
#' @param times Integer vector containing time of adoption of the innovation.
#' @param recode Logical. When TRUE recodes time (see details).
#' @param ... Ignored.
#' @details
#'
#' By construction this function requires time units to be between 1 and T, where
#' T is the length (number) of time periods. By default a recoding of time is
#' performed as \eqn{time'=time - min(time) + 1}{time'=time - min(time) + 1}.
#'
#' @export
#' @return A list of two n x T matrices with times of adoption.
#' \code{cumadopt} has 1's for all years in which a node indicates having the innovation.
#' \code{adopt} has 1's only for the year of adoption and 0 for the rest.
adopt_mat <- function(times, recode=TRUE, ...) UseMethod("adopt_mat")

#' @describeIn adopt_mat Integers
#' @export
adopt_mat.integer <- function(times, recode=TRUE, ...) {
  # Rescaling
  if (recode) times <- times - min(times) + 1L
  adopt_mat_cpp(times)
}

#' @describeIn adopt_mat Numeric
#' @export
adopt_mat.numeric <- function(times, recode=TRUE, ...) {
  if (inherits(times, 'numeric')) warning('-x- numeric. will be coersed to integer.')
  times <- as.integer(times)
  NextMethod("adopt_mat")
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

#' Difference in Time of Adoption (TOA) between nodes
#'
#' Creates n x n matrix indicating the difference in times of adoption between
#' each pair of nodes
#' @inheritParams adopt_mat
#' @return An n x n matrix indicating the difference in times of adoption between
#' each pair of nodes.
#' @export
#' @examples
#' # Generating a random vector of time
#' set.seed(123)
#' times <- sample(2000:2005, 10, TRUE)
#'
#' # Computing the TOA differences
#' toa_mat(times)
toa_mat <- function(times, recode=TRUE, ...) UseMethod("toa_mat")

#' @describeIn toa_mat Integer vectors
#' @export
toa_mat.integer <- function(times, recode=TRUE, ...) {
  # Rescaling
  if (recode) times <- times - min(times) + 1L
  toa_mat_cpp(times)
}

#' @describeIn toa_mat Numeric vectors
#' @export
toa_mat.numeric <- function(times, recode=TRUE, ...) {
  times <- as.integer(times)
  NextMethod("toa_mat")
}

# set.seed(123)
# x <- sample(2000:2005, 10, TRUE)
# toa_mat(x)
#
# microbenchmark(toaMat(x), toa_mat_cpp(x), times=1000)
# Unit: microseconds
#           expr     min      lq       mean   median       uq      max neval cld
# toaMat(x)      227.279 247.679 291.940566 272.4290 283.6845 3667.118  1000   b
# toa_mat_cpp(x)   3.539   4.623   6.887954   6.3755   7.4645   54.817  1000  a
# > 291.940566/6.887954
# [1] 42.38422


#' Manages isolated nodes
#' @param adjmat Square matrix. An graph as an adjacency matrix.
#' @param undirected Logical. TRUE when the graph is undirected.
#' @export
#' @return
#' In the case of \code{isolated}, an integer vector of size n with 1's where a
#' node is isolated. Otherwise a modified adjacency matrix excluding the
#' isolated nodes.
#' @examples
#' \dontrun{
#' # Generating random graph
#' set.seed(123)
#' adjmat <- rand_graph()
#'
#' # Making nodes 1 and 4 isolated
#' adjmat[c(1,4),] <- 0
#' adjmat[,c(1,4)] <- 0
#'
#' # Finding isolated nodes
#' isolated(adjmat)
#'
#' # Removing isolated nodes
#' drop_isolated(adjmat)
#' }
isolated <- function(adjmat, undirected=TRUE) {
  isolated_cpp(adjmat, undirected)
}

#' @rdname isolated
#' @export
drop_isolated <- function(adjmat, undirected=TRUE) {
  drop_isolated_cpp(adjmat, undirected)
}
