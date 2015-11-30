# Rcpp::sourceCpp("/home/george/Documents/usc/software/netdiffuseR/playground/adjmat.cpp")
# library(microbenchmark)
# library(netdiffuseR)

# Important difference with the previous version, this one accounts for duplicate
# dyads and also for self edges.

#' Conversion between adjacency matrix and edgelist
#'
#' Generates adjacency matrix from an edgelist and viceversa.
#'
#' @param edgelist Two column matrix/data.frame in the form of ego -source- and
#' alter -target- (see details).
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @param weights Numeric vector. Strength of ties (optional).
#' @param times Integer vector. Periodicity of the ties (optional).
#' @param t Integer scalar. If \code{times} but want to repeat the network \code{t} times.
#' @param simplify Logical scalar. When TRUE and \code{times=NULL} it will return an adjacency
#' matrix, otherwise an array of adjacency matrices.
#' @param undirected Logical scalar. TRUE when the graph is undirected.
#' @param self Logical scalar. TRUE when self edges are excluded.
#' @param multiple Logical scalar. TRUE when multiple edges should not be included
#' (see details).
#' @param use.incomplete Logical scalar. When FALSE, rows with \code{NA/NULL} values will be droped
#' and will not be considered in the graph, which may reduce the size of the
#' adjacency matrix (see
#' details).
#' @param recode.ids Logical scalar. When TRUE ids are recoded using \code{\link{as.factor}}
#' (see details).
#' @param ... Further arguments to be passed to the matrix method.
#' @details The edgelist must be coded from 1:n (otherwise it may cause an error).
#' By default, the function will \code{\link{recode}} the edgelist before starting.
#'
#' When multiple edges are included, each vertex between \eqn{\{i,j\}}{{i,j}} will be counted
#' as many times it appears in the edgelist. So if a vertex \eqn{\{i,j\}}{{i,j}} appears 2
#' times, the adjacency matrix element \code{(i,j)} will be 2.
#'
#' Including incomplete cases, \code{use.incomplete=TRUE}, can lead to an adjacency matrix
#' with isolated vertices. Otherwise, when \code{use.incomplete=FALSE}, if all the
#' edges in which a vertex participates have incomplete information in any of the
#' variables (a NA, NULL or NaN value, see \code{\link{complete.cases}}), it
#' will be dropped from the graph, thus, reducing the size of the adjacency
#' matrix by not including isolated vertices.
#'
#' The function performs several checks before starting to create the adjacency
#' matrix. These are:
#' \itemize{
#'  \item{Dimensions of the inputs, such as number of columns and length of vectors}
#'  \item{Having complete cases. If anly edge has a non-numeric value such as NAs or
#'  NULL in any variable (including \code{times} and \code{weights}), it will be
#'  removed. A full list of such edges can be retrieved from the attribute
#'  \code{incomplete}}
#'  \item{Nodes and times ids coding}
#' }
#'
#' \code{recode.ids=FALSE} is useful when the vertices ids have already been
#' coded. For example, after having use \code{adjmat_to_edgelist}, ids are
#' correctly encoded, so when going back (using \code{edgelist_to_adjmat})
#' \code{recode.ids} should be FALSE.
#'
#' @return In the case of \code{edgelist_to_adjmat} either an adjacency matrix
#' (if times is NULL) or an array of these (if times is not null). For
#' \code{adjmat_to_edgelist} the output is an edgelist.
#' @export
#' @examples
#' # Base data
#' set.seed(123)
#' n <- 5
#' edgelist <- rgraph_er(n, as.edgelist=TRUE)
#' times <- sample.int(3, nrow(edgelist), replace=TRUE)
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
#' @keywords manip
#' @family data management functions
edgelist_to_adjmat <- function(
  edgelist, weights=NULL,
  times=NULL, t=NULL, simplify=TRUE,
  undirected=FALSE, self=FALSE, multiple=FALSE,
  use.incomplete=TRUE, recode.ids=TRUE, ...) {

  switch (class(edgelist),
    data.frame = edgelist_to_adjmat.data.frame(
      edgelist, weights, times, t, simplify, undirected, self, multiple,
      use.incomplete, recode.ids, ...
    ),
    matrix = edgelist_to_adjmat.matrix(
      edgelist, weights, times, t, simplify, undirected, self, multiple,
      use.incomplete, recode.ids, ...)
  )
}

# @rdname edgelist_to_adjmat
# @export
edgelist_to_adjmat.data.frame <- function(
  edgelist, weights=NULL,
  times=NULL, t=NULL, simplify=TRUE,
  undirected=FALSE, self=FALSE, multiple=FALSE,
  use.incomplete=TRUE, recode.ids=TRUE, ...) {
  edgelist_to_adjmat.matrix(as.matrix(edgelist), ...)
}

# @rdname edgelist_to_adjmat
# @export
edgelist_to_adjmat.matrix <- function(
  edgelist, weights=NULL,
  times=NULL, t=NULL, simplify=TRUE,
  undirected=FALSE, self=FALSE, multiple=FALSE,
  use.incomplete=TRUE, recode.ids=TRUE, ...) {

  # Step 0: Checking dimensions
  if (ncol(edgelist) !=2) stop("Edgelist must have 2 columns")
  if (length(times) && nrow(edgelist) != length(times)) stop("-times- should have the same length as number of rows in -edgelist-")
  if (length(weights) && nrow(edgelist) != length(weights)) stop("-weights- should have the same length as number of rows in -edgelist-")

  ##############################################################################
  # Step 1: Incomplete cases.
  # Finding incomplete cases. This is always done since we need to provide a
  # complete list of variables to the C++ function, otherwise it will throw
  # an error.
  complete <- complete.cases(cbind(edgelist, times, weights))
  incomplete <- which(!complete)

  # If the user chooses to drop incomplete, then what changes is the selection
  # of ids in the graph.
  # Times and weights MUST be removed since wont be used in the C++ function
  if (length(incomplete))
    warning("Some vertices had NA/NULL values:\n\t",
            paste0(head(incomplete,20), collapse = ", "),
            ifelse(length(incomplete)>20,", ...", ""),
            "\nThe complete list will be stored as an attribute of the resulting",
            " adjacency matrix, namely, -incomplete-.")

  # If no.incomplete is activated, then the vectors should be fixed
  if (!use.incomplete) edgelist <- edgelist[complete,,drop=FALSE]

  ##############################################################################
  # Step 2: Recoding nodes ids
  # Recoding nodes ids
  if (recode.ids) dat <- recode(edgelist)
  else dat <- edgelist

  n <- max(dat, na.rm = TRUE)

  ##############################################################################
  # Step 3: Preparing -times- and -weights- considering complete cases.
  # Times + recoding
  m <- nrow(edgelist)
  if (length(times)) times <- times[complete]
  else times <- rep(1, m)

  oldtimes <- range(times)
  oldtimes <- oldtimes[1]:oldtimes[2]
  times    <- times - min(times, na.rm = TRUE) + 1L

  if (!length(t)) t <- max(times, na.rm = TRUE)

  # Weights
  if (length(weights)) weights <- weights[complete]
  else weights <- rep(1, m)

  ##############################################################################
  # Computing the graph
  graph <- vector("list", t)

  if (recode.ids) labs <- attr(dat, "recode")[["label"]]
  else labs <- 1:n

  for(i in 1:t) {
    index <- which(times <= i)
    graph[[i]] <- edgelist_to_adjmat_cpp(
      dat[index,,drop=FALSE], weights[index], n, undirected, self, multiple)

    # Naming
    dimnames(graph[[i]]) <- list(labs, labs)
  }

  # Times naming
  if (t>1 && (length(oldtimes) == 1))
    oldtimes <- 1:t

  names(graph) <- oldtimes

  if (t==1 & simplify) graph <- graph[[1]]

  attr(graph, "incomplete") <- incomplete

  return(graph)
}

#' @rdname edgelist_to_adjmat
#' @export
adjmat_to_edgelist <- function(graph, undirected=TRUE) {

  switch (class(graph),
          list      = adjmat_to_edgelist.list(graph, undirected),
          array     = adjmat_to_edgelist.array(graph, undirected),
          dgCMatrix = adjmat_to_edgelist.dgCMatrix(graph, undirected),
          matrix    = adjmat_to_edgelist.matrix(graph, undirected)
  )
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.matrix <- function(graph, undirected=TRUE) {
  adjmat_to_edgelist_cpp(methods::as(graph, "dgCMatrix"), undirected)
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.dgCMatrix <- function(graph, undirected=TRUE) {
  adjmat_to_edgelist_cpp(graph, undirected)
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.array <- function(graph, undirected=TRUE) {
  edgelist <- matrix(ncol=2,nrow=0)
  times <- vector('integer',0L)
  for (i in 1:dim(graph)[3]) {
    x <- adjmat_to_edgelist.matrix(graph[,,i], undirected)
    edgelist <- rbind(edgelist, x)
    times <- c(times, rep(i,nrow(x)))
  }

  return(list(edgelist, times))
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.list <- function(graph, undirected=TRUE) {
  edgelist <- matrix(ncol=2,nrow=0)
  times <- vector('integer',0L)
  for (i in 1:length(graph)) {
    x <- adjmat_to_edgelist.dgCMatrix(graph[[i]], undirected)
    edgelist <- rbind(edgelist, x)
    times <- c(times, rep(i,nrow(x)))
  }

  return(list(edgelist, times))
}

# # Benchmark with the previous version
# library(microbenchmark)
# library(netdiffuseR)
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

#' Time of adoption matrix
#'
#' Creates two matrices indicating the time of adoption (TOA) of the innovation, and times during
#' which an individual had already adopted the innovation.
#'
#' @param times Integer vector of size \eqn{n} containing time of adoption of the innovation.
#' @param recode Logical scalar. When TRUE recodes time (see details).
#' @param labels Character vector of size \eqn{n}. Labels (ids) of the vertices.
#' @param ... Ignored.
#' @details
#'
#' By construction this function requires time units to be between 1 and T, where
#' T is the length (number) of time periods. By default a recoding of time is
#' performed as \eqn{time'=time - min(time) + 1}{time'=time - min(time) + 1}.
#'
#' @export
#' @return A list of two \eqn{n \times T}{n x T}
#'  \item{\code{cumadopt}}{has 1's for all years in which a node indicates having the innovation.}
#'  \item{\code{adopt}}{has 1's only for the year of adoption and 0 for the rest.}
#' @keywords manip
toa_mat <- function(times, recode=TRUE, labels=NULL, ...) {
  switch(class(times),
    numeric = toa_mat.numeric(times, recode, labels, ...),
    integer = toa_mat.integer(times, recode, labels, ...)
  )
  # UseMethod("toa_mat")
}

# @rdname toa_mat
# @export
toa_mat.numeric <- function(times, recode=TRUE, labels=NULL, ...) {
  if (inherits(times, 'numeric')) warning('-x- numeric. will be coersed to integer.')
  times <- as.integer(times)
  toa_mat.integer(times, recode)
}

# @rdname toa_mat
# @export
toa_mat.integer <- function(times, recode=TRUE, labels=NULL, ...) {
  # Rescaling
  oldtimes <- range(times)
  oldtimes <- oldtimes[1]:oldtimes[2]
  if (recode) times <- times - min(times) + 1L
  output <- toa_mat_cpp(times)

  # Naming
  if (length(labels)) {
    rownames(output$cumadopt) <- labels
    rownames(output$adopt) <- labels
  }
  colnames(output$cumadopt) <- oldtimes
  colnames(output$adopt) <- oldtimes
  output
}

# set.seed(123)
# x <- sample(2000:2005, 10, TRUE)
# y <- as.numeric(as.factor(x))
#
# new <- toa_mat(x)
# old <- adoptMat(y)
#
# sum(new[[1]] - old[[1]])
# sum(new[[2]] - old[[2]])
#
# microbenchmark(adoptMat(y), toa_mat_cpp(x), times=1000)
# Unit: microseconds
#             expr    min     lq      mean median      uq      max neval cld
# adoptMat(y)      43.876 51.010 61.133262 53.002 55.9400 4070.201 10000   b
# toa_mat_cpp(x)  4.620  6.226  7.921307  7.374  8.2605  114.874 10000  a

#' Difference in Time of Adoption (TOA) between individuals
#'
#' Creates \eqn{n \times n}{n * n} matrix indicating the difference in times of adoption between
#' each pair of nodes
#' @inheritParams toa_mat
#' @details Each cell ij of the resulting matrix is calculated as \eqn{toa_j - toa_i}{%
#' toa(j) - toa(i)}, so that whenever its positive it means that the j-th individual (alter)
#' adopted the innovation sonner.
#' @return An \eqn{n \times n}{n * n} symmetric matrix indicating the difference in times of
#' adoption between each pair of nodes.
#' @export
#' @examples
#' # Generating a random vector of time
#' set.seed(123)
#' times <- sample(2000:2005, 10, TRUE)
#'
#' # Computing the TOA differences
#' toa_diff(times)
#' @keywords manip
toa_diff <- function(times, recode=TRUE, labels=NULL, ...) {
  switch (class(times),
    integer = toa_diff.integer(times, recode, labels, ...),
    numeric = toa_diff.numeric(times, recode, labels, ...)
  )
  # UseMethod("toa_diff")
}

# @rdname toa_diff
# @export
toa_diff.integer <- function(times, recode=TRUE, labels=NULL,...) {
  # Rescaling
  if (recode) times <- times - min(times) + 1L
  toa_diff_cpp(times)
}

# @rdname toa_diff
# @export
toa_diff.numeric <- function(times, recode=TRUE, labels=NULL,...) {
  times <- as.integer(times)
  toa_diff.integer(times, recode)
}

# set.seed(123)
# x <- sample(2000:2005, 10, TRUE)
# toa_diff(x)
#
# microbenchmark(toaMat(x), toa_diff_cpp(x), times=1000)
# Unit: microseconds
#           expr     min      lq       mean   median       uq      max neval cld
# toaMat(x)      227.279 247.679 291.940566 272.4290 283.6845 3667.118  1000   b
# toa_diff_cpp(x)   3.539   4.623   6.887954   6.3755   7.4645   54.817  1000  a
# > 291.940566/6.887954
# [1] 42.38422


#' Find and remove isolated vertices
#'
#' Find and remove unconnected vertices from the graph.
#'
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @param undirected Logical. TRUE when the graph is undirected.
#' @export
#' @return
#' When \code{graph} is an adjacency matrix:
#'  \item{isolated}{an matrix of size \eqn{n\times 1}{n*1} with 1's where a node is isolated}
#'  \item{drop_isolated}{a modified graph excluding isolated vertices.}
#'
#' Otherwise, when \code{graph} is a list
#'  \item{isolated}{an matrix of size \eqn{n\times T}{n*T} with 1's where a node is isolated}
#'  \item{drop_isolated}{a modified graph excluding isolated vertices.}
#' @examples
#' # Generating random graph
#' set.seed(123)
#' adjmat <- rgraph_er()
#'
#' # Making nodes 1 and 4 isolated
#' adjmat[c(1,4),] <- 0
#' adjmat[,c(1,4)] <- 0
#' adjmat
#'
#' # Finding isolated nodes
#' iso <- isolated(adjmat)
#' iso
#'
#' # Removing isolated nodes
#' drop_isolated(adjmat)
#'
#'
#' # Now with a dynamic graph
#' graph <- rgraph_er(n=10, t=3)
#'
#' # Making 1 and 5 isolated
#' graph <- lapply(graph, "[<-", i=c(1,5), j=1:10, value=0)
#' graph <- lapply(graph, "[<-", i=1:10, j=c(1,5), value=0)
#' graph
#'
#' isolated(graph)
#' drop_isolated(graph)
#' @keywords manip
#' @family data management functions
isolated <- function(graph, undirected=getOption("diffnet.undirected")) {
  switch (class(graph),
    matrix = isolated.matrix(graph, undirected),
    dgCMatrix = isolated.dgCMatrix(graph, undirected),
    array = isolated.array(graph, undirected),
    list = isolated.list(graph, undirected)
  )
  UseMethod("isolated")
}

# @export
# @rdname isolated
isolated.matrix <- function(graph, undirected=getOption("diffnet.undirected")) {
  out <- isolated_cpp(methods::as(graph, "dgCMatrix"), undirected)
  dimnames(out) <- list(rownames(graph), "isolated")
  out
}

# @export
# @rdname isolated
isolated.dgCMatrix <- function(graph, undirected=getOption("diffnet.undirected")) {
  out <- isolated_cpp(graph, undirected)
  dimnames(out) <- list(rownames(graph), "isolated")
  out
}

# @export
# @rdname isolated
isolated.array <- function(graph, undirected=getOption("diffnet.undirected")) {
  nper <- dim(graph)[3]
  n    <- dim(graph)[2]

  # Creating output list and anciliary vector (to see if is isolated or not!)
  iso  <- Matrix::Matrix(0, ncol=nper, nrow=n, sparse=TRUE)
  for(i in 1:nper)
    iso[,i] <- isolated_cpp(methods::as(graph[,,i], "dgCMatrix"), undirected)

  isolated <- ifelse(apply(iso, 1, sum)==nper, 1, 0)

  # Naming
  dimnames(iso) <- list(rownames(graph), dimnames(graph)[[3]])
  names(isolated) <- rownames(graph)

  list(
    isolated_t=iso,
    isolated=isolated
  )
}

# @export
# @rdname isolated
isolated.list <- function(graph, undirected=getOption("diffnet.undirected")) {
  nper<- length(graph)
  n   <- nrow(graph[[1]])

  # Creating output list and anciliary vector (to see if is isolated or not!)
  iso  <- Matrix::Matrix(0, ncol=nper, nrow=n, sparse=TRUE)
  for(i in 1:nper)
    iso[,i] <- isolated_cpp(graph[[i]], undirected)
  isolated <- structure(
    ifelse(apply(iso, 1, sum)==nper, 1, 0),
    dim=c(n,1), dimnames=list(rownames(graph[[1]]), "isolated")
  )

  # Naming
  dimnames(iso) <- list(rownames(graph[[1]]), names(graph))

  list(
    isolated_t=iso,
    isolated=isolated
  )
}

#' @export
#' @rdname isolated
drop_isolated <- function(graph, undirected=getOption("diffnet.undirected")) {
  switch (class(graph),
    matrix = drop_isolated.matrix(graph, undirected),
    list = drop_isolated.list(graph, undirected),
    dgCMatrix = drop_isolated.dgCMatrix(graph, undirected),
    array = drop_isolated.array(graph, undirected)
  )
  # UseMethod("drop_isolated")
}

# @rdname isolated
# @export
drop_isolated.matrix <- function(graph, undirected=getOption("diffnet.undirected")) {
  iso <- isolated(graph, undirected)
  out <- drop_isolated_cpp(methods::as(graph, "dgCMatrix"), iso, undirected)

  # Indexing the set of non-zero elements
  iso <- rownames(iso[which(iso==0),,drop=FALSE])
  dimnames(out) <- list(iso, iso)
  out
}

# @rdname isolated
# @export
drop_isolated.dgCMatrix <- function(graph, undirected=getOption("diffnet.undirected")) {
  iso <- isolated(graph, undirected)
  out <- drop_isolated_cpp(graph, iso, undirected)

  # Indexing the set of non-zero elements
  iso <- rownames(iso[which(iso==0),,drop=FALSE])
  dimnames(out) <- list(iso, iso)
  out
}

# @rdname isolated
# @export
drop_isolated.array <- function(graph, undirected=getOption("diffnet.undirected")) {
  # Getting isolated vecs
  iso <- isolated.array(graph, undirected)[[2]]
  ison <- rownames(iso[which(iso==0),,drop=FALSE])

  m   <- sum(iso)
  n   <- dim(graph)[1]
  t   <- dim(graph)[3]
  out <- vector("list", t)
  names(out) <- names(graph)

  # Removing
  for(i in 1:t) {
    out[[i]] <- drop_isolated_cpp(methods::as(graph[,,i], "dgCMatrix"), iso, undirected)
    dimnames(out[[i]]) <- list(ison, ison)
  }

  out
}

# @rdname isolated
# @export
drop_isolated.list <- function(graph, undirected=getOption("diffnet.undirected")) {
  # Getting isolated vecs
  iso <- isolated.list(graph, undirected)[[2]]
  ison <- rownames(iso[which(iso==0),,drop=FALSE])
  m   <- sum(iso)
  n   <- nrow(graph[[1]])
  t   <- length(graph)
  out <- vector("list", t)
  names(out) <- names(graph)

  # Removing
  for(i in 1:t) {
    out[[i]] <- drop_isolated_cpp(graph[[i]], iso, undirected)
    dimnames(out[[i]]) <- list(ison,ison)
  }

  out
}
