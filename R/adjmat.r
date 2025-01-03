#' Conversion between adjacency matrix and edgelist
#'
#' Generates adjacency matrix from an edgelist and vice versa.
#'
#' @param edgelist Two column matrix/data.frame in the form of ego -source- and
#' alter -target- (see details).
#' @templateVar undirected TRUE
#' @templateVar self TRUE
#' @templateVar multiple TRUE
#' @template graph_template
#' @param w Numeric vector. Strength of ties (optional).
#' @param t0 Integer vector. Starting time of the ties (optional).
#' @param t1 Integer vector. Finishing time of the ties (optional).
#' @param t Integer scalar. Repeat the network \code{t} times (if no \code{t0,t1} are provided).
#' @param simplify Logical scalar. When TRUE and \code{times=NULL} it will return an adjacency
#' matrix, otherwise an array of adjacency matrices.
#' (see details).
#' @param keep.isolates Logical scalar. When FALSE, rows with \code{NA/NULL} values
#' (isolated vertices unless have autolink) will be droped (see details).
#' @param recode.ids Logical scalar. When TRUE ids are recoded using \code{\link{as.factor}}
#' (see details).
#' @details
#'
#' When converting from edglist to adjmat the function will \code{\link{recode}} the
#' edgelist before starting. The user can keep track after the recording by checking
#' the resulting adjacency matrices' \code{\link{row.names}}. In the case that the
#' user decides skipping the recoding (because wants to keep vertices index numbers,
#' implying that the resulting graph will have isolated vertices), he can override
#' this by setting \code{recode.ids=FALSE} (see example).
#'
#' When multiple edges are included, \code{multiple=TRUE},each vertex between \eqn{\{i,j\}}{{i,j}} will be counted
#' as many times it appears in the edgelist. So if a vertex \eqn{\{i,j\}}{{i,j}} appears 2
#' times, the adjacency matrix element \code{(i,j)} will be 2.
#'
#' Edges with incomplete information (missing data on \code{w} or \code{times}) are
#' not included on the graph. Incomplete cases are tagged using \code{\link{complete.cases}}
#' and can be retrieved by the user by accessing the attribute \code{incomplete}.
#'
#' Were the case that either ego or alter are missing (i.e. \code{NA} values), the
#' function will either way include the non-missing vertex. See below for an example
#' of this.
#'
#' The function performs several checks before starting to create the adjacency
#' matrix. These are:
#' \itemize{
#'  \item{Dimensions of the inputs, such as number of columns and length of vectors}
#'  \item{Having complete cases. If anly edge has a non-numeric value such as NAs or
#'  NULL in either \code{times} or \code{w}, it will be
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
#'
#'
#' @return In the case of \code{edgelist_to_adjmat} either an adjacency matrix
#' (if times is NULL) or an array of these (if times is not null). For
#' \code{adjmat_to_edgelist} the output is an edgelist with the following columns:
#' \item{ego}{Origin of the tie.}
#' \item{alter}{Target of the tie.}
#' \item{value}{Value in the adjacency matrix.}
#' \item{time}{Either a 1 (if the network is static) or the time stamp of the tie.}
#' @export
#' @examples
#' # Base data
#' set.seed(123)
#' n <- 5
#' edgelist <- rgraph_er(n, as.edgelist=TRUE, p=.2)[,c("ego","alter")]
#' times <- sample.int(3, nrow(edgelist), replace=TRUE)
#' w <- abs(rnorm(nrow(edgelist)))
#'
#' # Simple example
#' edgelist_to_adjmat(edgelist)
#' edgelist_to_adjmat(edgelist, undirected = TRUE)
#'
#' # Using w
#' edgelist_to_adjmat(edgelist, w)
#' edgelist_to_adjmat(edgelist, w, undirected = TRUE)
#'
#' # Using times
#' edgelist_to_adjmat(edgelist, t0 = times)
#' edgelist_to_adjmat(edgelist, t0 = times, undirected = TRUE)
#'
#' # Using times and w
#' edgelist_to_adjmat(edgelist, t0 = times, w = w)
#' edgelist_to_adjmat(edgelist, t0 = times, undirected = TRUE, w = w)
#'
#' # Not recoding ----------------------------------------------------
#' # Notice that vertices 3, 4 and 5 are not present in this graph.
#' graph <- matrix(c(
#'  1,2,6,
#'  6,6,7
#' ), ncol=2)
#'
#' # Generates an adjmat of size 4 x 4
#' edgelist_to_adjmat(graph)
#'
#' # Generates an adjmat of size 7 x 7
#' edgelist_to_adjmat(graph, recode.ids=FALSE)
#'
#' # Dynamic with spells -------------------------------------------------------
#' edgelist <- rbind(
#'    c(1,2,NA,1990),
#'    c(2,3,NA,1991),
#'    c(3,4,1991,1992),
#'    c(4,1,1992,1993),
#'    c(1,2,1993,1993)
#' )
#'
#' graph <- edgelist_to_adjmat(edgelist[,1:2], t0=edgelist[,3], t1=edgelist[,4])
#'
#' # Creating a diffnet object with it so we can apply the plot_diffnet function
#' diffnet <- as_diffnet(graph, toa=1:4)
#' plot_diffnet(diffnet, label=rownames(diffnet))
#'
#' # Missing alter in the edgelist ---------------------------------------------
#' data(fakeEdgelist)
#'
#' # Notice that edge 202 is isolated
#' fakeEdgelist
#'
#' # The function still includes vertex 202
#' edgelist_to_adjmat(fakeEdgelist[,1:2])
#'
#' edgelist
#'
#' @keywords manip
#' @family data management functions
#' @include graph_data.r
#' @author George G. Vega Yon & Thomas W. Valente
edgelist_to_adjmat <- function(
  edgelist, w=NULL,
  t0=NULL, t1=NULL, t=NULL, simplify=TRUE,
  undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self"), multiple=getOption("diffnet.multiple"),
  keep.isolates=TRUE, recode.ids=TRUE) {

  cls <- class(edgelist)
  if ("data.frame" %in% cls)
    edgelist_to_adjmat.data.frame(
      edgelist, w, t0, t1, t, simplify, undirected, self, multiple,
      keep.isolates, recode.ids
    )
  else if ("matrix" %in% cls)
    edgelist_to_adjmat.matrix(
      edgelist, w, t0, t1, t, simplify, undirected, self, multiple,
      keep.isolates, recode.ids)
  else stop("-edgelist- should be either a data.frame, or a matrix.")
}

# @rdname edgelist_to_adjmat
# @export
edgelist_to_adjmat.data.frame <- function(
  edgelist, w,
  t0,t1, t, simplify,
  undirected, self, multiple,
  keep.isolates, recode.ids) {

  edgelist_to_adjmat.matrix(as.matrix(edgelist), w, t0,t1, t, simplify,
                            undirected, self, multiple, keep.isolates,
                            recode.ids)
}

# @rdname edgelist_to_adjmat
# @export
edgelist_to_adjmat.matrix <- function(
  edgelist, w,
  t0, t1,
  t, simplify,
  undirected, self, multiple,
  keep.isolates, recode.ids) {

  # Step 0: Checking dimensions
  if (ncol(edgelist) !=2) stop("Edgelist must have 2 columns")
  if (length(t0) && nrow(edgelist) != length(t0))
    stop("-t0- should have the same length as number of rows in -edgelist-.",
         " Currently they are ",length(t0), " and ", nrow(edgelist),
         " respectively.")
  if (length(t1) && nrow(edgelist) != length(t1))
    stop("-t1- should have the same length as number of rows in -edgelist-.",
         " Currently they are ",length(t1), " and ", nrow(edgelist),
         " respectively.")
  if (length(w) && nrow(edgelist) != length(w))
    stop("-w- should have the same length as number of rows in -edgelist-.",
         " Currently they are ",length(w), " and ", nrow(edgelist),
         " respectively.")

  ##############################################################################
  # Step 1: Incomplete cases.
  # Finding incomplete cases. This is always done since we need to provide a
  # complete list of variables to the C++ function, otherwise it will throw
  # an error.
  if (length(w)) complete <- complete.cases( w)
  else complete <- rep(TRUE, nrow(edgelist))
  edgelist <- edgelist[complete,]

  incomplete <- which(!complete)

  # Getting the list of isolated vertices
  not.isolated <- which(!is.na(edgelist[,1]) & !is.na(edgelist[,2]))

  # If the user chooses to drop incomplete, then what changes is the selection
  # of ids in the graph.
  # Times and w MUST be removed since wont be used in the C++ function
  if (length(incomplete))
    warning("Some edges a had NA/NULL value on either -times- or -w-:\n\t",
            paste0(head(incomplete,20), collapse = ", "),
            ifelse(length(incomplete)>20,", ...", ""),
            "\nThese won't be included in the adjacency matrix. The complete list will be stored as an attribute of the resulting",
            " adjacency matrix, namely, -incomplete-.")


  # If no.incomplete is activated, then the vectors should be fixed
  if (!keep.isolates) edgelist <- edgelist[not.isolated,,drop=FALSE]

  ##############################################################################
  # Step 2: Recoding nodes ids
  # Recoding nodes ids
  if (recode.ids) {
    # Recoding
    dat <- recode(edgelist)
    n <- max(dat, na.rm = TRUE)

    # Retrieving codes + labels. If use.incomplete == FALSE, then the edgelist
    # has already been filtered
    recodedids <- attr(dat, "recode")
    if (keep.isolates) dat <- dat[not.isolated,]
    attr(dat, "recode") <- recodedids
    rm(recodedids)
  }
  else {
    dat <- edgelist
    n <- max(dat, na.rm = TRUE)
  }

  ##############################################################################
  # Step 3: Preparing -times- and -w- considering complete cases.
  # Times + recoding
  m <- nrow(dat)
  if (length(t0)) t0 <- t0[complete][not.isolated]
  else t0 <- rep(NA, m)

  if (length(t1)) t1 <- t1[complete][not.isolated]
  else t1 <- rep(NA, m)

  suppressWarnings(oldtimes <- range(c(t0,t1), na.rm=TRUE))
  if (all(!is.finite(oldtimes))) oldtimes <- rep(1,2)
  oldtimes <- oldtimes[1]:oldtimes[2]

  t0    <- t0 - oldtimes[1] + 1L
  t1    <- t1 - oldtimes[1] + 1L

  # Replacing NAs
  t0[is.na(t0)] <- 1
  if (!length(t)) t <- max(c(t0,t1), na.rm=TRUE)
  t1[is.na(t1)] <- t

  # Weights
  if (length(w)) w <- w[complete][not.isolated]
  else w <- rep(1, m)

  ##############################################################################
  # Computing the graph
  graph <- vector("list", t)

  if (recode.ids) labs <- attr(dat, "recode")[["label"]]
  else labs <- 1:n

  for(i in 1:t) {
    index <- which((t0 <= i) & (i <= t1))
    graph[[i]] <- edgelist_to_adjmat_cpp(
      dat[index,,drop=FALSE], w[index], n, undirected, self, multiple)

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
adjmat_to_edgelist <- function(
  graph,
  undirected=getOption("diffnet.undirected", FALSE),
  keep.isolates = getOption("diffnet.keep.isolates", TRUE)
  ) {

  cls <- class(graph)
  out <- if ("list" %in% cls) {
      adjmat_to_edgelist.list(graph, undirected, keep.isolates)
    } else if ("matrix" %in% cls) {
      cbind(adjmat_to_edgelist.matrix(graph, undirected, keep.isolates),1)
    } else if ("array" %in% cls) {
      adjmat_to_edgelist.array(graph, undirected, keep.isolates)
    } else if ("dgCMatrix" %in% cls) {
      cbind(adjmat_to_edgelist.dgCMatrix(graph, undirected, keep.isolates), 1)
    } else
      stopifnot_graph(graph)

  colnames(out) <- c("ego", "alter", "value", "time")
  out
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.matrix <- function(graph, undirected, keep.isolates) {
  # as.integer(factor(diffnet$meta$ids[fakeDynEdgelist[,1]], diffnet$meta$ids))
  out <- adjmat_to_edgelist_cpp(methods::as(graph, "dgCMatrix"), undirected)

  # If keep isolates
  if (keep.isolates) {
    N <- 1:nvertices(graph)
    test <- which(!(N %in% unlist(out[,1:2])))

    # If there are isolates
    if (length(test)) out <- rbind(out, cbind(test, NA, NA))
  }

  return(out)
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.dgCMatrix <- function(graph, undirected, keep.isolates) {
  out <- adjmat_to_edgelist_cpp(graph, undirected)

  # If keep isolates
  if (keep.isolates) {
    N <- 1:nvertices(graph)
    test <- which(!(N %in% unlist(out[,1:2])))

    # If there are isolates
    if (length(test)) out <- rbind(out, cbind(test, NA, NA))
  }

  return(out)
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.array <- function(graph, undirected, keep.isolates) {
  edgelist <- matrix(ncol=3,nrow=0)
  times <- vector('integer',0L)

  # Getting the names
  tnames <- as.integer(dimnames(graph)[[3]])
  if (!length(tnames)) tnames <- 1:dim(graph)[3]

  for (i in 1:dim(graph)[3]) {
    x <- adjmat_to_edgelist.matrix(graph[,,i], undirected, keep.isolates)
    edgelist <- rbind(edgelist, x)
    times <- c(times, rep(tnames[i],nrow(x)))
  }

  return(cbind(edgelist, times=times))
}

# @rdname edgelist_to_adjmat
# @export
adjmat_to_edgelist.list <- function(graph, undirected, keep.isolates) {
  edgelist <- matrix(ncol=3,nrow=0)
  times <- vector('integer',0L)

  # Getting the names
  tnames <- as.integer(names(graph))
  if (!length(tnames)) tnames <- 1:dim(graph)[3]

  for (i in 1:length(graph)) {
    x <- adjmat_to_edgelist.dgCMatrix(graph[[i]], undirected, keep.isolates)
    edgelist <- rbind(edgelist, x)
    times <- c(times, rep(tnames[i],nrow(x)))
  }

  # # Adjusting the length
  # ids <- apply(edgelist, 1, paste0, collapse="")
  # times <- as.integer(unname(tapply(times, ids, min)))
  #
  # edgelist <- unique(edgelist)
  # edgelist <- edgelist[order(edgelist[,1],edgelist[,2]),]

  return(cbind(edgelist, times=times))
}

#' Time of adoption matrix
#'
#' For a single behavior, creates two matrices recording times of adoption of the innovation. One matrix
#' records the time period of adoption for each node with zeros elsewhere. The
#' second records the cumulative time of adoption such that there are ones for
#' the time of adoption and every time period thereafter. For \eqn{Q} behaviors,
#' creates a list of length \eqn{Q}, where each element contains those two
#' matrices for each behavior.
#'
#' @param obj Either an integer vector of length \eqn{n} containing time of adoption
#' of the innovation, a matrix of size \eqn{n \times Q} (for multiple \eqn{Q} behaviors), or
#' a \code{\link{diffnet}} object (both for single or multiple behaviors).
#' @param labels Character vector of length \eqn{n}. Labels (ids) of the vertices.
#' @param t0 Integer scalar. Sets the lower bound of the time window (e.g. 1955).
#' @param t1 Integer scalar. Sets the upper bound of the time window (e.g. 2000).
#' @details
#'
#' In order to be able to work with time ranges other than \eqn{1,\dots, T}{1,..., T}
#' the function receives as input the boundary labels of the time windows through
#' the variables \code{t0} and \code{t}. While by default the function assumes that
#' the the boundaries are given by the range of the \code{times} vector, the user
#' can set a personalized time range exceeding the one given by the \code{times}
#' vector. For instance, times of adoption may range between 2001 and 2005 but the
#' actual data, the network, is observed between 2000 and 2005 (so there is not
#' left censoring in the data), hence, the user could write:
#'
#' \preformatted{
#' adopmats <- toa_mat(times, t0=2000, t1=2005)
#' }
#'
#' That way the resulting \code{cumadopt} and \code{adopt} matrices would have
#' 2005 - 2000 + 1 = 6 columns instead of 2005 - 2001 + 1 = 5 columns, with the
#' first column of the two matrices containing only zeros (as the first adoption
#' happend after the year 2000).
#'
#' For multiple behaviors, the input can be a matrix or a \code{diffnet} object.
#' In this case, the output will be a list, with each element replicating the output
#' for a single diffusion: a matrix recording the time period of adoption for
#' each node, and a second matrix with ones from the moment the node adopts the behavior.
#'
#' @examples
#' # Random set of times of adoptions
#' times <- sample(c(NA, 2001:2005), 10, TRUE)
#'
#' toa_mat(times)
#'
#' # Now, suppose that we observe the graph from 2000 to 2006
#' toa_mat(times, t0=2000, t1=2006)
#'
#' # For multiple behaviors, the input can be a matrix..
#' times_1 <- c(2001L, 2004L, 2003L, 2008L)
#' times_2 <- c(2001L, 2005L, 2006L, 2008L)
#' times <- matrix(c(times_1, times_2), nrow = 4, ncol = 2)
#'
#' toa <- toa_mat(times)
#' toa[[1]]$adopt         # time period of adoption for the first behavior
#'
#' #.. or a diffnet object
#' graph <- lapply(2001:2008, function(x) rgraph_er(4))
#' diffnet <- new_diffnet(graph, times)
#'
#' toa <- toa_mat(diffnet)
#' toa[[1]]$cumadopt      # cumulative adoption matrix for the first behavior

#'
#' @export
#' @return For a single behavior, a list of two \eqn{n \times T}{n x T}:
#'  \item{\code{cumadopt}}{ has 1's for all years in which a node indicates having the innovation.}
#'  \item{\code{adopt}}{ has 1's only for the year of adoption and 0 for the rest.}
#'  For \eqn{Q} behaviors, a list of length \eqn{Q}, each element containing
#'  \code{cumadopt} ans \code{adopt} matrices.
#' @keywords manip
#' @include graph_data.r
#' @author George G. Vega Yon, Thomas W. Valente, and Aníbal Olivera M.
toa_mat <- function(obj, labels=NULL, t0=NULL, t1=NULL) {

  if (inherits(obj, "matrix")) {
    num_of_behaviors <- dim(obj)[2]
  } else if (inherits(obj, "diffnet")){
    if (inherits(obj$toa, "matrix")) {
      num_of_behaviors <- dim(obj$toa)[2]}
    else {num_of_behaviors <- 1}
  } else {num_of_behaviors <- 1}

  if (!inherits(obj, "diffnet")) {
    if (!length(t0)) t0 <- min(obj, na.rm = TRUE)
    if (!length(t1)) t1 <- max(obj, na.rm = TRUE)
  }

  ans <- list()
  if (num_of_behaviors == 1) {
    cls <- class(obj)
    ans[[1]] <- if ("numeric" %in% cls) {
            toa_mat.numeric(obj, labels, t0, t1)
            } else if ("integer" %in% cls) {
            toa_mat.integer(obj, labels, t0, t1)
            } else if  ("diffnet" %in% cls) {
            with(obj, list(adopt=adopt,cumadopt=cumadopt))
            } else {
              stopifnot_graph(obj)
            }
  } else {
    for (q in 1:num_of_behaviors) {
      ans[[q]] <- if ("matrix" %in% class(obj)) {
              if ("integer" %in% class(obj[,q])){
                toa_mat.integer(obj[,q], labels, t0, t1)
              } else if ("numeric" %in% class(obj[,q])) { # Why included?
                toa_mat.numeric(obj[,q], labels, t0, t1)
              }
            } else if  ("diffnet" %in% class(obj)) { # Why included?
              with(obj, list(adopt=adopt[[q]],cumadopt=cumadopt[[q]]))
            } else {
              stopifnot_graph(obj[,q])
            }
    }
  }

  for (q in 1:num_of_behaviors) {
    if (inherits(obj, "diffnet")) {
      dimnames(ans[[q]]$adopt) <- with(obj$meta, list(ids,pers))
      dimnames(ans[[q]]$cumadopt) <- with(obj$meta, list(ids,pers))
    }
  }

  if (num_of_behaviors==1) {
    return(ans[[1]])
  } else {
    return(ans)
  }
}

toa_mat.default <- function(per, t0, t1) {
  ans <- matrix(0L, ncol=t1-t0+1L, nrow=length(per))
  ans[cbind(1L:nrow(ans), per - t0 + 1L)] <- 1L

  list(
    adopt = ans,
    cumadopt = t(apply(ans, 1, cumsum))
  )
}

# @rdname toa_mat
# @export
toa_mat.numeric <- function(times, labels=NULL,
                            t0 = min(times, na.rm = TRUE),
                            t1 = max(times, na.rm=TRUE)) {
  if (inherits(times, 'numeric')) warning('-x- numeric. will be coersed to integer.')

  # Coercing into integer
  times <- as.integer(times)
  t0 <- as.integer(t0)
  t1 <- as.integer(t1)

  toa_mat.integer(times, labels, t0, t1)
}

# @rdname toa_mat
# @export
toa_mat.integer <- function(times, labels=NULL,
                            t0 = min(times, na.rm = TRUE),
                            t1 = max(times, na.rm=TRUE)) {
  # Rescaling
  output <- toa_mat.default(times, t0, t1)

  # Naming
  cn <- t0:t1
  if (length(labels)) rn <- labels
  else rn <- 1:length(times)

  dimnames(output[[1]]) <- list(rn, cn)
  dimnames(output[[2]]) <- list(rn, cn)

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
#' Creates an \eqn{n \times n}{n * n} matrix, or for \eqn{Q}{Q} behaviors, a list
#' of length \eqn{Q}{Q} containing \eqn{n \times n}{n * n} matrices, that indicates
#' the difference in adoption times between each pair of nodes.
#' @inheritParams toa_mat
#' @details Each cell \eqn{ij}{ij} of the resulting matrix is calculated as \eqn{toa_j - toa_i}{%
#' toa(j) - toa(i)}, so that whenever its positive it means that the j-th individual (alter)
#' adopted the innovation sooner.
#' @return An \eqn{n \times n}{n * n} anti-symmetric matrix (or a list of them,
#' for \eqn{Q}{Q} behaviors) indicating the difference in times of
#' adoption between each pair of nodes.
#' @export
#' @examples
#' # For a single behavior -----------------------------------------------------
#'
#' # Generating a random vector of time
#' set.seed(123)
#' times <- sample(2000:2005, 10, TRUE)
#'
#' # Computing the TOA differences
#' toa_diff(times)
#'
#' # For Q=2 behaviors ---------------------------------------------------------
#'
#' # Generating a matrix time
#'
#' times_1 <- c(2001L, 2004L, 2003L, 2008L)
#' times_2 <- c(2001L, 2005L, 2006L, 2008L)
#' times <- matrix(c(times_1, times_2), nrow = 4, ncol = 2)
#'
#' # Computing the TOA differences
#' toa_diff(times)
#'
#' # Or, from a diffnet object
#'
#' graph <- lapply(2001:2008, function(x) rgraph_er(4))
#' diffnet <- new_diffnet(graph, times)
#'
#' # Computing the TOA differences
#' toa_diff(diffnet)
#'

#'
#' @keywords manip
#' @include graph_data.r
#' @author George G. Vega Yon, Thomas W. Valente, and Aníbal Olivera M.
toa_diff <- function(obj, t0=NULL, labels=NULL) {

  # Calculating t0 (if it was not provided)
  if (!inherits(obj, "diffnet") && !length(t0)){
    t0 <- as.integer(min(obj, na.rm = TRUE))
  } else {
    t0 <- obj$meta$pers[1]}

  # determining num_of_behavior and prepare for multi-diffusion
  num_of_behavior <- 1
  multiple <- FALSE

  if (inherits(obj, "matrix")) { # multiple
    num_of_behavior <- ncol(obj)
    obj <- lapply(asplit(obj, MARGIN = 2), as.integer)
    multiple <- TRUE
  } else if (inherits(obj, "diffnet")) {
    if (inherits(obj$toa, "matrix")) { # multiple
      num_of_behavior <- ncol(obj$toa)
      obj <- split_behaviors(obj)
      multiple <- TRUE
    }
  }

  if (multiple) {
    out_list <- lapply(seq_len(num_of_behavior), function(q) toa_diff.unique(obj[[q]], t0))
    return(out_list)
  } else {
    return(toa_diff.unique(obj, t0))
  }
}

#
#
#   if (multiple) {
#     for (q in 1:ncol(obj$toa)) {
#
#
#       # Calculating t0 (if it was not provided)
#       if (!inherits(obj, "diffnet") && !length(t0)) {
#         t0 <- as.integer(min(obj[,q], na.rm = TRUE))
#       } else {
#         t0 <- obj$meta$pers[1]}
#
#       # Computing the difference
#       if (inherits(obj, "integer")) {
#         out <- toa_diff_cpp(obj - t0 + 1L)
#       } else if (inherits(obj, "numeric")) {
#         warning("coercing -obj- to integer.")
#         out <- toa_diff_cpp(as.integer(obj) - t0 + 1L)
#       } else if (inherits(obj, "diffnet")) {
#         out <- toa_diff_cpp(obj$toa - t0 + 1L)
#       } else stop("No method defined for class -",class(obj),"-")
#
#       out
#
#     }
#
#
#   } else {
#     # Calculating t0 (if it was not provided)
#     if (!inherits(obj, "diffnet") && !length(t0))
#       t0 <- as.integer(min(obj, na.rm = TRUE))
#     else
#       t0 <- obj$meta$pers[1]
#
#     # Computing the difference
#     if (inherits(obj, "integer")) {
#       out <- toa_diff_cpp(obj - t0 + 1L)
#     } else if (inherits(obj, "numeric")) {
#       warning("coercing -obj- to integer.")
#       out <- toa_diff_cpp(as.integer(obj) - t0 + 1L)
#     } else if (inherits(obj, "diffnet")) {
#       out <- toa_diff_cpp(obj$toa - t0 + 1L)
#     } else stop("No method defined for class -",class(obj),"-")
#
#     return(out)
#   }
# }

toa_diff.unique <- function(obj, t0) {
  # Computing the difference
  if (inherits(obj, "integer")) {
    out <- toa_diff_cpp(obj - t0 + 1L)
  } else if (inherits(obj, "numeric")) {
    warning("coercing -obj- to integer.")
    out <- toa_diff_cpp(as.integer(obj) - t0 + 1L)
  } else if (inherits(obj, "diffnet")) {
    out <- toa_diff_cpp(obj$toa - t0 + 1L)
  } else stop("No method defined for class -",class(obj),"-")

  return(out)
}

# @rdname toa_diff
# @export
toa_diff.integer <- function(times, t0, labels) {
  # Rescaling
  times <- times - t0 + 1L
  toa_diff_cpp(times)
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
#' @templateVar undirected TRUE
#' @template graph_template
#' @templateVar undirected 1
#' @templateVar self 1
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
#' @author George G. Vega Yon
isolated <- function(
  graph,
  undirected = getOption("diffnet.undirected", FALSE),
  self = getOption("diffnet.self", FALSE)
) {
  cls <- class(graph)
  ans <- if ("matrix" %in% cls) {
      isolated.default(methods::as(graph, "dgCMatrix"), undirected, self)
    } else if ("dgCMatrix" %in% cls) {
      isolated.default(graph, undirected, self)
    } else if ("array" %in% cls) {
      lapply(apply(graph, 3, methods::as, Class="dgCMatrix"), isolated.default,
                       undirected=undirected, self=self)
    } else if ("list" %in% cls) {
      lapply(graph, isolated.default, undirected=undirected, self=self)
    } else if ("diffnet" %in% cls) {
      lapply(graph$graph, isolated.default, undirected=undirected, self=self)
    } else
      stopifnot_graph(graph)

  if (any(class(graph) %in% c("list", "diffnet", "array")))
    apply(do.call(cbind, ans), 1, all)
  else ans
  # UseMethod("isolated")
}

# @export
# @rdname isolated
isolated.default <- function(
  graph,
  undirected = getOption("diffnet.undirected", FALSE),
  self = getOption("diffnet.self", FALSE)
) {

  graph@x <- rep(1.0, length(graph@x))
  d <- Matrix::rowSums(graph)
  if (undirected)
    d <- d + Matrix::colSums(graph)

  if (!self)
    d <- d - Matrix::diag(graph)

  unname(d == 0)
}

isolated.list <- function(
  graph,
  undirected = getOption("diffnet.undirected", FALSE),
  self = getOption("diffnet.self", FALSE)
  ) {

  ids <- lapply(graph, isolated.default, undirected=undirected, self=self)

  apply(do.call(cbind, ids), 1, all)

}


#' @export
#' @rdname isolated
drop_isolated <- function(
  graph,
  undirected = getOption("diffnet.undirected", FALSE),
  self = getOption("diffnet.self", FALSE)
) {

  # Find isolates
  ids <- which(!isolated(graph))

  if (inherits(graph, "list"))
    lapply(graph, "[", i=ids, j=ids, drop=FALSE)
  else if (inherits(graph, "diffnet"))
    graph[ids,]
  else if (inherits(graph, "dgCMatrix"))
    graph[ids,,drop=FALSE][,ids,drop=FALSE]
  else if (inherits(graph, "array"))
    graph[ids,,,drop=FALSE][,ids,,drop=FALSE]

}


simmelian_mat <- function(graph, ...) {
#  tmethod <- if(isS4(graph)) getMethod("t", class(graph)) else t
  tmp <- graph & t(graph) #tmethod(graph)
  methods::as(tmp & (tmp %*% tmp), "dgCMatrix")
}



#' Approximate Geodesic Distances
#'
#' Computes approximate geodesic distance matrix using graph powers and keeping
#' the amount of memory used low.
#'
#' @template graph_template
#' @param n Integer scalar. Degree of approximation. Bigger values increase
#' precision (see details).
#' @param warn Logical scalar. When \code{TRUE}, it warns if the algorithm
#' performs less steps than required.
#'
#' @details
#'
#' While both \pkg{igraph} and \pkg{sna} offer very good and computationally
#' efficient routines for computing geodesic distances, both functions return
#' dense matrices, i.e. not sparse, which can be troublesome. Furthermore,
#' from the perspective of social network analysis, path lengths of more than 6 steps,
#' for example, may not be meaningful, or at least, relevant for the researcher.
#' In such cases, \code{approx_geodesic} serves as a solution to this problem,
#' computing geodesics up to the number of steps, \code{n}, desired, hence,
#' if \code{n = 6}, once the algorithm finds all paths of 6 or less steps it
#' will stop, returning a sparse matrix with zeros for those pairs of
#' vertices for which it was not able to find a path with less than \code{n}
#' steps.
#'
#' Depending on the graph size and density, \code{approx_geodesic}'s performance
#' can be compared to that of \code{\link[sna:geodist]{sna::geodist}}. Although,
#' as \code{n} increases, \code{geodist} becomes a better alternative.
#'
#' The algorithm was implemented using power graphs. At each itereation i the
#' power graph of order \code{i} is computed, and its values are compared
#' to the current values of the geodesic matrix (which is initialized in zero).
#'
#' \enumerate{
#' \item Initialize the output \code{ans(n, n)}
#' \item For \code{i=1} to \code{i < n} do
#' \enumerate{
#'   \item Iterate through the edges of \code{G^i}, if \code{ans} has a zero
#'   value in the corresponding row+column, replace it with \code{i}
#'   \item next
#' }
#' \item Replace all diagonal elements with a zero and return.
#' }
#'
#' This implementation can be more memory efficient that the aforementioned ones,
#' but at the same time it can be significant slower.
#'
#' \code{approx_geodist} is just an allias for \code{approx_geodesic}.
#'
#' @return A sparse matrix of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} of size
#' \code{nnodes(graph)^2} with geodesic distances up to \code{n}.
#'
#' @examples
#' # A very simple example -----------------------------------------------------
#' g <- ring_lattice(10, 3)
#' approx_geodesic(g, 6)
#' sna::geodist(as.matrix(g))[[2]]
#' igraph::distances(
#'   igraph::graph_from_adjacency_matrix(g, mode = "directed"),
#'   mode = "out"
#' )
#'
#' @aliases Geodesic Shortest-Path
#' @export
approx_geodesic <- function(graph, n = 6L, warn=FALSE) {
  cls <- class(graph)

  if ("dgCMatrix" %in% cls) {
    approx_geodesicCpp(graph, n, warn)
  } else if ("matrix" %in% cls) {
    approx_geodesicCpp(methods::as(graph, "dgCMatrix"), n, warn)
  } else if ("list" %in% cls) {
    lapply(graph, approx_geodesicCpp, n = n, warn = warn)
  } else if ("diffnet" %in% cls) {
    lapply(graph$graph, approx_geodesicCpp, n = n, warn = warn)
  } else if ("array" %in% cls) {
    apply(graph, 3, approx_geodesicCpp, n = n, warn = warn)
  } else stopifnot_graph(graph)
}

#' @rdname approx_geodesic
#' @export
approx_geodist <- approx_geodesic
