#' Indexing diffnet objects (on development)
#'
#' Access and assign (replace) elements from the adjacency matrices or the vertex
#' attributes data frames.
#'
#' @param x A diffnet class object.
#' @param name String vector. Names of the vertices attributes.
#' @param i Index of the i-th row of the adjacency matrix (see details).
#' @param j Index of the j-th column of the adjacency matrix (see details)
#' @param k Index of the k-th slice of the adjacency matrix (see details).
#' @param value Value to assign (see details)
#' @param drop Logical scalar. When \code{TRUE} returns an adjacency matrix, otherwise
#' a filtered diffnet object.
#' @param as.df Logical scalar. When \code{TRUE} returns a data frame, otherwise
#' a list of length \eqn{T}.
#' @param ... Further argumnets to be passed to the method (on development)
#' @details
#' The \code{[[.diffnet} methods provides access to the diffnet attributes
#' data frames, static and dynamic. By providing the \code{name} of the corresponding
#' attribute, depending on whether it is static or dynamic the function will return
#' either a data frame--static attributes--or a list of these--dynamic attributes.
#' For the assigning method, \verb{[[<-.diffnet}, the function will infer what
#' kind of attribute is by analyzing the dimmensions of \code{value}, in particular
#' we have the following possible cases:
#'
#' \tabular{llr}{
#' \strong{Class}    \tab \strong{Dimension}   \tab \strong{Inferred} \cr
#' \code{matrix}     \tab \eqn{n\times T}{n*T} \tab Dynamic           \cr
#' \code{matrix}     \tab \eqn{n\times 1}{n*1} \tab Static            \cr
#' \code{matrix}     \tab \eqn{(n\times T)\times 1}{(n*T)*1} \tab Dynamic            \cr
#' \code{data.frame} \tab \eqn{n\times T}{n*T} \tab Dynamic           \cr
#' \code{data.frame} \tab \eqn{n\times 1}{n*1} \tab Static            \cr
#' \code{data.frame} \tab \eqn{(n\times T)\times 1}{(n*T)*1} \tab Dynamic            \cr
#' \code{vector}     \tab \eqn{n}              \tab Static            \cr
#' \code{vector}     \tab \eqn{n\times T}{n*T} \tab Dynamic           \cr
#' \code{list}*      \tab \eqn{T} data.frames/matrices/vectors\tab Dynamic    \cr
#' }
#' *: With \eqn{n\times 1}{n * 1} \code{data.frame}/\code{matrix} or \eqn{n} length vector.
#'
#' Other cases will return with error.
#'
#' In the case of the slices index \code{k}, either an
#' integer vector with the positions, a character vector with the labels of the
#' time periods or a logical vector of length \code{T} can be used to specify
#' which slices to retrieve. Likewise, indexing vertices works in the same way
#' with the only difference that, instead of time period labels and a logical vector
#' of length \code{T}, vertices ids labels and a logical vector of length \code{n}
#' should be provided.
#'
#' When subsetting slices, the function modifies the \code{toa} vector as well as the
#' \code{adopt} and \code{cumadopt} matrices collapsing network tinmming. For example,
#' if a network goes from time 1 to 20 and we set \code{k=3:10}, all individuals
#' who adopted prior to time 3 will be set as adopters at time 3, and all individuals
#' who adopted after time 10 will be set as adopters at time 10, changing the
#' adoption and cumulative adoption matrices. Importantly, \code{k} have no
#' gaps, and it should be within the graph time period range.
#'
#' @return In the case of the asignning methods, a diffnet object. Otherwise,
#' for \code{[[.diffnet} a vector extracted from one of the attributes data frames,
#' and for \code{[.diffnet} a list of length \code{length(k)} with the corresponding
#' \code{[i,j]} elements from the adjacency matrix.
#' @name diffnet_index
#' @family diffnet methods
#' @author George G. Vega Yon
#' @examples
#'
#' # Creating a random diffusion network ---------------------------------------
#' set.seed(111)
#' graph <- rdiffnet(100,5)
#'
#' # Accessing to a static attribute
#' graph[["real_threshold"]]
#'
#' # Accessing to subsets of the adjacency matrix
#' graph[1,,1:3, drop=TRUE]
#' graph[,,1:3, drop=TRUE]
#'
#' # ... Now, as diffnet objects (the default)
#' graph[1,,1:3, drop=FALSE]
#' graph[,,1:3, drop=FALSE]
#'
#' # Changing values in the adjacency matrix
#' graph[1, , , drop=TRUE]
#' graph[1,,] <- -5
#' graph[1, , , drop=TRUE]
#'
#' # Adding attributes (dynamic) -----------------------------------------------
#' # Preparing the data
#' set.seed(1122)
#' x <- rdiffnet(30, 5, seed.p.adopt=.15)
#'
#' # Calculating exposure, and storing it diffe
#' expoM <- exposure(x)
#' expoL <- lapply(seq_len(x$meta$nper), function(x) expoM[,x,drop=FALSE])
#' expoD <- do.call(rbind, expoL)
#'
#' # Adding data (all these are equivalent)
#' x[["expoM"]] <- expoM
#' x[["expoL"]] <- expoL
#' x[["expoD"]] <- expoD
#'
#' # Lets compare
#' identical(x[["expoM"]], x[["expoL"]]) # TRUE
#' identical(x[["expoM"]], x[["expoD"]]) # TRUE
NULL

# @export
# @rdname diffnet_index
diffnet.subset.slices <- function(graph, k) {

  if (!inherits(graph, "diffnet")) stop("-graph- must be a 'diffnet' object")

  # Subset must be of length 2 (at least)
  if (length(k) < 2)
    stop("-k- is a vector of length ",length(k),
         ". It must be at least of length 2.")

  # Analyzing class
  uses_labels <- ifelse(inherits(k, "character"), TRUE, FALSE)
  if (uses_labels) {
    k <- as.integer(k)
  }

  # If logical vector
  if (inherits(k, "logical")) {
    k <- which(k)
  }

  # Subset must be continuous...
  test <- (k[-1] - k[-length(k)]) > 1
  if (any(test))
    stop("-k- must represent a range without gaps.")

  # Ordering
  k  <- sort(k)
  nslices <- length(k)
  pers    <- graph$meta$pers

  # Checking k
  test <- if (!uses_labels) !(k %in% seq_len(graph$meta$nper))
  else !(k %in% pers)

  if (any(test))
    stop("The specified -k- (",
         paste0(k[test], collapse = ", "),
         ") are invalid.")

  # Recomputing in terms of indexes
  if (uses_labels) k <- which(pers %in% as.character(k))

  # Removing not included k
  graph$graph            <- graph$graph[k]
  graph$vertex.dyn.attrs <- graph$vertex.dyn.attrs[k]

  # Changing adoption matrices

  beforeslice <- which(graph$toa < pers[k][1])
  afterslice  <- which(graph$toa > pers[k][nslices])

  graph$adopt[beforeslice,k[1]] <- 1
  graph$adopt[afterslice ,k[nslices]] <- 1
  graph$adopt    <- graph$adopt[,k]

  graph$cumadopt[beforeslice,k] <- 1
  graph$cumadopt[afterslice,k[nslices]] <- 1
  graph$cumadopt <- graph$cumadopt[,k]

  # Changing toa mat (truncating it)
  graph$toa[beforeslice] <- pers[k][1]
  graph$toa[afterslice]  <- pers[k][nslices]

  # Changing meta
  graph$meta$nper <- nslices
  graph$meta$pers <- pers[k]

  graph
}

#' Infer whether \code{value} is dynamic or static.
#'
#' Intended for internal use only, this function is used in \code{\link{diffnet_index}}
#' methods.
#'
#' @param value Either a matrix, data frame or a list. Attribute values.
#' @param meta A list. A diffnet object's meta data.
#'
#' @return The value object either as a data frame (if static) or as a list
#' of data frames (if dynamic). If \code{value} does not follows the permitted
#' types of \code{\link{diffnet_index}}, then returns with error.
#'
#' @export
diffnet_check_attr_class <- function(value, meta) {
  # Processing meta
  n <- meta$n
  t <- meta$nper

  if (inherits(value, "matrix") | inherits(value, "data.frame")) {
    # Static case
    vdim <- dim(value)
    if (all(vdim == c(n, t)))    { # Dynamic Matrix

      # Coercing into a list
      value <- lapply(seq_len(t), function(x) value[, x, drop=TRUE])
      return(value)

    } else if (all(vdim == c(n, 1L))) { # Static Matrix

      # Coercing into a vector
      return(value[,1,drop=TRUE])

    } else if (all(vdim == c(n*t, 1L))) { # Dynamic Matrix

      # Coercing into a list
      value <- lapply(seq_len(t), function(x)
        value[((x-1)*n + 1):(x*n),,drop=TRUE]
        )
      return(value)

    } else {
      stop("-value- data.frame/matrix has incorrect size (", vdim[1], " x ", vdim[2],"). ",
           "Please refer to the manual to see accepted values.")
    }

  } else if (inherits(value, "list")) {

    # Checking classes of the elements
    isdf <- sapply(value, inherits, what="data.frame")
    # Which ones aren't either
    test <- !isdf & !sapply(value, inherits, what="matrix")

    # If no data.frame/matrix, then no vector?
    test <- which(ifelse(!test, test, !is.vector(value)))

    if (length(test)) {
      stop("Not all the elements in the list are data.frame/matrix or vector:\n\t",
           paste0(test, collapse=", "), ".")
    }

    # Checking the dimensions of the elements
    test <- which(sapply(value, function(x) {
      if (is.vector(x)) {
        length(x) != n
      }
      else {
        !all(dim(x) == c(n,1))
      }
      }))

    if (length(test)) {
      stop("Not all elements in the list have the right dimension (",n," elements):\n\t",
           paste0(test, collapse=", "), ".")
    }

    # Corecing non df to dfs, and then coercing into vectors
    value[!isdf] <- lapply(value[!isdf], as.data.frame)
    value <- lapply(value, "[[", 1)

    return(value)

  } else if (is.vector(value)) {
    # Checking the length
    vdim <- length(value)

    if (vdim == n) {# Checking static

      # Returning asis (ideal case?)
      return(value)

    } else if (vdim == (n*t)) { # Checking dynamic

      # Coercing into a list
      value <- lapply(seq_len(t), function(x)
        value[((x-1)*n + 1):(x*n)]
      )

      return(value)

    } else {
      stop("-value- vector has incorrect size (", vdim, "). ",
           "Please refer to the manual to see accepted values.")
    }

  } else stop("-value- must be either a list, a data frame or a matrix.")
}

#' @export
#' @rdname diffnet_index
`[[.diffnet` <- function(x, name, as.df=FALSE) {

  # Checking names
  if (name %in% colnames(x$vertex.static.attrs)) {
    x <- x$vertex.static.attrs[[name]]
    if (as.df) x <- as.data.frame(setNames(list(x), list(name)))
    return(x)
  } else if (name %in% colnames(x$vertex.dyn.attrs[[1]])) {

    x <- lapply(x$vertex.dyn.attrs, "[[", name)
    if (as.df) {
      x <- as.data.frame(do.call(c, x))
      colnames(x) <- name
      return(x)
    }
    else return(x)

  } else {
    stop("No dynamic or static attribute with such name.")
  }

}

#' @export
#' @rdname diffnet_index
`[[<-.diffnet` <- function(x, i, j, value) {
  # If j index is specified, then the addition is made to a subset
  if (missing(j)) j <- seq_len(x$meta$n)

  # Checking and preparing the data
  meta <- x$meta
  if (length(j) < meta$n) meta$n <- length(j)
  value <- diffnet_check_attr_class(value, meta)

  # Adding the attribute
  if (!inherits(value, "list")) {
    x$vertex.static.attrs[[i]][j] <- value
  } else {
    # Checking if is empty or not
    for (l in 1:meta$nper)
      x$vertex.dyn.attrs[[l]][[i]][j] <- value[[l]]
  }
  x
}

#' @export
#' @rdname diffnet_index
`[.diffnet` <- function(x, i, j, k, drop=FALSE) {

  # Checking drop
  if (drop) { # So it requires the adjmat
    # Checking ids
    if (missing(i)) i <- seq_len(x$meta$n)
    if (missing(j)) j <- seq_len(x$meta$n)
  } else { # So its subsetting the diffnet
    if      (missing(i)  & missing(j) ) i <- j <- seq_len(x$meta$n)
    else if (missing(i)  & !missing(j)) i <- j
    else if (!missing(i) & missing(j) ) j <- i
    else if (!missing(i) & !missing(j))
      if (!identical(i, j))
        stop("Whe subsetting a diffnet and -i- and -j- are provided these should,",
             "be identical.")
  }

  # Slices
  if (missing(k)) k <- seq_len(x$meta$nper)

  # Subsetting
  if (drop) return(lapply(x$graph[k], "[", i=i, j=j, drop=FALSE))
  else {
    x <- diffnet.subset.slices(x, k)
    x$graph <- lapply(x$graph, "[", i=i, j=j, drop=FALSE)

    # Subsetting
    # 1.0: graph and attributes
    for (l in 1:x$meta$nper)
       x$vertex.dyn.attrs[[l]] <- x$vertex.dyn.attrs[[l]][i,,drop=FALSE]


    # 2.0: Matrices
    x$adopt               <- x$adopt[i,,drop=FALSE]
    x$cumadopt            <- x$cumadopt[i,,drop=FALSE]
    x$vertex.static.attrs <- x$vertex.static.attrs[i,,drop=FALSE]
    x$toa                 <- x$toa[i]

    # 3.0: Attrubytes
    x$meta$ids <- rownames(x$adopt)#x$meta$ids[i]
    x$meta$n   <- length(x$meta$ids)

    return(x)
  }
}

#' @export
#' @rdname diffnet_index
`[<-.diffnet` <- function(x, i, j, k, value) {
  # Checking ids
  if (missing(i)) i <- seq_len(x$meta$n)
  if (missing(j)) j <- seq_len(x$meta$n)
  if (missing(k)) k <- seq_len(x$meta$nper)

  for (l in k)
    x$graph[[l]][i,j] <- value

  x
}


