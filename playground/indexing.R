#' Indexing diffnet objects
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
#' @return In the case of the assignning methods, a diffnet object. Otherwise,
#' for \code{[[.diffnet} a vector extracted from one of the attributes data frames,
#' and for \code{[.diffnet} a list of length \code{length(k)} with the corresponding
#' \code{[i,j]} elements from the adjacency matrix.
#'
#' @details The \code{[[.diffnet} methods provides access to the diffnet attributes
#' data frames, static and dynamic. By providing the \code{name} of the corresponding
#' attribute, depending on whether it is static or dynamic the function will return
#' either a data frame--static attributes--or a list of these--dynamic attributes.
#' For the assigning method, \code{[[<-.diffnet}, the function will infer what
#' kind of attribute is by analyzing the dimensions of \code{value}, in particular
#' we have the following possible cases:
#'
#' \tabular{llr}{
#' \strong{Class}    \tab \strong{Dimension}   \tab \strong{Inferred} \cr
#' \code{matrix}     \tab \eqn{n\times T}{n*T} \tab Dynamic           \cr
#' \code{matrix}     \tab \eqn{n\times 1}{n*T1 \tab Static            \cr
#' \code{data.frame} \tab \eqn{n\times T}{n*T} \tab Dynamic           \cr
#' \code{data.frame} \tab \eqn{n\times 1}{n*1} \tab Static            \cr
#' \code{list}*      \tab \eqn{T} data.frames/matrices\tab Dynamic    \cr
#' }
#' *: With \eqn{n\times 1}{n\times 1} \code{data.frame} or \code{matrix}.
#'
#' Other cases will return with error.
#'
#' \enumerate{
#'  \item A matrix of size \code{n\times T}{n * T} will be consider a dynamic atribute
#'  \item
#' }
#'
#' In the case of the slices index \code{k}, either an
#' integer vector with the positions, a character vector with the labels of the
#' time periods or a logical vector of length \code{T} can be used to specify
#' which slices to retrieve. Likewise, indexing vertices works in the same way
#' with the only difference that, instead of time period labels and a logical vector
#' of length \code{T}, vertices ids labels and a logical vector of length \code{n}
#' should be provided.
#' @name diffnet-index
#'
#' @examples
#'
#' # Creating a random diffusion network ---------------------------------------
#' set.seed(111)
#' graph <- rdiffnet(10,5)
#'
#' # Accessing to a static attribute
#' graph[["real_threshold"]]
#'
#' # Accessing to subsets of the adjacency matrix
#' graph[1,,1]
#' graph[,,1]
#'
#' # Changing values in the adjacency matrix
#' graph[1,,]
#' graph[1,,] <- -5
#' graph[1,,]
#' @export

#' @export
#' @rdname diffnet-index
`[[.diffnet` <- function(x, name) {
    x$vertex.static.attrs[[name]]
}

#' @export
#' @rdname diffnet-index
`[[<-.diffnet` <- function(x, i, j, value) {
  x$vertex.static.attrs[[i]][j] <- value
  x
}

#' @export
#' @rdname diffnet-index
`[.diffnet` <- function(x, i, j, k) {
  # Checking ids
  if (missing(i)) i <- seq_len(x$meta$n)
  if (missing(j)) j <- seq_len(x$meta$n)
  if (missing(k)) k <- seq_len(x$meta$nper)

  lapply(x$graph[k], "[", i=i, j=j, drop=FALSE)
}

#' @export
#' @rdname diffnet-index
`[<-.diffnet` <- function(x, i, j, k, value) {
  # Checking ids
  if (missing(i)) i <- seq_len(x$meta$n)
  if (missing(j)) j <- seq_len(x$meta$n)
  if (missing(k)) k <- seq_len(x$meta$nper)

  for (l in k)
    x$graph[[l]][i,j] <- value

  x
}

