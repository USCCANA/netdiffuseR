#' Recodes an edgelist such that ids go from 1 to n
#' @param data Edgelist as either a matrix or dataframe with ego and alter
#' @param ... Further arguments for the method (ignored)
#' @return A recoded edgelist as a two-column matrix/data.frame depending
#' on the class of \code{data}. The output includes an attribute called "recode"
#' which contains a two column data.frame providing a mapping between the
#' previous code and the new code (see the examples)
#' @export
#' @details Required for using most of the package's functions, as ids are used
#' as a reference for accessing elements in adjacency matrices.
#' @seealso \code{\link{edgelist_to_adjmat}}
#' @examples
#' # Simple example
#' edgelist <- cbind(c(1,1,3,6),c(4,3,200,1))
#' edgelist
#' recoded_edgelist <- recode(edgelist)
#' recoded_edgelist
#'
#' # Retrieving the "recode" attribute
#' attr(recoded_edgelist, "recode")
#' @keywords misc
recode <- function(data, ...) UseMethod("recode")

#' @rdname recode
#' @export
recode.data.frame <- function(data, ...) {
  cn <- colnames(data)
  data <- recode.matrix(as.matrix(data), ...)
  output <- as.data.frame(data)
  colnames(data) <- cn
  attr(output, "recode") <- attr(data, "recode")
  output
}

#' @rdname recode
#' @export
recode.matrix <- function(data, ...) {

  # Checking the size of the matrix
  data <- as.factor(as.vector(data))
  n <- length(data)
  output <- cbind(data[1:(n/2)], data[(n/2+1):n])
  data <- unique(data)
  attr(output, "recode") <- data.frame(
    code=as.integer(data), label=as.character(data),
    stringsAsFactors = FALSE)
  output
}
