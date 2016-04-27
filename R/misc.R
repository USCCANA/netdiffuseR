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
#' @author George G. Vega Yon
recode <- function(data, ...) UseMethod("recode")

#' @rdname recode
#' @export
recode.data.frame <- function(data, ...) {
  dn <- dimnames(data)
  data <- recode.matrix(as.matrix(data), ...)
  output <- as.data.frame(data)
  dimnames(output) <- dn
  attr(output, "recode") <- attr(data, "recode")
  output
}

#' @rdname recode
#' @export
recode.matrix <- function(data, ...) {

  # Checking the size of the matrix
  dn <- dimnames(data)
  data <- as.factor(as.character(as.vector(data)))
  n <- length(data)
  output <- cbind(data[1:(n/2)], data[(n/2+1):n])
  data <- unique(data)

  # Previous order w/ codes
  rc <- data.frame(
    code=as.integer(data), label=as.character(data),
    stringsAsFactors = FALSE)

  rc <- rc[order(rc[,1]),]

  # Removing NA
  rc <- rc[!is.na(rc$code),]

  dimnames(output) <- dn

  attr(output, "recode") <- rc
  output
}
