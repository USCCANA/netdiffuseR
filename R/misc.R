#' Recodes an edgelist such that ids go from 1 to n
#' @param data Edgelist as either a matrix or dataframe with ego and alter
#' @param ... Further arguments for the method (ignored)
#' @return A recoded edgelist
#' @export
#' @details Recomended for ease of use
#' @examples
#' edgelist <- cbind(c(1,1,3,6),c(4,3,200,1))
#' edgelist
#' recode(edgelist)
recode <- function(data, ...) UseMethod("recode")

#' @describeIn recode Method for recoding data frames. It keeps the class of the
#' object as it returns a data.frame.
#' @export
recode.data.frame <- function(data, ...) {
  cn <- colnames(data)
  data <- as.data.frame(recode.matrix(as.matrix(data), ...))
  colnames(data) <- cn
  data
}

#' @describeIn recode Method for matrices. Likewise the method for data frames,
#' it returns an object of the same class as the input, a matrix.
#' @export
recode.matrix <- function(data, ...) {

  # Checking the size of the matrixRcppArmadilloForward.h

  data <- as.factor(as.vector(data))
  n <- length(data)
  cbind(data[1:(n/2)], data[(n/2+1):n])
}
