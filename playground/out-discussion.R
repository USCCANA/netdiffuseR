n <- 4  # Number of nodes
t <- 3  # Number of time steps

# static graph
graph <- matrix(c(0, 1, 0, 0,
                  1, 0, 1, 0,
                  0, 1, 0, 1,
                  0, 0, 1 ,0), nrow = n)

cumadopt <- matrix(c(1, 0, 0,
                     1, 1, 0,
                     0, 1, 1,
                     0, 0, 1), nrow = n)

attrs <- matrix(c(10, 20, 30,
                  40, 50, 60,
                  70, 80, 90,
                  100,110,120), nrow = n)

# Toy model of .exposure
.exposure <- function(graph, cumadopt, attrs,
                      outgoing = TRUE, valued = TRUE, normalized = FALSE, self = FALSE) {

  ans <- ( graph %*% (attrs * cumadopt) )
  ans_norm <- ans/( graph %*% attrs + 1e-20 )
  return(as.vector(ans/ans_norm))
}

# for static graphs it returns `1L`
nslices <- function(graph) {
  if ("matrix" %in% class(graph)) {
    return(1L)
  } else if ("list" %in% class(graph)) {
    return(length(graph))
  }
}

out <- matrix(nrow = nrow(cumadopt), ncol = ncol(cumadopt))

lags <- 0

if (lags >= 0L) {
  for (i in 1:(nslices(graph) - lags)) {
    out[, i + lags] <- .exposure(graph,
                                 cumadopt[, i, drop = FALSE],
                                 attrs[, i, drop = FALSE],
                                 outgoing = TRUE,
                                 valued = TRUE,
                                 normalized = FALSE,
                                 self = FALSE)
  }
} else {
  for (i in (1 - lags):nslices(graph)) {
    out[, i + lags] <- .exposure(graph[[i]],
                                 cumadopt[, i],
                                 attrs[, i],
                                 outgoing = TRUE,
                                 valued = TRUE,
                                 normalized = FALSE,
                                 self = FALSE)
  }
}

# View the result of 'out'
print(out)
