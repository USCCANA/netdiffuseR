n <- 4  # Number of nodes
t <- 3  # Number of time steps
q <- 2  # pathogens

# NEW:      Graph -> n x n
#           attrs -> n x T
#           cumadopt-> n x T x q
#           ans   -> n x T
#           out   -> n x q x T

graph_array <- array(c(
  # First time slice (graph1)
  c(0, 1, 0, 0,
    1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0),
  # Second time slice (graph2)
  c(0, 1, 1, 0,
    1, 0, 0, 1,
    1, 0, 0 ,1,
    0 ,1 ,1 ,0),
  # Third time slice (graph3)
  c(0 ,0 ,1 ,1,
    0 ,0 ,1 ,1,
    1 ,1 ,0 ,0,
    1 ,1 ,0 ,0)),
  dim = c(n,n ,t))

graph <- as_spmat(graph_array)

cumadopt <- array(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,11,12,
                    # Second contagion
                    2,4,6,
                    8,10,12,
                    14,16,18,
                    20,22,24), dim = c(n,t,q))

attrs <- matrix(c(10, 20, 30,
                  40, 50, 60,
                  70, 80, 90,
                  100,110,120), nrow = n)

# Toy model of .exposure

.exposure <- function(graph_slice, cumadopt_slice, attrs_slice,
                      outgoing = TRUE, valued = TRUE, normalized = FALSE, self = FALSE) {

  norm <- graph_slice %*% attrs_slice + 1e-20

  if (!is.na(dim(cumadopt)[3])) {
    ans <- array(0, dim = c(dim(cumadopt)[1],dim(cumadopt)[3]))

    for (q in 1:dim(cumadopt)[3]) {
      if (normalized) {
        ans[,q] <- as.vector(graph_slice %*% (attrs_slice * cumadopt_slice[,,k]) / norm)
      } else {
        ans[,q] <- as.vector(graph_slice %*% (attrs_slice * cumadopt_slice[,,k]))
      }
    }
  } else {
    ans <- graph_slice %*% (attrs_slice * cumadopt_slice)

    if (normalized) {
      ans <- ans/ norm
    }
  }

  #as.vector(ans)
  return(as.vector(ans))
}

# for static graphs it returns `1L`
# nslices  --> from diffnet-methods

lags = 0

if (!is.na(dim(cumadopt)[3])) {
  out <- array(NA, dim = c(dim(cumadopt)[1], dim(cumadopt)[2], dim(cumadopt)[3]))

  if (lags >= 0L) {
    for (i in 1:(nslices(graph) - lags)) {
      out[, i + lags, ] <- .exposure(graph[[i]],
                                     cumadopt[, i, , drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = TRUE,
                                     valued = TRUE,
                                     normalized = FALSE,
                                     self = FALSE)
    }
  } else {
    for (i in (1 - lags):nslices(graph)) {
      out[, i + lags, ] <- .exposure(graph[[i]],
                                     cumadopt[, i, , drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = TRUE,
                                     valued = TRUE,
                                     normalized = FALSE,
                                     self = FALSE)
    }
  }
} else {
  out <- array(NA, dim = c(dim(cumadopt)[1], dim(cumadopt)[2]))

  if (lags >= 0L) {
    for (i in 1:(nslices(graph) - lags)) {
      out[, i + lags] <- .exposure(graph[[i]],
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
                                     cumadopt[, i, drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = TRUE,
                                     valued = TRUE,
                                     normalized = FALSE,
                                     self = FALSE)
    }
  }
}

# out it's working perfectly for multiple diff processes!!
#
# --
#
# --
#
# Now what if there's only one diff process:

# 1 pathogen ONLY

cumadopt <- array(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,11,12
                    ), dim = c(n,t))

attrs <- matrix(c(10, 20, 30,
                  40, 50, 60,
                  70, 80, 90,
                  100,110,120), nrow = n)

if (!is.na(dim(cumadopt)[3])) {
  out <- array(NA, dim = c(dim(cumadopt)[1], dim(cumadopt)[2], dim(cumadopt)[3]))

  if (lags >= 0L) {
    for (i in 1:(nslices(graph) - lags)) {
      out[, i + lags, ] <- .exposure(graph[[i]],
                                     cumadopt[, i, , drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = TRUE,
                                     valued = TRUE,
                                     normalized = FALSE,
                                     self = FALSE)
    }
  } else {
    for (i in (1 - lags):nslices(graph)) {
      out[, i + lags, ] <- .exposure(graph[[i]],
                                     cumadopt[, i, , drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = TRUE,
                                     valued = TRUE,
                                     normalized = FALSE,
                                     self = FALSE)
    }
  }
} else {
  out <- array(NA, dim = c(dim(cumadopt)[1], dim(cumadopt)[2]))

  if (lags >= 0L) {
    for (i in 1:(nslices(graph) - lags)) {
      out[, i + lags] <- .exposure(graph[[i]],
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
                                   cumadopt[, i, drop = FALSE],
                                   attrs[, i, drop = FALSE],
                                   outgoing = TRUE,
                                   valued = TRUE,
                                   normalized = FALSE,
                                   self = FALSE)
    }
  }
}

# Works well..
