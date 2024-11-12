n <- 4  # Number of nodes
t <- 3  # Number of time steps
q <- 2  # pathogens

# NEW:      Graph -> n x n
#           attrs -> n x T
#           cumadopt-> n x T x q
#           ans   -> n x T
#           out   -> n x T x q

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

# one
cumadopt_one <- matrix(c(
  0, 1, 1,
  1, 1, 1,
  0, 0, 1,
  0, 0, 0
), nrow = n, byrow = TRUE)

#two
cumadopt_two <- array(0, dim = c(n, t, q))

cumadopt_two[,,1] <- matrix(c(
  0, 1, 1,
  1, 1, 1,
  0, 0, 1,
  0, 0, 0
), nrow = n, byrow = TRUE)

cumadopt_two[,,2] <- matrix(c(
  0, 1, 1,
  0, 1, 1,
  0, 0, 1,
  0, 0, 1
), nrow = n, byrow = TRUE)

# attributes between [0-1]
attrs <- matrix(runif(n * t), nrow = n)


# Toy model of .exposure
.exposure <- function(graph_slice, cumadopt_slice, attrs_slice,
                      outgoing = TRUE, valued = TRUE, normalized = FALSE, self = FALSE) {

  norm <- graph_slice %*% attrs_slice + 1e-20

  if (!is.na(dim(cumadopt)[3])) {
    ans <- array(0, dim = c(dim(cumadopt)[1],dim(cumadopt)[3]))

    for (q in 1:dim(cumadopt)[3]) {
      if (normalized) {
        ans[,q] <- as.vector(graph_slice %*% (attrs_slice * cumadopt_slice[,,q]) / norm)
      } else {
        ans[,q] <- as.vector(graph_slice %*% (attrs_slice * cumadopt_slice[,,q]))
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

cumadopt <- cumadopt_two

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

cumadopt <- cumadopt_one

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
