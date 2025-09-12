# library(netdiffuseR)

#' Network Bootstrapping
#'
#' Implements the bootstrapping method described in Snijders and Borgatti (1999).
#' This function is essentially a wrapper of \code{\link[boot:boot]{boot}}.
#'
#' @templateVar undirected FALSE
#' @templateVar self TRUE
#' @templateVar valued FALSE
#' @template graph_template
#' @param statistic A function that returns a vector with the statistic(s) of interest.
#' The first argument must be the graph, and the second argument a vector of indices
#' (see details)
#' @param resample.args List. Arguments to be passed to \code{\link{resample_graph}}
#' @param R Number of reps
#' @param x A \code{diffnet_bootnet} class object.
#' @param ... Further arguments passed to the method (see details).
#' @param main Character scalar. Title of the histogram.
#' @param xlab Character scalar. x-axis label.
#' @param breaks Passed to \code{\link{hist}}.
#' @param annotated Logical scalar. When TRUE marks the observed data average and the simulated data average.
#' @details
#' Just like the \code{boot} function of the \pkg{boot} package, the \code{statistic}
#' that is passed must have as arguments the original data (the graph in this case),
#' and a vector of indicides. In each repetition, the graph that is passed is a
#' resampled version generated as described in Snijders and Borgatti (1999).
#'
#' When \code{self = FALSE}, for pairs of individuals that haven been drawn more than
#' once the algorithm, in particular, \code{resample_graph}, takes care of filling
#' these pseudo autolinks that are not in the diagonal of the network. By default
#' it is assumed that these pseudo-autolinks depend on whether the original graph
#' had any, hence, if the diagonal has any non-zero value the algorithm assumes that
#' \code{self = TRUE}, skiping the 'filling algorithm'. It is important to notice
#' that, in order to preserve the density of the original network, when
#' assigning an edge value to a pair of the form \eqn{(i,i)} (pseudo-autolinks),
#' such is done with probabilty proportional to the density of the network, in
#' other words, before choosing from the existing list of edge values, the
#' algorithm decides whether to set a zero value first.
#'
#' The vector of indices that is passed to \code{statistic}, an integer vector with range
#' 1 to \eqn{n}, corresponds to the drawn sample of nodes, so the user can, for
#' example, use it to get a subset of a \code{data.frame} that will be used with
#' the \code{graph}.
#'
#' @export
#' @family Functions for inference
#' @references Snijders, T. A. B., & Borgatti, S. P. (1999). Non-Parametric
#' Standard Errors and Tests for Network Statistics. Connections, 22(2), 1–10.
#' Retrieved from \url{https://www.stats.ox.ac.uk/~snijders/Snijders_Borgatti.pdf}
#' @return A list of class \code{diffnet_bootnet} containing the following:
#' \item{graph}{The graph passed to \code{bootnet}.}
#' \item{p.value}{The resulting p-value of the test (see details).}
#' \item{t0}{The observed value of the statistic.}
#' \item{mean_t}{The average value of the statistic applied to the simulated networks.}
#' \item{var_t}{A vector of length \code{length(t0)}. Bootstrap variances.}
#' \item{R}{Number of simulations.}
#' \item{statistic}{The function \code{statistic} passed to \code{bootnet}.}
#' \item{boot}{A \code{boot} class object as return from the call to \code{boot}.}
#' \item{resample.args}{The list \code{resample.args} passed to \code{bootnet}.}
#' @name bootnet
#' @examples
#' # Computing edgecount -------------------------------------------------------
#' set.seed(13)
#' g <- rgraph_ba(t=99)
#'
#' ans <- bootnet(g, function(w, ...) length(w@x), R=100)
#' ans
#'
#' # Generating
NULL

bootnet_fillselfR <- function(graph, index, E) {
  n <- nrow(graph)

  # Finding repeated ones
  reps  <- vector("list", n)
  for (i in 1:n) {
    reps[[index[i]]] <- c(reps[[index[i]]], i)
  }

  # Adding samples
  m    <- length(E)
  dens <- length(E)/(n*n)
  for (r in reps) {

    # If it wasn't drawn more than once
    if (length(r) < 2) next

    # Sampling edges
    for (j in 1:length(r)) {
      for (k in j:length(r)) {

        if (k == j) next

        # Add accordingly to density
        if (runif(1) <= dens) {
          graph[r[j],r[k]] <- E[floor(runif(1)*m) + 1]
        }

        if (runif(1) <= dens) {
          graph[r[k],r[j]] <- E[floor(runif(1)*m) + 1]
        }
      }

    }
  }

  return(graph)
}

#' @rdname bootnet
#' @param useR Logical scalar. When \code{TRUE}, autolinks are filled using an
#' \R based rutine. Otherwise it uses the \pkg{Rcpp} implementation (default).
#' This is intended for testing only.
#' @export
resample_graph <- function(graph, self=NULL, useR=FALSE,...) {
  # Checking undirected (if exists)
  checkingUndirected(graph)

  cls <- class(graph)
  out <- if ("dgCMatrix" %in% cls) {
    resample_graph.dgCMatrix(graph, self, useR, ...)
  } else if ("list" %in% cls) {
    resample_graph.list(graph, self, useR, ...)
  } else if ("matrix" %in% cls) {
    resample_graph.dgCMatrix(methods::as(graph, "dgCMatrix"), self, useR, ...)
  } else if ("diffnet" %in% cls) {
    resample_graph.list(graph$graph, self, useR, ...)
  } else if ("array" %in% cls) {
    resample_graph.array(graph, self, useR, ...)
  } else stopifnot_graph(graph)

  attr(out, "undirected") <- FALSE

  return(out)
}

resample_graph.dgCMatrix <- function(graph, self=NULL, useR=FALSE,
                                     index=NULL, ...) {

  # Checking if self or not
  if (!length(self))
    self <- sum(Matrix::diag(graph)) > 0

  # Storing edge values to be resampled
  E <- graph@x
  if (!length(E))
    stop("The graph is empty.")

  # Sample graph
  if (!length(index)) # Random sample (with replacement)
    index   <-sample(1L:nrow(graph), nnodes(graph), TRUE)

  graph_i <- graph[index,][,index]

  # Getting ids
  original_ids <- dimnames(graph_i)

  # Checking self-edges
  if (!self && useR) {
    graph_i <- bootnet_fillselfR(graph_i, index, E)
  } else if (!self & !useR) {
    graph_i <- bootnet_fillself(graph_i, index, E)
  }

  # Adding names
  attr(graph_i, "sample_indices") <- index
  dimnames(graph_i) <- original_ids

  return(graph_i)
}

resample_graph.list <- function(graph, self, useR, ...) {
  t   <- length(graph)
  out <- graph

  # Names
  tn <- names(graph)
  if (!length(tn)) tn <- 1:t
  names(out) <- tn

  # Generating the first resample (since its dynamic, the sample_indices
  # attribute must be outside)
  out[[1]] <- resample_graph.dgCMatrix(graph[[1]], self, useR, ...)
  attr(out, "sample_indices") <- attr(out[[1]], "sample_indices")

  # If a short list, then return its only element
  if (t == 1) return(out)

  # Repeating
  newids <- attr(out, "sample_indices")
  for (i in 2:t) {
    # Filling blanks and removing sample_indices
    out[[i]] <- resample_graph.dgCMatrix(out[[i]], self, useR, newids, ...)
    attr(out[[i]], "sample_indices") <- NULL
  }

  return(out)
}

resample_graph.array <-function(graph, self, useR, ...) {
  t   <- dim(graph)[3]
  out <- apply(graph, 3, methods::as, Class="dgCMatrix")

  # Checking time names
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t
  names(out) <- tn

  return(resample_graph.list(out, self, useR, ...))
}

bootnet_pval <- function(meanobs, meansim) {
  n <- length(meanobs)
  ans <- vector("numeric",n)
  for (i in 1:n)
    ans[i] <- 2*min(mean(meansim[,i] < 0),
                    mean(meansim[,i] > 0))

  ans
}

#' @rdname bootnet
#' @export
bootnet <- function(
  graph,
  statistic,
  R,
  resample.args = list(self=FALSE),
  ...
) {

  # # Checking class
  # if (!inherits(graph, "diffnet"))
  #   stop("-graph- must be of class diffnet.")

  # Preparing the call to boot
  resample.args$graph <- graph
  statisticpll <- function(d, i, fn, resample.args, ...) {
    g <- do.call(resample_graph, resample.args)
    fn(g, attr(g, "sample_indices"),...)
  }

  # Calling boot
  boot_res <- boot::boot(1, statisticpll, R=R, fn=statistic, resample.args=resample.args,
                         ...)

  # The t0 must be applied with no rewiring!
  boot_res$t0 <- statistic(graph, 1:nnodes(graph), ...)

  # Calc pval
  # To be conservative, in a two tail test we use the min of the two
  # So, following davidson & mckinnon Confidence intrval section,
  # p(tau) = 2 * min[F(tau), 1-F(tau)]
  p.value <- bootnet_pval(boot_res$t0, boot_res$t)

  # Creating the object
  structure(list(
    graph       = graph,
    p.value     = p.value,
    t0          = boot_res$t0,
    mean_t      = colMeans(boot_res$t, na.rm = TRUE),
    var_t       = apply(boot_res$t, 2, var, na.rm=TRUE),
    R           = R,
    statistic   = statistic,
    boot        = boot_res,
    resample.args = resample.args
  ), class="diffnet_bootnet")

}
#
# #' @export
# #' @rdname bootnet
# boot <- bootgraph

#' @export
#' @param recursive Ignored
#' @rdname bootnet
c.diffnet_bootnet <- function(..., recursive=FALSE) {
  # Checking arguments names
  args <- list(...)
  nm <- lapply(args, names)
  if (!all(sapply(nm, function(x) identical(x, nm[[1]]))))
    stop("arguments are not all the same type of \"diffnet_bootnet\" object")

  # Checking

  # Checking graph dim
  res         <- args[[1]]
  res$boot    <- do.call(c, lapply(args, "[[", "boot"))
  # res$p.value <- with(res, 2*min(mean(boot$t < boot$t0),
  #                                mean(boot$t > boot$t0)))
  res$p.value <- bootnet_pval(res$boot$t0, res$boot$t)
  res$mean_t  <- colMeans(res$boot$t, na.rm=TRUE)
  res$var_t   <- apply(res$boot$t, 2, var, na.rm=TRUE)
  res$R       <- res$boot$R

  res
}

#' @export
#' @rdname bootnet
print.diffnet_bootnet <- function(x, ...) {

  # Neat column printing
  netcol <- function(obs, bias, sd, pval) {
    txt <- paste0(sprintf("%10.4f  %10.4f  %10.4f  %10.4f", obs, bias, sd, pval), collapse="\n")
    paste(
      sprintf("%10s  %10s  %10s  %10s", "original", "bias", "std. error", "p.val"),
      txt, sep = "\n"
    )
  }

  with(x,  {
    nsim <- ifelse(!is.na(R), R, 0)

    cat("\nNetwork Bootstrap (Snijders and Borgatti, 1999)\n",
        "# Num. of draws   : ", formatC(nsim, digits = 0, format = "f", big.mark = ","),"\n",
        "# nodes           : ", formatC(nnodes(x$graph), digits = 0, format = "f", big.mark = ","),"\n",
        "# of time periods : ", formatC(nslices(x$graph), digits = 0, format = "f", big.mark = ","),"\n",
        paste(rep("-",options()$width), collapse=""),"\n",
        netcol(t0, t0-apply(boot$t, 2, mean), sqrt(var_t), p.value),"\n",
        sep="")
  })
  invisible(x)
}

#' @export
#' @param b0 Character scalar. When \code{annotated=TRUE}, label for the value of \code{b0}.
#' @param b Character scalar. When \code{annotated=TRUE}, label for the value of \code{b}.
#' @param ask Logical scalar. When \code{TRUE}, asks the user to type \code{<Enter>} to
#' see each plot (as many as statistics where computed).
#' @rdname bootnet
hist.diffnet_bootnet <- function(
  x,
  main      = "Empirical Distribution of Statistic",
  xlab      = expression(Values~of~t),
  breaks    = 20,
  annotated = TRUE,
  b0        = expression(atop(plain("") %up% plain("")), t[0]),
  b         = expression(atop(plain("") %up% plain("")), t[]),
  ask       = TRUE,
  ...) {

  if (ask) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(ask=ask)
  }

  for (i in 1:ncol(x$boot$t)) {
    out <- hist(x$boot$t[,i],  breaks=breaks, plot=FALSE)
    ran <- range(out$mids)

    if (annotated) {
      mt <- mean(x$boot$t[,i], na.rm=TRUE)
      ran <- range(c(ran, mt, x$boot$t0[i]))
      hist(x$boot$t[,i], breaks=breaks, main=main, xlab=xlab, xlim = ran, ...)
    } else {
      hist(x$boot$t[,i], breaks=breaks, main=main, xlab=xlab, ...)
    }

    # Adding margin note
    if (annotated) {
      mtext(b0, side = 1, at=x$boot$t0[i])
      mtext(b, side = 1, at=mt)
    }
  }

  invisible(out)
}

#' @export
#' @param y Ignored.
#' @rdname bootnet
#' @details The `plot.diffnet_bootnet` method is a wrapper for the
#' `hist` method.
plot.diffnet_bootnet <- function(x, y, ...) {
  hist.diffnet_bootnet(x = x, ...)
}
