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
#' @param statistic An statistic
#' @param resample.args List. Arguments to be passed to \code{\link{resample_graph}}
#' @param R Number of reps
#' @param x A \code{diffnet_bootnet} class object.
#' @param ... Further arguments passed to the method (see details).
#' @param main Character scalar. Title of the histogram.
#' @param xlab Character scalar. x-axis label.
#' @param breaks Passed to \code{\link{hist}}.
#' @param annotated Logical scalar. When TRUE marks the observed data average and the simulated data average.
#' @export
#' @family Functions for inference
#' @references Snijders, T. A. B., & Borgatti, S. P. (1999). Non-Parametric
#' Standard Errors and Tests for Network Statistics. Connections, 22(2), 1–10.
#' Retrieved from \url{https://insna.org/PDF/Connections/v22/1999_I-2_61-70.pdf}
#' @return A list of class \code{diffnet_bootnet} containing the following:
#' \item{graph}{The graph passed to \code{bootnet}.}
#' \item{p.value}{The resulting p-value of the test (see details).}
#' \item{t0}{The observed value of the statistic.}
#' \item{mean_t}{The average value of the statistic applied to the simulated networks.}
#' \item{R}{Number of simulations.}
#' \item{statistic}{The function \code{statistic} passed to \code{bootnet}.}
#' \item{boot}{A \code{boot} class object as return from the call to \code{boot}.}
#' \item{resample.args}{The list \code{resample.args} passed to \code{bootnet}.}
#' @name bootnet
#' @examples
#' #
#' set.seed(13)
#' g <- rgraph_ba(t=99)
#'
#' ans <- bootnet(g, function(w, ...) length(w@x), R=100)
#' ans
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
        rand <- runif(2)

        # Add accordingly to density
        if (rand[1] <= dens) {
          graph[r[j],r[k]] <- E[floor(rand[1]/dens*m) + 1]
        }

        if (rand[2] <= dens) {
          graph[r[k],r[j]] <- E[floor(rand[2]/dens*m) + 1]
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
                                     only_fill=FALSE, ...) {

  # Parameters
  n   <- nrow(graph)
  ids <- if (!only_fill) {
    1L:n
  } else {
    as.integer(rownames(graph))
  }

  # Checking if self or not
  if (!length(self))
    self <- sum(Matrix::diag(graph)) > 0

  # Storing edge values to be resampled
  E <- graph@x
  if (!length(E))
    stop("The graph is empty.")

  # Sample graph
  if (!only_fill) {
    # Random sample (with replacement)
    index   <-sample(ids, n, TRUE)
    graph_i <- graph[index,][,index]
  } else {
    index   <- ids
    graph_i <- graph
  }

  # Checking self-edges
  if (!self && useR) {
    graph_i <- bootnet_fillselfR(graph_i, index, E)
  } else if (!self & !useR) {
    graph_i <- bootnet_fillself(graph_i, index, E)
  }

  # Adding names
  dimnames(graph_i) <- list(index,index)

  return(graph_i)
}

resample_graph.list <- function(graph, self, useR, ...) {
  t   <- length(graph)
  out <- graph

  # Names
  tn <- names(graph)
  if (!length(tn)) tn <- 1:t
  names(out) <- tn

  # Generating the first resample
  out[[1]] <- resample_graph.dgCMatrix(graph[[1]], self, useR, ...)

  # If a short list, then return its only element
  if (t == 1) return(out)

  # Repeating
  n      <- nnodes(out)
  newids <- as.integer(rownames(out[[1]]))
  for (i in 2:t) {
    dimnames(out[[i]]) <- list(1:n, 1:n)
    out[[i]] <- resample_graph.dgCMatrix(
      out[[i]][newids,][,newids], self, useR, only_fill = TRUE, ...)
  }

  return(out)
}

resample_graph.array <-function(graph, self, useR, ...) {
  n   <- dim(graph)[1]
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
    ans[i] <- 2*min(mean(meansim[,i] < meanobs[i]),
                    mean(meansim[,i] > meanobs[i]))

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
    fn(do.call(resample_graph, resample.args),...)
  }

  # Calling boot
  boot_res <- boot::boot(1, statisticpll, R=R, fn=statistic, resample.args=resample.args,
                         ...)

  # The t0 must be applied with no rewiring!
  boot_res$t0 <- statistic(graph, ...)

  # Calc pval
  # To be conservative, in a two tail test we use the min of the two
  # So, following davidson & mckinnon Confidence intrval section,
  # p(tau) = 2 * min[F(tau), 1-F(tau)]
  p.value <- bootnet_pval(boot_res$t0, boot_res$t)

  # Creating the object
  out <- list(
    graph       = graph,
    p.value     = p.value,
    t0          = boot_res$t0,
    mean_t      = colMeans(boot_res$t, na.rm = TRUE),
    var_t       = apply(boot_res$t, 2, var, na.rm = TRUE),
    R           = R,
    statistic   = statistic,
    boot        = boot_res,
    resample.args = resample.args
  )

  return(structure(out, class="diffnet_bootnet"))
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

  # Checking graph dim
  res         <- args[[1]]
  res$boot    <- do.call(c, lapply(args, "[[", "boot"))
  # res$p.value <- with(res, 2*min(mean(boot$t < boot$t0),
  #                                mean(boot$t > boot$t0)))
  res$p.value <- bootnet_pval(res$boot$t0, res$boot$t)
  res$mean_t  <- colMeans(res$boot$t, na.rm=TRUE)
  res$R       <- res$boot$R

  res
}

#' @export
#' @rdname bootnet
print.diffnet_bootnet <- function(x, ...) {

  with(x,  {
    nsim <- ifelse(!is.na(R), R, 0)

    cat("\nNetwork Bootstrap (Snijders and Borgatti, 1999)\n",
        "# Num. of drawns  : ", formatC(nsim, digits = 0, format = "f", big.mark = ","),"\n",
        "# nodes           : ", formatC(nnodes(x$graph), digits = 0, format = "f", big.mark = ","),"\n",
        "# of time periods : ", formatC(nslices(x$graph), digits = 0, format = "f", big.mark = ","),"\n",
        paste(rep("-",80), collapse=""),"\n",
        "  Observed mean         = ", paste0(sprintf("%12.4f",t0), collapse=", "), "\n",
        "  Bootstrap Std. Errors = ", paste0(sprintf("%12.4f",var_t), collapse=", "), "\n",
        "  p-value               = ", paste0(sprintf("%12.4f", p.value), collapse=", "),"\n",
        sep="")
  })
  invisible(x)
}

#' @export
#' @param b0 Character scalar. When \code{annotated=TRUE}, label for the value of \code{b0}.
#' @param b Character scalar. When \code{annotated=TRUE}, label for the value of \code{b}.
#' @rdname bootnet
hist.diffnet_bootnet <- function(
  x,
  main="Empirical Distribution of Statistic",
  xlab=expression(Values~of~t),
  breaks=20,
  annotated=TRUE,
  b0=expression(atop(plain("") %up% plain("")), t[0]),
  b =expression(atop(plain("") %up% plain("")), t[]),
  ...) {

  out <- hist(x$boot$t,  breaks=breaks, plot=FALSE)
  ran <- range(out$mids)
  if (annotated) {
    mt <- mean(x$boot$t, na.rm=TRUE)
    ran <- range(c(ran, mt, x$boot$t0))
    hist(x$boot$t, breaks=breaks, main=main, xlab=xlab, xlim = ran, ...)
  } else {
    hist(x$boot$t, breaks=breaks, main=main, xlab=xlab, ...)
  }

  # Adding margin note
  if (annotated) {
    mtext(b0, side = 1, at=x$boot$t0)
    mtext(b, side = 1, at=mt)
  }
  invisible(out)
}
