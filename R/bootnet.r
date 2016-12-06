# library(netdiffuseR)

#' Network Bootstrapping
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
"bootnet"

#' @rdname bootnet
#' @export
resample_graph <- function(graph, self=NULL, ...) {

  # Parameters
  n   <- nrow(graph)
  ids <- 1L:n

  # Checking if self or not
  if (!length(self))
    self <- sum(Matrix::diag(graph)) > 0

  # Storing edge values to be resampled
  E <- graph@x
  if (!length(E))
    stop("The graph is empty.")

  # Random sample (with replacement)
  index <- sample(ids, n, TRUE)

  # Sample graph
  graph_i <- graph[index,][,index]

  # Checking self-edges
  if (!self) {

    # Finding repeated ones
    reps  <- vector("list", n)
    for (j in ids) {
      k <- index[j]
      reps[[k]] <- c(reps[[k]], j)
    }

    # Adding samples
    for (r in reps) {

      # If it wasn't drawn more than once
      if (length(r) <= 1) next

      # readline("(pause)")

      # Sampling edges
      for (j in seq_along(r))
        for (k in seq_along(r)[-j])
          graph_i[r[j],r[k]] <- sample(E, 1)
    }
  }

  return(list(graph=graph_i, indicess=index))
}
#
# set.seed(123)
# g <- rgraph_ba(t = 5-1, self=FALSE)
# bootnet(g, mean, 5)
#
# g
#
# set.seed(123)
# G <- methods::as(matrix(0, ncol=4, nrow=4, dimnames = list(1:4,1:4)), "dgCMatrix")
# G[1,2] <- 12
# G[3,2] <- 32
# G[4,2] <- 42
# bootnet_sample(G, mean, 3, self = FALSE)

# G[c(2,4,2,4),][,c(2,4,2,4)]

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
  boot_res$t0 <- statistic(graph)

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
    var_t       = var(colMeans(boot_res$t, na.rm = TRUE)),
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
        "  Bootstrap Std. Errors = ", paste0(sprintf("%12.4f",mean_t), collapse=", "), "\n",
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
