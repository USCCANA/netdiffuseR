#' Structure dependence test
#'
#' Test whether or not a network estimates can be considered as structure dependent, i.e.
#' a function of the network structure. By rewiring the graph and calculating
#' a particular statistic \eqn{t}, the test compares the observed mean of \eqn{t}
#' against the empirical distribution of it obtained from rewiring the network.
#'
#' @param graph A \code{\link{diffnet}} graph.
#' @param statistic A function that returns either a scalar or a vector.
#' @param R Integer scalar. Number of repetitions.
#' @param rewire.args List. Arguments to be passed to \code{\link{rewire_graph}}
#' @param ncpus Integer scalar. Number of CPU, passed to \code{\link[parallel:makeCluster]{makeCluster}}
#' from the \pkg{parallel} package.
#' @param cl An object of class \code{c("SOCKcluster", "cluster")}
#' @param x A \code{diffnet_boot} class object.
#' @param ... Ignored in the \code{print} method. Otherwise, arguments passed to \code{\link{hist}}.
#' @param main Character scalar. Title of the histogram.
#' @param xlab Character scalar. x-axis label.
#' @param breaks Passed to \code{\link{hist}}.
#' @param annotated Logical scalar. When TRUE marks the observed data average and the simulated data average.
#' @return An object of class \code{diffnet_bot}. As in \code{boot}
#' function of the \pkg{boot} package, a list with the following elements:
#' \item{t0}{The observed value of statistic applied to data.}
#' \item{t}{A matrix with R rows each of which is a bootstrap replicate of the result of calling statistic.}
#' \item{R}{The value of R as passed to \code{boot_net}.}
#' \item{graph}{The graph passed to \code{boot_net}.}
#' \item{seed}{The value of \code{.Random.seed} wheb \code{boot_net} started to work.}
#' \item{statistic}{The function \code{statistic} passed to \code{boot_net}}
#'
#' The output from the \code{hist} method is the same as \code{\link{hist.default}}.
#' @details
#' From the \code{print} method, p-value for the null of the statistic been
#' equal between graph and its rewired versions is computed following Davidson
#' & MacKinnon
#'
#' \deqn{p(\tau)=2\times\min\left(EDF(\tau), 1-EDF(\tau)\right)}{p(t) = 2*min[EDF(t), 1-EDF(t)]}
#'
#' Where \eqn{EDF} is the Empirical Distribution Function.
#' @export
#' @references
#' On development.
#' @author Vega Yon
boot_net <- function(
  graph,
  statistic,
  R,
  rewire.args=list(p=1, undirected=getOption("diffnet.undirected")),
  ncpus=getOption("boot.ncpus", 1L),
  cl=NULL) {

  # Checking class
  if (!inherits(graph, "diffnet"))
    stop("-graph- must be of class diffnet.")

  # parallel setup
  envir <- environment()
  curseed <- .Random.seed

  if (!length(cl)) cl <- parallel::makeCluster(ncpus)
  parallel::clusterExport(cl, c("graph", "ncpus", "R", "rewire.args"),
                          envir = envir)
  on.exit(parallel::stopCluster(cl))

  rewired_graphs <- parallel::clusterEvalQ(cl, {
    library(netdiffuseR)
    with(rewire.args, lapply(1:(R/ncpus), rewire_graph,
                             graph=graph, p=p, undirected=undirected))
  })

  # Coercing into a single list and computing thresholds
  rewired_graphs <- unlist(rewired_graphs, recursive = FALSE)
  rewired_stat   <- sapply(rewired_graphs, statistic)

  if (is.vector(rewired_stat)) rewired_stat <- cbind(rewired_stat)

  # Creating the object
  out <- list(
    t0 = statistic(graph),
    t  = rewired_stat,
    R  = R,
    graph = graph,
    seed = curseed,
    statistic = statistic
  )

  return(structure(out, class="diffnet_boot"))
}

#' @export
#' @rdname boot_net
print.diffnet_boot <- function(x, ...) {
  with(x,  {
    tmean <- colMeans(t, na.rm = TRUE)

    # Calc pval
    # To be conservative, in a two tail test we use the min of the two
    # So, following davidson & mckinnon Confidence intrval section,
    # p(tau) = 2 * min[F(tau), 1-F(tau)]
    test <- 2*min(mean(t < t0), mean(t > t0))
    cat("Network Rewiring graph (",nrow(t),"simulations)\n",paste(rep("-",80), collapse=""),"\n",
        " H0: t - t0 = 0\n",
        " t0      = ", t0, "\n",
        " t       = ", tmean, "\n",
        " p-value = ", sprintf("%.5f",test), sep="")
  })
  invisible(x)
}

#' @export
#' @rdname boot_net
hist.diffnet_boot <- function(
  x,
  main="Distribution of Statistic on rewired network",
  xlab="Values of t",
  breaks=20,
  annotated=TRUE,
  ...) {
  out <- hist(x$t, main=main, xlab=xlab, breaks=breaks, ...)

  # Adding margin note
  if (annotated) {
    mt <- mean(x$t, na.rm=TRUE)
    mtext(expression(atop(plain("") %up% plain("")), t[0]),
          side = 1, at=x$t0)
    mtext(expression(atop(plain("") %up% plain("")), t[]),
          side = 1, at=mt)
  }
  out
}
# #' @rdname boot_net
# boot_thr <- boot_net(net, function(x) {
#   t <- threshold(x)
#   cbind(mean=mean(t, na.rm=TRUE), sd=sd(t, na.rm = TRUE))
# }
# , nsim)
#
# ttest<-(boot_thr$t[1,] - mean(threshold(net), na.rm=TRUE))/boot_thr$t[2,]
