#' Bootstrapping
#' @param graph A \code{\link{diffnet}} graph.
#' @param statistic A function that returns either a scalar or a vector.
#' @param R Integer scalar. Number of repetitions.
#' @param rewire.args List. Arguments to be passed to \code{\link{rewire_graph}}
#' @param ncpus Integer scalar. Number of CPU, passed to \code{\link[parallel:makeCluster]{makeCluster}}
#' from the \pkg{parallel} package.
#' @param cl An object of class \code{c("SOCKcluster", "cluster")}
#' @param x A \code{diffnet_boot} class object.
#' @param ... Ignored.
#' @return An object of class \code{diffnet_bot}. As in \code{boot}
#' function of the \pkg{boot} package, a list with the following elements:
#' \item{t0}{The observed value of statistic applied to data.}
#' \item{t}{A matrix with R rows each of which is a bootstrap replicate of the result of calling statistic.}
#' \item{R}{The value of R as passed to \code{boot_net}.}
#' \item{graph}{The graph passed to \code{boot_net}.}
#' \item{seed}{The value of \code{.Random.seed} wheb \code{boot_net} started to work.}
#' \item{statistic}{The function \code{statistic} passed to \code{boot_net}}
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
    test  <- sum(t < t0)/R
    if (test > .5) test  <- 1-test
    cat("boot graph\n",
        "t0         = ", t0, "\n",
        "t          = ", tmean, "\n",
        "P(t != t0) = ", sprintf("%.5f",test))
  })
  invisible(x)
}
