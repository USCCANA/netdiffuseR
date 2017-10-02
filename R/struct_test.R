#' Structure dependence test
#'
#' Test whether or not a network estimates can be considered structurally dependent, i.e.
#' a function of the network structure. By rewiring the graph and calculating
#' a particular statistic \eqn{t}, the test compares the observed mean of \eqn{t}
#' against the empirical distribution of it obtained from rewiring the network.
#'
#' @param graph A \code{\link{diffnet}} graph.
#' @param statistic A function that returns either a scalar or a vector.
#' @param R Integer scalar. Number of repetitions.
#' @param rewire.args List. Arguments to be passed to \code{\link{rewire_graph}}
#' @param x A \code{diffnet_struct_test} class object.
#' @param ... Further arguments passed to the method (see details).
#' @param main Character scalar. Title of the histogram.
#' @param xlab Character scalar. x-axis label.
#' @param breaks Passed to \code{\link{hist}}.
#' @param annotated Logical scalar. When TRUE marks the observed data average and the simulated data average.
#' @return A list of class \code{diffnet_struct_test} containing the following:
#' \item{graph}{The graph passed to \code{struct_test}.}
#' \item{p.value}{The resulting p-value of the test (see details).}
#' \item{t0}{The observed value of the statistic.}
#' \item{mean_t}{The average value of the statistic applied to the simulated networks.}
#' \item{R}{Number of simulations.}
#' \item{statistic}{The function \code{statistic} passed to \code{struct_test}.}
#' \item{boot}{A \code{boot} class object as return from the call to \code{boot}.}
#' \item{rewire.args}{The list \code{rewire.args} passed to \code{struct_test}.}
#'
#' @details
#'
#' \code{struct_test} computes the test by generating the null distribution using
#' Monte Carlo simulations (rewiring). \code{struct_test_asymp} computes the
#' test using an asymptotic approximation. While available, we do not recommend
#' using the asymptotic approximation since it has not shown good results when
#' compared to the MC approximation. Furthermore, the asymptotic version has only
#' been implemented for \code{graph} as static graph.
#'
#' The output from the \code{hist} method is the same as \code{\link{hist.default}}.
#'
#' \code{struct_test} is a wrapper for the function \code{\link[boot:boot]{boot}} from the
#' \pkg{boot} package. Instead of resampling data--vertices or edges--in each iteration the function
#' rewires the original graph using \code{\link{rewire_graph}} and applies
#' the function defined by the user in \code{statistic}.
#'
#' The default values to \code{rewire_graph} via \code{rewire.args} are:
#' \tabular{ll}{
#' \code{p}          \tab Number or Integer with default \code{n_rewires(graph)}. \cr
#' \code{undirected} \tab Logical scalar with default \code{getOption("diffnet.undirected", FALSE)}. \cr
#' \code{copy.first} \tab Logical scalar with \code{TRUE}. \cr
#' \code{algorithm}  \tab Character scalar with default \code{"swap"}.
#' }
#'
#' In \code{struct_test} \code{\dots} are passed to \code{boot}, otherwise are passed
#' to the corresponding method (\code{\link{hist}} for instance).
#'
#' From the \code{print} method, p-value for the null of the statistic been
#' equal between graph and its rewired versions is computed as follows
#'
#'
#' \deqn{%
#' p(\tau)=2\times\min\left(\mbox{Pr}(t\leq\tau), \mbox{Pr}(t\geq\tau)\right) %
#' }{ %
#' p(tau) = 2*min[Pr(t<=tau), Pr(t>=tau)] %
#' }
#'
#' Where \eqn{\mbox{Pr}\{\cdot\}}{Pr(.)} is approximated using the
#' Empirical Distribution Function retrieved from the simulations.
#'
#' For the case of the asymptotic approximation, under the null we have
#'
#' \deqn{%
#' \sqrt{n}\left(\hat\beta(Y,G)-\mu_\beta\right)\sim^d\mbox{N}\left(0,\sigma_\beta^2\right)
#' }{%
#' sqrt(n)*[mean(beta) - E[beta]]/Var[beta] ~ N(0,1)
#' }
#'
#'
#' The test is actually on development by Vega Yon and Valente. A copy of the
#' working paper can be distributed upon request to \email{g.vegayon@gmail.com}.
#'
#' The function \code{n_rewires} proposes a vector of number of rewirings that
#' are performed in each iteration.
#'
#' @family Functions for inference
#' @references
#' Vega Yon, George G. and Valente, Thomas W. (On development).
#'
#' Davidson, R., & MacKinnon, J. G. (2004). Econometric Theory and Methods. New York:
#' Oxford University Press.
#' @examples
#' # Creating a random graph
#' set.seed(881)
#' diffnet <- rdiffnet(100, 5, seed.graph="small-world")
#'
#' # Testing structure-dependency of threshold
#' res <- struct_test(diffnet, function(g) mean(threshold(g), na.rm=TRUE), R=100)
#' res
#' hist(res)
#'
#' # Adding a legend
#' legend("topright", bty="n",
#'  legend=c(
#'    expression(t[0]:~Baseline),
#'    expression(t:~Rewired~average)
#'  )
#'  )
#'
#' # Concatenating results
#' c(res, res)
#'
#' # Running in parallel fashion
#' \dontrun{
#' res <- struct_test(diffnet, function(g) mean(threshold(g), na.rm=TRUE), R=100,
#' ncpus=4, parallel="multicore")
#' res
#' hist(res)
#' }
#' @author George G. Vega Yon
#' @name struct_test
NULL

#' @export
#' @param p Either a Numeric scalar or vector of length \code{nslices(graph)-1}
#' with the number of rewires per links.
#' @rdname struct_test
n_rewires <- function(graph, p=c(20L, rep(.1, nslices(graph) - 1))) {
  as.integer(round(unlist(nlinks(graph))*p))
}

struct_test_pval <- function(meanobs, meansim) {
  n <- length(meanobs)
  ans <- vector("numeric",n)
  for (i in 1:n)
    ans[i] <- 2*min(mean(meansim[,i] < meanobs[i]),
                    mean(meansim[,i] > meanobs[i]))

  ans
}

#' @rdname struct_test
#' @export
struct_test <- function(
  graph,
  statistic,
  R,
  rewire.args=list(),
  ...
  ) {

  # Checking defaults
  if (!length(rewire.args$p)) rewire.args$p <- n_rewires(graph)
  if (!length(rewire.args$undirected))
    rewire.args$undirected <- getOption("diffnet.undirected", FALSE)
  if (!length(rewire.args$copy.first)) rewire.args$copy.first <- TRUE
  if (!length(rewire.args$algorithm)) rewire.args$algorithm <- "swap"

  # Preparing the call to boot
  rewire.args$graph <- graph
  statisticpll <- function(d, i, fn, rewire.args, ...) {
    fn(do.call(rewire_graph, rewire.args),...)
  }

  # Calling boot
  boot_res <- boot::boot(1, statisticpll, R=R, fn=statistic, rewire.args=rewire.args,
                   ...)

  # The t0 must be applied with no rewiring!
  boot_res$t0 <- statistic(graph)

  # Calc pval
  # To be conservative, in a two tail test we use the min of the two
  # So, following davidson & mckinnon Confidence intrval section,
  # p(tau) = 2 * min[F(tau), 1-F(tau)]
  p.value <- struct_test_pval(boot_res$t0, boot_res$t)

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
    rewire.args = rewire.args
  )

  return(structure(out, class="diffnet_struct_test"))
}

#' @export
#' @param recursive Ignored
#' @rdname struct_test
c.diffnet_struct_test <- function(..., recursive=FALSE) {
  # Checking arguments names
  args <- list(...)
  nm <- lapply(args, names)
  if (!all(sapply(nm, function(x) identical(x, nm[[1]]))))
    stop("arguments are not all the same type of \"diffnet_struct_test\" object")

  # Checking graph dim
  res         <- args[[1]]
  res$boot    <- do.call(c, lapply(args, "[[", "boot"))
  # res$p.value <- with(res, 2*min(mean(boot$t < boot$t0),
  #                                mean(boot$t > boot$t0)))
  res$p.value <- struct_test_pval(res$boot$t0, res$boot$t)
  res$mean_t  <- colMeans(res$boot$t, na.rm=TRUE)
  res$R       <- res$boot$R

  res
}

#' @export
#' @rdname struct_test
print.diffnet_struct_test <- function(x, ...) {

  # Neat column printing
  netcol <- function(obs, expe, pval) {
    txt <- paste0(sprintf("  %10.4f  %10.4f  %10.4f", obs, expe, pval), collapse="\n")
    paste(
      sprintf("  %10s  %10s  %10s", "observed", "expected", "p.val"),
      txt, sep = "\n"
    )
  }

  with(x,  {
    nsim <- ifelse(!is.na(R), R, 0)

    cat("\nStructure dependence test\n",
        "# Simulations     : ", formatC(nsim, digits = 0, format = "f", big.mark = ","),"\n",
        "# nodes           : ", formatC(nnodes(x$graph), digits = 0, format = "f", big.mark = ","),"\n",
        "# of time periods : ", formatC(nslices(x$graph), digits = 0, format = "f", big.mark = ","),"\n",
        paste(rep("-",80), collapse=""),"\n",
        " H0: E[beta(Y,G)|G] - E[beta(Y,G)] = 0 (no structure dependency)\n",
        netcol(t0, mean_t, p.value),
        sep="")
  })
  invisible(x)
}

#' @export
#' @param b0 Character scalar. When \code{annotated=TRUE}, label for the value of \code{b0}.
#' @param b Character scalar. When \code{annotated=TRUE}, label for the value of \code{b}.
#' @param ask Logical scalar. When \code{TRUE}, asks the user to type \code{<Enter>} to see each plot (as
#' many as statistics where computed).
#' @rdname struct_test
hist.diffnet_struct_test <- function(
  x,
  main="Empirical Distribution of Statistic",
  xlab=expression(Values~of~t),
  breaks=20,
  annotated=TRUE,
  b0=expression(atop(plain("") %up% plain("")), t[0]),
  b =expression(atop(plain("") %up% plain("")), t[]),
  ask = TRUE,
  ...) {

  # Par parameters
  if (ask) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(ask=ask)
  }

  out <- vector("list", ncol(x$boot$t))
  for (i in 1:ncol(x$boot$t)) {

    # Computing the histogram and its range
    out[[i]] <- hist(x$boot$t[,i],  breaks=breaks, plot=FALSE)
    ran      <- range(out[[i]]$mids, na.rm = TRUE)

    if (annotated) {

      # Annotated version (pretty)
      mt <- mean(x$boot$t[,i], na.rm=TRUE)
      ran <- range(c(ran, mt, x$boot$t0[i]), na.rm=TRUE)
      hist(x$boot$t[,i], breaks=breaks, main=main, xlab=xlab, xlim = ran, ...)

    } else {

      # Not annotated version
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
# #' @rdname struct_test
# boot_thr <- struct_test(net, function(x) {
#   t <- threshold(x)
#   cbind(mean=mean(t, na.rm=TRUE), sd=sd(t, na.rm = TRUE))
# }
# , nsim)
#
# ttest<-(boot_thr$t[1,] - mean(threshold(net), na.rm=TRUE))/boot_thr$t[2,]

#' @export
#' @rdname struct_test
#' @param Y Numeric vector of length \eqn{n}.
#' @param statistic_name Character scalar. Name of the metric to compute. Currently
#' this can be either \code{"distance"},\code{">"},\code{"<"},\code{"=="}, \code{">="},
#' or \code{"<="}.
struct_test_asymp <- function(
  graph, Y,
  statistic_name="distance", p=2.0, ...) {

  # Distance metric
  D <- vertex_covariate_compare(graph, Y, statistic_name)

  # Computing observed mean
  m_obs <- mean(Matrix::rowSums(D*graph)/(Matrix::rowSums(graph) + 1e-15), na.rm=TRUE)

  # Checking NA
  isna <- which(is.na(Y))
  if (length(isna)) {
    warning("There are some NA values in Y: ",paste(isna, collapse=", "))
    Y <- Y[-isna]
  }

  # Computing theoreticals
  m_null <- struct_test_mean(Y, statistic_name)
  v_null <- struct_test_var(Y, statistic_name)

  p.value <- pnorm(sqrt(nnodes(graph))*(m_obs-m_null)/sqrt(v_null))
  p.value <- ifelse(p.value > .5, 1-p.value, p.value)*2

  # Creating the object
  out <- list(
    graph       = graph,
    p.value     = p.value,
    t0          = m_obs,
    mean_t      = m_null,
    var_t       = v_null,
    R           = NA,
    statistic   = statistic_name,
    boot        = NA,
    rewire.args = NA
  )

  return(structure(out, class="diffnet_struct_test"))
}

# library(mvtnorm)
#
# rm(list=ls())
# set.seed(123)
# n <- 400
# G <- rgraph_ws(n,8,.2)
# G <- G/(Matrix::rowSums(G) + 1e-15)
# X <- cbind(rchisq(n, 4))
#
# # Case 1
# ape::Moran.I(as.vector(X), as.matrix(G))
# struct_test_approx(G, X, "distance")
#
# # Case 2
# X <- as.vector(
#   solve(Matrix::Diagonal(n) - .4*G) %*% matrix(rmvnorm(1, rep(1,n)), ncol=1)
# )
# ape::Moran.I(as.vector(X), as.matrix(G))
# struct_test_approx(G, X, "distance")

