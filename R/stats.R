#' Indegree, outdegree and degree of the vertices
#'
#' Computes the requested degree measure for each node in the graph.
#'
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @param cmode Character scalar. Either "indegree", "outdegree" or "degree".
#' @param undirected Logical scalar. TRUE when the graph is undirected.
#' @param self Logical scalar.. TRUE when self edges should not be considered.
#' @param valued Logical scalar. When FALSE sets every non-zero entry of \code{graph} to one.
#' @return Either a numeric vector of size \eqn{n}{n} with the degree of each node (if graph is
#' a matrix), or a matrix of size \eqn{n\times T}{n * T}.
#' @export
#' @family statistics
#' @keywords univar
#' @aliases degree indegree outdegree
#' @examples
#'
#' # Comparing degree measurements ---------------------------------------------
#' # Creating an undirected graph
#' graph <- rgraph_ba()
#' graph
#'
#' data.frame(
#'    In=dgr(graph, "indegree", undirected = FALSE),
#'    Out=dgr(graph, "outdegree", undirected = FALSE),
#'    Degree=dgr(graph, "degree", undirected = FALSE)
#'  )
#'
#' # Testing on Korean Family Planning (weighted graph) ------------------------
#' data(kfamilyDiffNet)
#' d_unvalued <- dgr(kfamilyDiffNet, valued=FALSE)
#' d_valued   <- dgr(kfamilyDiffNet, valued=TRUE)
#'
#' any(d_valued!=d_unvalued)
#'
#' @author George G. Vega Yon
dgr <- function(graph, cmode="degree",
                undirected=getOption("diffnet.undirected", FALSE),
                self=getOption("diffnet.self",FALSE),
                valued=getOption("diffnet.valued", FALSE)) {

  switch (class(graph),
    matrix = dgr.matrix(graph, cmode, undirected, self, valued),
    array = dgr.array(graph, cmode, undirected, self, valued),
    dgCMatrix = dgr.dgCMatrix(graph, cmode, undirected, self, valued),
    list = dgr.list(graph, cmode, undirected, self, valued),
    diffnet = dgr.list(graph$graph, cmode, undirected = graph$meta$undirected, self, valued),
    stopifnot_graph(graph)
  )

}

# @rdname dgr
# @export
dgr.matrix <- function(
  graph, cmode, undirected, self, valued) {

  # Checking dimensions
  dm <- dim(graph)
  if (dm[1] != dm[2]) stop("-graph- must be a square matrix.")

  # Retrieving the number
  if      (cmode=="indegree")  cmode <- 0
  else if (cmode=="outdegree") cmode <- 1
  else if (cmode=="degree")    cmode <- 2
  else stop('Invalid -cmode- ',cmode,'. Should be either ',
            '"indegree", "outdegree" or "degree".')

  # Computing degree
  output <- degree_cpp(methods::as(graph, "dgCMatrix"), cmode, undirected, self,
                       valued)

  # Naming
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:nrow(graph)
  rownames(output) <- rn

  output
}

# @rdname dgr
# @export
dgr.dgCMatrix <- function(graph, cmode, undirected, self, valued) {

  # Checking dimensions
  dm <- dim(graph)
  if (dm[1] != dm[2]) stop("-graph- must be a square matrix.")

  # Retrieving the number
  if      (cmode=="indegree")  cmode <- 0
  else if (cmode=="outdegree") cmode <- 1
  else if (cmode=="degree")    cmode <- 2
  else stop('Invalid -cmode- ',cmode,'. Should be either ',
            '"indegree", "outdegree" or "degree".')

  # Computing degree
  output <- degree_cpp(graph, cmode, undirected, self, valued)

  # Naming
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:nrow(graph)
  rownames(output) <- rn

  output
}

# @rdname dgr
# @export
dgr.list <- function(graph, cmode, undirected, self, valued) {
  n <- ncol(graph[[1]])
  t <- length(graph)
  output <- matrix(ncol=t, nrow=n)

  for(i in 1:t)
    output[,i] <- dgr(graph[[i]], cmode, undirected, self, valued)

  # Adding names
  cn <- names(graph)
  if (!length(cn)) cn <- 1:length(graph)
  colnames(output) <- cn

  # Naming
  rn <- rownames(graph[[1]])
  if (!length(rn)) rn <- 1:nrow(graph[[1]])
  rownames(output) <- rn

  output
}

# @rdname dgr
# @export
dgr.array <- function(graph, cmode, undirected, self, valued) {
  n <- dim(graph)[1]
  t <- dim(graph)[3]
  output <- matrix(ncol=t, nrow=n)

  for(i in 1:t)
    output[,i] <- dgr(methods::as(graph[,,i], "dgCMatrix"), cmode, undirected, self, valued)

  # Adding names
  cn <- dimnames(graph)[[3]]
  if (!length(cn)) cn <- 1:dim(graph)[3]
  colnames(output) <- cn

  # Naming
  rn <- dimnames(graph)[[1]]
  if (!length(rn)) rn <- 1:nrow(graph)
  rownames(output) <- rn

  output
}

#' Ego exposure
#'
#' Calculates exposure to adoption over time via multiple different types of weight
#' matrices.  The basic  model is exposure to adoption by immediate neighbors
#' (outdegree) at the time period prior to ego’s adoption. This exposure can also be
#' based on (1) incoming ties, (2) structural equivalence, (3) indirect ties, (4)
#' attribute weighted (5) network-metric weighted (e.g., central nodes have more
#' influence), and attribute-weighted (e.g., based on homophily or tie strength).
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param cumadopt nxT matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @param attrs Either a character scalar (if \code{graph} is diffnet),
#' or a numeric matrix of size \eqn{n\times T}{n * T}. Weighting for each time, period (see details).
#' @param alt.graph Either a dynamic graph that should be used instead of \code{graph},
#' or \code{"se"} (see details).
#' @param outgoing Logical scalar. When \code{TRUE}, computed using outgoing ties.
#' @param valued Logical scalar. When \code{FALSE}, values of \code{graph} are set to one.
#' @param normalized Logical scalar. When \code{TRUE}, the exposure will be between zero
#' and one (see details).
#' @param ... Further arguments passed to \code{\link{struct_equiv}} (only used when
#' \code{alt.graph="se"}).
#' @param groupvar Passed to \code{\link{struct_equiv}}.
#' @details
#' Exposure is calculated as follows:
#'
#' \deqn{ %
#' E_t = \left(S_t \times \left[x_t \circ A_t\right]\right) / (S_t \times x_t) %
#' }{%
#' E(t) = (S(t) \%*\% [x(t) * A(t)]) / [S(t) \%*\% x(t)]
#' }
#'
#' Where \eqn{S_t}{S(t)} is the graph in time \eqn{t}, \eqn{x_t}{x(t)} is an attribute
#' vector of size \eqn{n} at time \eqn{t}, \eqn{A_t}{A(t)} is the t-th column of
#' the cumulative adopters matrix (a vector of length \eqn{n} with \eqn{a_{ti}=1}{a(t,i)=1}
#' if \eqn{i} has adopted at or prior to \eqn{t}), \eqn{\circ}{*} is the kronecker
#' product (element-wise), and \eqn{\times}{\%*\%} is the matrix product.
#'
#' By default the graph used for this calculation, \eqn{S}, is the social network. Alternatively,
#' in the case of \code{diffnet} objects, the user can provide an alternative
#' graph using \code{alt.graph}. An example of this would be using \eqn{1/SE},
#' the element-wise inverse of the structural equivalence matrix (see example below).
#' Furthermore, if \code{alt.graph="se"}, the inverse of the structural equivalence
#' is computed via \code{\link{struct_equiv}} and used instead of the provided
#' graph. Notice that when using a valued graph the option \code{valued} should
#' be equal to \code{TRUE}, this check is run automatically when running the
#' model using structural equivalence.
#'
#' \bold{An important remark} is that when calculating \bold{structural equivalence} the
#' function \bold{assumes that this is to be done to the entire graph} regardless of
#' disconnected communities (as in the case of the medical innovations
#' data set). Hence, structural equivalence for individuals for two different
#' communites may not be zero. If the user wants to calculate structural
#' equivalence separately by community, he should create different diffnet
#' objects and do so (see example below). Alternatively, for the case of
#' diffnet objects, by using the option \code{groupvar} (see \code{\link{struct_equiv}}), the user can provide
#' the function with the name of a grouping variable--which should one in the
#' set of static vertex attributes--so that the algorithm is done by group
#' (or community) instead of in an aggregated way.
#'
#' If the user does not specifies a particular weighting attribute in \code{attrs},
#' the function sets this as a matrix of ones. Otherwise the function will return
#' an attribute weighted exposure. When \code{graph} is of class \code{diffnet},
#' \code{attrs} can be a character scalar specifying the name of any of the graph's
#' attributes, both dynamic and static. See the examples section for a demonstration using
#' degree.
#'
#' When \code{outgoing=FALSE}, \eqn{S} is replaced by its transposed, so in the
#' case of a social network exposure will be computed based on the incomming ties.
#'
#' If \code{normalize=FALSE} then denominator, \eqn{S_t \times x_t}{S(t) \%*\% x(t)},
#' is not included. This can be useful when, for example, exposure needs to be
#' computed as a count instead of a proportion. A good example of this can be
#' found at the examples section of the function \code{\link{rdiffnet}}.
#'
#' @references
#' Burt, R. S. (1987). "Social Contagion and Innovation: Cohesion versus Structural
#' Equivalence". American Journal of Sociology, 92(6), 1287.
#' \url{http://doi.org/10.1086/228667}
#'
#' Valente, T. W. (1995). "Network models of the diffusion of innovations"
#'  (2nd ed.). Cresskill N.J.: Hampton Press.
#'
#' @examples
#' # Calculating the exposure based on Structural Equivalence ------------------
#' set.seed(113132)
#' graph <- rdiffnet(100, 10)
#'
#' SE <- lapply(struct_equiv(graph), "[[", "SE")
#' SE <- lapply(SE, function(x) {
#'    x <- 1/x
#'    x[!is.finite(x)] <- 0
#'    x
#' })
#'
#' # Recall setting valued equal to TRUE!
#' expo_se <- exposure(graph, alt.graph=SE , valued=TRUE)
#'
#' # These three lines are equivalent to:
#' expo_se2 <- exposure(graph, alt.graph="se", valued=TRUE)
#' # Notice that we are setting valued=TRUE, but this is not necesary since when
#' # alt.graph = "se" the function checks this to be setted equal to TRUE
#'
#' # Weighted Exposure using degree --------------------------------------------
#' eDE <- exposure(graph, attrs=dgr(graph))
#'
#' # Which is equivalent to
#' graph[["deg"]] <- dgr(graph)
#' eDE2 <- exposure(graph, attrs="deg")
#'
#' # Comparing using incomming edges -------------------------------------------
#' eIN <- exposure(graph, outgoing=FALSE)
#'
#' # Structral equivalence for different communities ---------------------------
#' data(medInnovationsDiffNet)
#'
#' # METHOD 1: Using the c.diffnet method:
#'
#' # Creating subsets by city
#' cities <- unique(medInnovationsDiffNet[["city"]])
#'
#' diffnet <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == cities[1]]
#' diffnet[["expo_se"]] <- exposure(diffnet, alt.graph="se", valued=TRUE)
#'
#' for (v in cities[-1]) {
#'    diffnet_v <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == v]
#'    diffnet_v[["expo_se"]] <- exposure(diffnet_v, alt.graph="se", valued=TRUE)
#'    diffnet <- c(diffnet, diffnet_v)
#' }
#'
#' # We can set the original order (just in case) of the data
#' diffnet <- diffnet[medInnovationsDiffNet$meta$ids]
#' diffnet
#'
#' # Checking everything is equal
#' test <- summary(medInnovationsDiffNet, no.print=TRUE) ==
#'    summary(diffnet, no.print=TRUE)
#'
#' stopifnot(all(test))
#'
#' # METHOD 2: Using the 'groupvar' argument
#' # Further, we can compare this with using the groupvar
#' diffnet[["expo_se2"]] <- exposure(diffnet, alt.graph="se",
#'    groupvar="city", valued=TRUE)
#'
#' # These should be equivalent
#' test <- diffnet[["expo_se", as.df=TRUE]] == diffnet[["expo_se2", as.df=TRUE]]
#' stopifnot(all(test))
#'
#' # METHOD 3: Computing exposure, rbind and then adding it to the diffnet object
#' expo_se3 <- NULL
#' for (v in unique(cities))
#'    expo_se3 <- rbind(
#'      expo_se3,
#'      exposure(
#'        diffnet[diffnet[["city"]] == v],
#'        alt.graph = "se", valued=TRUE
#'      ))
#'
#' # Just to make sure, we sort the rows
#' expo_se3 <- expo_se3[diffnet$meta$ids,]
#'
#' diffnet[["expo_se3"]] <- expo_se3
#'
#' test <- diffnet[["expo_se", as.df=TRUE]] == diffnet[["expo_se3", as.df=TRUE]]
#' stopifnot(all(test))
#'
#'
#' # METHOD 4: Using the groupvar in struct_equiv
#' se <- struct_equiv(diffnet, groupvar="city")
#' se <- lapply(se, "[[", "SE")
#' se <- lapply(se, function(x) {
#'    x <- 1/x
#'    x[!is.finite(x)] <- 0
#'    x
#' })
#'
#' diffnet[["expo_se4"]] <- exposure(diffnet, alt.graph=se, valued=TRUE)
#'
#' test <- diffnet[["expo_se", as.df=TRUE]] == diffnet[["expo_se4", as.df=TRUE]]
#' stopifnot(all(test))
#'
#'
#'
#' @family statistics
#' @keywords univar
#' @return A matrix of size \eqn{n\times T}{n * T} with exposure for each node.
#' @export
#' @author George G. Vega Yon, Stephanie R. Dyal, Timothy B. Hayes & Thomas W. Valente
exposure <- function(graph, cumadopt, attrs = NULL, alt.graph=NULL,
                     outgoing=getOption("diffnet.outgoing", TRUE),
                     valued=getOption("diffnet.valued", FALSE), normalized=TRUE,
                     groupvar=NULL,
                     ...) {

  # Checking diffnet attributes
  if (length(attrs) == 1 && inherits(attrs, "character")) {
    if (!inherits(graph, "diffnet"))
      stop("Specifying -attrs- as a character scalar is only valid for -diffnet- objects.")

    # Retrieving attribute
    attrs <- graph[[attrs]]

    # Coercing into a matrix
    attrs <- if (inherits(attrs, "list")) do.call(cbind, attrs)
    else matrix(attrs, ncol=nslices(graph), nrow=nvertices(graph))
  }

  # Checking groupvar
  if (length(groupvar) == 1 && inherits(graph, "diffnet"))
    groupvar <- graph[[groupvar]]

  # Checking cumadopt mat
  if (missing(cumadopt))
    if (!inherits(graph, "diffnet")) {
      stop("-cumadopt- should be provided when -graph- is not of class 'diffnet'")
    } else {
      cumadopt <- graph$cumadopt
    }

  # Checking diffnet graph
  if (inherits(graph, "diffnet")) graph <- graph$graph

  # Checking attrs
  if (!length(attrs)) {
    attrs <- matrix(1, ncol=ncol(cumadopt), nrow=nrow(cumadopt))
  }

  # Checking alt graph
  if (length(alt.graph)) {
    graph <- if (inherits(alt.graph, "character")) {
      if (alt.graph != "se") stop("Only character -alt.graph- value allowed is \"se\".")

      se <- lapply(struct_equiv(graph, groupvar=groupvar, ...), "[[", "SE")
      se <- lapply(se, function(x) {
        x <- 1/x
        x[!is.finite(x)] <- 0
        x
        })

      # Changing valued
      if (!valued) {
        warning("To use alt.graph=\"se\" -valued- has been switched to TRUE.")
        valued <- TRUE
      }

      se

    } else alt.graph


  }

  switch (class(graph),
    array   = exposure.array(graph, cumadopt, attrs, outgoing, valued, normalized),
    list    = exposure.list(graph, cumadopt, attrs, outgoing, valued, normalized),
    diffnet = exposure.list(graph, cumadopt, attrs, outgoing, valued, normalized),
    stopifnot_graph(graph)
  )
}

# @rdname exposure
# @export
exposure.array <- function(
  graph, cumadopt, attrs,
  outgoing, valued, normalized) {

  # Preparing the data
  n <- nrow(graph)
  t <- dim(graph)[3]
  graphl <- vector("list", t)
  for (i in 1:t)
    graphl[[i]] <- methods::as(graph[,,i], "dgCMatrix")

  # attrs can be either
  #  degree, indegree, outdegree, or a user defined vector.
  #  by default is user equal to 1
  da <- dim(attrs)
  if (!length(da)) stop("-attrs- must be a matrix of size n by T.")
  if (any(da != dim(cumadopt))) stop("Incorrect size for -attrs-. ",
                                     "It must be of size that -cumadopt-.")

  # Dimnames
  rn <- rownames(cumadopt)
  if (!length(rn)) rn <- 1:nrow(cumadopt)

  tn <- colnames(cumadopt)
  if (!length(tn)) tn <- 1:ncol(cumadopt)

  # Calculating the exposure, and asigning names
  output <- exposure_for(graphl, cumadopt, attrs, outgoing, valued, normalized)
  dimnames(output) <- list(rn, tn)
  output
}

# @rdname exposure
# @export
exposure.list <- function(
  graph, cumadopt, attrs,
  outgoing, valued, normalized) {

  # attrs can be either
  #  degree, indegree, outdegree, or a user defined vector.
  #  by default is user equal to 1
  da <- dim(attrs)
  if (!length(da)) stop("-attrs- must be a matrix of size n by T.")
  if (any(da != dim(cumadopt))) stop("Incorrect size for -attrs-. ",
                                     "It must be of size that -cumadopt-.")

  n <- nrow(graph[[1]])
  t <- length(graph)

  # Coercing into dgCMatrices
  test <- !sapply(graph, inherits, what="dgCMatrix")
  if (any(test))
    graph[which(test)] <- lapply(graph[which(test)],
                                 function(x) methods::as(x, "dgCMatrix"))

  output <- exposure_for(graph, cumadopt, attrs, outgoing, valued, normalized)

  rn <- rownames(cumadopt)
  if (!length(rn)) rn <- 1:nrow(cumadopt)

  tn <- colnames(cumadopt)
  if (!length(tn)) tn <- 1:ncol(cumadopt)

  dimnames(output) <- list(rn, tn)
  output
}

exposure_for <- function(graph, cumadopt, attrs, outgoing, valued, normalized) {
  out <- matrix(nrow = nrow(cumadopt), ncol = ncol(cumadopt))
  for (i in 1:nslices(graph))
    out[,i]<-exposure_cpp(graph[[i]], cumadopt[,i,drop=FALSE], attrs[,i,drop=FALSE],
                 outgoing, valued, normalized)
  return(out)
}

#' Cummulative count of adopters
#'
#' For each time period, calculates the number of adopters, the proportion of adopters,
#' and the adoption rate.
#'
#' @param obj A \eqn{n\times T}{n * T} matrix (Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}) or a \code{\link{diffnet}} object.
#' @details
#'
#' The rate of adoption--returned in the 3rd row out the resulting
#' matrix--is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{q_{t-1}}}{[q(t) - q(t-1)]/q(t-1)}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}. Note that
#' it is only calculated fot \eqn{t>1}.
#' @return A \eqn{3\times T}{3 * T} matrix, where its rows contain the number of adoptes, the proportion of
#' adopters and the rate of adoption respectively, for earch period of time.
#' @family statistics
#' @keywords univar
#' @export
#' @author George G. Vega Yon, Stephanie R. Dyal, Timothy B. Hayes & Thomas W. Valente
cumulative_adopt_count <- function(obj) {

  if (inherits(obj, "diffnet")) obj <- obj$cumadopt

  x <- cumulative_adopt_count_cpp(obj)
  dimnames(x) <- list(c("num", "prop", "rate"), colnames(obj))
  return(x)
}

#' Network Hazard Rate
#'
#' The hazard rate is the instantaneous probability of adoption at each time
#' representing the likelihood members will adopt at that time (Allison 1984).
#' The shape of the hazard rate indicates the pattern of new adopters over time.
#' Rapid diffusion with convex cumulative adoption curves will have hazard functions
#' that peak early and decay over time whereas slow concave cumulative adoption
#' curves will have hazard functions that are low early and rise over time.
#' Smooth hazard curves indicate constant adoption whereas those that oscillate
#' indicate variability in adoption behavior over time.
#' @aliases plot_hazarrate
#' @param obj A \eqn{n\times T}{n * T} matrix (Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}) or a \code{\link{diffnet}} object.
#' @param x An object of class \code{diffnet_hr}.
#' @param y ignored.
#' @param main Character scalar. Title of the plot
#' @param xlab Character scalar. x-axis label.
#' @param ylab Character scalar. y-axis label.
#' @param include.grid Logical scalar. When TRUE includes a grid on the plot.
#' @param bg Character scalar. Color of the points.
#' @param pch Integer scalar. See \code{\link{par}}.
#' @param type Character scalar. See \code{\link{par}}.
#' @param no.plot Logical scalar. When TRUE, suppress plotting (only returns hazard rates).
#' @param add Logical scalar. When TRUE it adds the hazard rate to the current plot.
#' @param ylim Numeric vector. See \code{\link{plot}}.
#' @param ... further arguments to be passed to \code{\link{par}}
#' @details
#'
#' This function computes hazard rate, plots it and returns the hazard rate vector
#' invisible (so is not printed on the console). For \eqn{t>1}, hazard rate is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{n - q_{t-1}}}{[q(t) - q(t-1)]/[n - q(t-1)]}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}, and \eqn{n} is
#' the number of vertices in the graph.
#'
#' In survival analysis, hazard rate is defined formally as
#'
#' \deqn{%
#' \lambda(t)=\lim_{h\to +0}\frac{F(t+h)-F(t)}{h}\frac{1}{1-F(t)} %
#' }{%
#' \lambda(t-1)= lim (t -> +0) [F(t+h)-F(t)]/h * 1/[1-F(t)] %
#' }
#'
#' Then, by approximating \eqn{h=1}, we can rewrite the equation as
#'
#' \deqn{%
#' \lambda(t)=\frac{F(t+1)-F(t)}{1-F(t)} %
#' }{%
#' \lambda(t-1)= [F(t+1)-F(t)]/[1-F(t)] %
#' }
#'
#' Furthermore, we can estimate \eqn{F(t)}, the probability of not having adopted
#' the innovation in time \eqn{t}, as the proportion of adopters in that time, this
#' is \eqn{F(t) \sim q_t/n}{F(t) ~ q(t)/n}, so now we have
#'
#' \deqn{%
#' \lambda(t)=\frac{q_{t+1}/n-q_t/n}{1-q_t/n} = \frac{q_{t+1} - q_t}{n - q_t} %
#' }{%
#' \lambda(t-1)= [q(t+1)/n-q(t)/n]/[1-q(t)/n] = [q(t+1) - q(t)]/[n - q(t)] %
#' }
#'
#' As showed above.
#'
#' The \code{plot_hazard} function is an alias for the \code{plot.diffnet_hr} method.
#' @return A row vector of size \eqn{T} with hazard rates for \eqn{t>1} of class \code{diffnet_hr}.
#' The class of the object is only used by the S3 plot method.
#' @family statistics
#' @family visualizations
#' @keywords univar
#' @examples
#' # Creating a random vector of times of adoption
#' toa <- sample(2000:2005, 20, TRUE)
#'
#' # Computing cumulative adoption matrix
#' cumadopt <- toa_mat(toa)$cumadopt
#'
#' # Visualizing the hazard rate
#' hazard_rate(cumadopt)
#' @references
#' Allison, P. (1984). Event history analysis regression for longitudinal event
#' data. Beverly Hills: Sage Publications.
#'
#' Wooldridge, J. M. (2010). Econometric Analysis of Cross Section and Panel Data
#' (2nd ed.). Cambridge: MIT Press.
#' @export
#' @author George G. Vega Yon, Stephanie R. Dyal, Timothy B. Hayes & Thomas W. Valente
hazard_rate <- function(obj, no.plot=FALSE, include.grid=TRUE, ...) {

  if (inherits(obj, "diffnet")) obj <- obj$cumadopt

  x <- hazard_rate_cpp(obj)
  dimnames(x) <- list("hazard", colnames(obj))
  class(x) <- c("diffnet_hr", class(x))
  if (!no.plot) plot.diffnet_hr(x, include.grid=include.grid, ...)
  invisible(x)
}

#' @rdname hazard_rate
#' @export
plot_hazard <- function(x,main="Hazard Rate", xlab="Time", ylab="Hazard Rate", type="b",
                        include.grid=TRUE, bg="lightblue", add=FALSE, ylim=c(0,1), pch=21,
                        ...) {
  hr <- hazard_rate(x, no.plot = TRUE)
  plot.diffnet_hr(x=hr, main=main, xlab=xlab, ylab=ylab, type=type, include.grid=include.grid, bg=bg,
                  add=add, ylim=ylim, pch=pch, ...)
}

#' @rdname hazard_rate
#' @export
plot.diffnet_hr <- function(x,y=NULL, main="Hazard Rate", xlab="Time",
                            ylab="Hazard Rate", type="b",
                            include.grid=TRUE, bg="lightblue", pch=21, add=FALSE, ylim=c(0,1),
                            ...) {

  if (add) {
    lines(y=t(x), x=colnames(x), type=type, bg=bg, pch=pch, ...)
  } else {
    plot(y=t(x), x=colnames(x), type=type, main=main, xlab=xlab, ylab=ylab,
         ylim=ylim, bg=bg, pch=pch,...)
    if (include.grid) grid()
  }

  invisible(x)
}

#' Retrive threshold levels from the exposure matrix
#'
#' Thresholds are each vertexes exposure at the time of adoption.
#' Substantively it is the proportion of adopters required for each ego to adopt. (see \code{\link{exposure}}).
#'
#' @param obj Either a \eqn{n\times T}{n * T} matrix (eposure to the innovation obtained from
#' \code{\link{exposure}}) or a \code{diffnet} object.
#' @param toa Integer vector. Indicating the time of adoption of the innovation.
#' @param t0 Integer scalar. See \code{\link{toa_mat}}.
#' @param include_censored Logical scalar. When \code{TRUE} (default), threshold
#' levels are not reported for observations adopting in the first time period.
#' @param ... Further arguments to be passed to \code{\link{exposure}}.
#' @return A vector of size \eqn{n} indicating the threshold for each node.
#' @family statistics
#' @seealso Threshold can be visualized using \code{\link{plot_threshold}}
#' @keywords univar
#' @details By default exposure is not computed for vertices adopting at the
#' first time period, \code{include_censored=FALSE}, as estimating threshold for
#' left censored data may yield biased outcomes.
#' @examples
#' # Generating a random graph with random Times of Adoption
#' set.seed(783)
#' toa <- sample.int(4, 5, TRUE)
#' graph <- rgraph_er(n=5, t=max(toa) - min(toa) + 1)
#'
#' # Computing exposure using Structural Equivalnece
#' adopt <- toa_mat(toa)
#' se <- struct_equiv(graph)
#' se <- lapply(se, function(x) methods::as((x$SE)^(-1), "dgCMatrix"))
#' expo <- exposure(graph, adopt$cumadopt, alt.graph=se)
#'
#' # Retrieving threshold
#' threshold(expo, toa)
#'
#' # We can do the same by creating a diffnet object
#' diffnet <- as_diffnet(graph, toa)
#' threshold(diffnet, alt.graph=se)
#' @export
#' @author George G. Vega Yon, Stephanie R. Dyal, Timothy B. Hayes & Thomas W. Valente
threshold <- function(obj, toa, t0=min(toa, na.rm = TRUE), include_censored=FALSE, ...) {

  if (inherits(obj, "diffnet")) {
    t0 <- min(obj$meta$pers)
    toa <- obj$toa
    obj <- exposure(obj, ...)
  } else {
    if (missing(toa))
      stop("-toa- should be provided when -obj- is not of class 'diffnet'")
  }

  toa <- toa - t0 + 1L
  output <- threshold_cpp(obj, toa, include_censored)
  dimnames(output) <- list(rownames(obj), "threshold")

  # Correcting weird cases
  if (!include_censored) output[toa==1] <- NA
  output[is.na(toa)] <- NA
  output
}
