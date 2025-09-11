#' Indegree, outdegree and degree of the vertices
#'
#' Computes the requested degree measure for each node in the graph.
#' @templateVar undirected TRUE
#' @templateVar self TRUE
#' @templateVar valued TRUE
#' @template graph_template
#' @param cmode Character scalar. Either "indegree", "outdegree" or "degree".
#' @return A numeric matrix of size \eqn{n\times T}{n * T}. In the case of \code{plot},
#'  returns an object of class \code{\link[graphics:hist]{histogram}}.
#' @export
#' @family statistics
#' @family visualizations
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
#' # Classic Scale-free plot ---------------------------------------------------
#' set.seed(1122)
#' g <- rgraph_ba(t=1e3-1)
#' hist(dgr(g))
#'
#' # Since by default uses logscale, here we suppress the warnings
#' # on points been discarded for <=0.
#' suppressWarnings(plot(dgr(g)))
#'
#' @author George G. Vega Yon
dgr <- function(graph, cmode="degree",
                undirected=getOption("diffnet.undirected", FALSE),
                self=getOption("diffnet.self",FALSE),
                valued=getOption("diffnet.valued", FALSE)) {

  cls <- class(graph)
  ans <- if ("matrix" %in% cls) {
    dgr.matrix(graph, cmode, undirected, self, valued)
  } else if ("array" %in% cls) {
    dgr.array(graph, cmode, undirected, self, valued)
  } else if ("dgCMatrix" %in% cls) {
    dgr.dgCMatrix(graph, cmode, undirected, self, valued)
  } else if ("list" %in% cls) {
    dgr.list(graph, cmode, undirected, self, valued)
  } else if ("diffnet" %in% cls) {
    dgr.list(graph$graph, cmode, undirected = graph$meta$undirected, self, valued)
  } else if ("igraph" %in% cls) {
    graph <- as_generic_graph.igraph(graph)
    dgr.dgCMatrix(graph$graph[[1]], cmode, graph$meta$undirected, self, valued)
  } else if ("network" %in% cls) {
    graph <- as_generic_graph.network(graph)
    dgr.dgCMatrix(graph$graph[[1]], cmode, graph$meta$undirected, self, valued)
  } else if ("matrix.csc" %in% cls) {
    dgr.dgCMatrix(methods::as(graph, "dgCMatrix"))
  } else
    stopifnot_graph(graph)

  return(structure(ans, class=c("diffnet_degSeq", class(ans))))
}

# cmode:
#  0: in
#  1: out
#  2: total
.dgr <- function(graph, cmode, undirected, self, valued) {

  if (cmode < 0 || cmode > 2) stop("Invalid degree.")

  # Checking if it is valued or not
  n <- nrow(graph)
  if (ncol(graph) != n)
    stop("-graph- should be a square matrix.")

  if (!valued)
    graph@x <- rep(1, length(graph@x))

  # Computing sums
  ans <- if (cmode == 2) {
    if (undirected)  Matrix::rowSums(graph)
    else (Matrix::rowSums(graph) + Matrix::colSums(graph))
  } else if (cmode == 0) {
    Matrix::colSums(graph)
  } else {
    Matrix::rowSums(graph)
  }

  # Checking if self or not
  ans <- if (!self) {
    if (cmode == 2) {
      if (!undirected) ans - Matrix::diag(graph)*2
      else ans - Matrix::diag(graph)
    } else {
      ans - Matrix::diag(graph)
    }
  } else ans

  return(matrix(ans, ncol=1))
}

#' @export
#' @rdname dgr
#' @param x An \code{diffnet_degSeq object}
#' @param breaks Passed to \code{\link{hist}}.
#' @param log Passed to \code{\link{plot}} (see \code{\link{par}}).
#' @param hist.args Arguments passed to \code{\link{hist}}.
#' @param xlab Character scalar. Passed to \code{\link{plot}}.
#' @param ylab Character scalar. Passed to \code{\link{plot}}.
#' @param ... Further arguments passed to \code{\link{plot}}.
#' @param slice Integer scalar. In the case of dynamic graphs, number of time
#'  point to plot.
#' @param y Ignored
#' @param freq Logical scalar. When \code{TRUE} the y-axis will reflex counts,
#'  otherwise densities.
plot.diffnet_degSeq <- function(
  x,
  breaks = min(100L, nrow(x)/5),
  freq=FALSE,
  y=NULL,
  log="xy",
  hist.args=list(),
  slice=ncol(x),
  xlab="Degree",
  ylab="Freq",
  ...
  ) {

  ans <- do.call(hist, c(hist.args, list(x=x[,slice], breaks = breaks, plot=FALSE)))
  with(ans, plot(x=mids,y=if (freq) counts else density,log=log, xlab=xlab, ylab=ylab,...))
  invisible(ans)

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
  output <- .dgr(methods::as(graph, "dgCMatrix"), cmode, undirected, self,
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
  output <- .dgr(graph, cmode, undirected, self, valued)

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
#' based on (1) incoming ties, (2) structural equivalence, (3) indirect ties,
#' (4) network-metric weighted (e.g., central nodes have more
#' influence), and (5) attribute-weighted (e.g., based on homophily or tie strength).
#'
#' @templateVar valued TRUE
#' @templateVar dynamic TRUE
#' @templateVar self TRUE
#' @template graph_template
#' @param cumadopt \eqn{n\times T}{n * T} matrix for single diffusion.
#' \eqn{n\times T \times Q}{n * T * Q} array for \eqn{Q}{Q} diffusion processes.
#' Cumulative adoption matrix obtained from \code{\link{toa_mat}}
#' @param attrs Either a character scalar (if \code{graph} is diffnet),
#' a numeric matrix of size \eqn{n\times T}{n * T}, or an array of size
#' \eqn{n\times T \times Q}{n * T * Q} (only for multi diffusion).
#'  Weighting for each time period (see details).
#' @param alt.graph Either a graph that should be used instead of \code{graph},
#' or \code{"se"} (see details).
#' @param outgoing Logical scalar. When \code{TRUE}, computed using outgoing ties.
#' @param normalized Logical scalar. When \code{TRUE}, the exposure will be between zero
#' and one (see details).
#' @param ... Further arguments passed to \code{\link{struct_equiv}} (only used when
#' \code{alt.graph="se"}).
#' @param groupvar Passed to \code{\link{struct_equiv}}.
#' @param lags Integer scalar. When different from 0, the resulting exposure
#' matrix will be the lagged exposure as specified (see examples).
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
#' If the \code{alt.graph} is static, then the function will warn about it
#' and will recycle the graph to compute exposure at each time point.
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
#' case of a social network exposure will be computed based on the incoming ties.
#'
#' If \code{normalize=FALSE} then denominator, \eqn{S_t \times x_t}{S(t) \%*\% x(t)},
#' is not included. This can be useful when, for example, exposure needs to be
#' computed as a count instead of a proportion. A good example of this can be
#' found at the examples section of the function \code{\link{rdiffnet}}.
#'
#' @references
#' Burt, R. S. (1987). "Social Contagion and Innovation: Cohesion versus Structural
#' Equivalence". American Journal of Sociology, 92(6), 1287.
#' \doi{10.1086/228667}
#'
#' Valente, T. W. (1995). "Network models of the diffusion of innovations"
#'  (2nd ed.). Cresskill N.J.: Hampton Press.
#'
#' @examples
#' # Calculating lagged exposure -----------------------------------------------
#'
#' set.seed(8)
#' graph <- rdiffnet(20, 4)
#'
#' expo0 <- exposure(graph)
#' expo1 <- exposure(graph, lags = 1)
#'
#' # These should be equivalent
#' stopifnot(all(expo0[, -4] == expo1[, -1])) # No stop!
#'
#'
#' # Calculating the exposure based on Structural Equivalence ------------------
#' set.seed(113132)
#' graph <- rdiffnet(100, 4)
#'
#' SE <- lapply(struct_equiv(graph), "[[", "SE")
#' SE <- lapply(SE, function(x) {
#'    x <- 1/x
#'    x[!is.finite(x)] <- 0
#'    x
#' })
#'
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
#' # Comparing using incoming edges -------------------------------------------
#' eIN <- exposure(graph, outgoing=FALSE)
#'
#' # Structral equivalence for different communities ---------------------------
#' data(medInnovationsDiffNet)
#'
#' # Only using 4 time slides, this is for convenience
#' medInnovationsDiffNet <- medInnovationsDiffNet[, , 1:4]
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
#' stopifnot(all(test[!is.na(test)]))
#'
#' # METHOD 2: Using the 'groupvar' argument
#' # Further, we can compare this with using the groupvar
#' diffnet[["expo_se2"]] <- exposure(diffnet, alt.graph="se",
#'    groupvar="city", valued=TRUE)
#'
#' # These should be equivalent
#' test <- diffnet[["expo_se", as.df=TRUE]] == diffnet[["expo_se2", as.df=TRUE]]
#' stopifnot(all(test[!is.na(test)]))
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
#' stopifnot(all(test[!is.na(test)]))
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
#' stopifnot(all(test[!is.na(test)]))
#'
#'
#' # Examples for multi-diffusion ---------------------------
#'
#' # Running a multi-diffusion simulation, with q=2 behaviors
#' set.seed(999)
#' n <- 40; t <- 5; q <- 2;
#' graph <- rgraph_ws(n, t, p=.3)
#' seed_prop_adopt <- rep(list(0.1), q)
#'
#' diffnet <- rdiffnet(seed.graph = graph, t = t, seed.p.adopt = seed_prop_adopt)
#'
#' # Getting the cumulative adoption array of dims n x T x q
#' cumadopt_2 <- diffnet$cumadopt  # list of matrices
#' cumadopt_2 <- array(unlist(cumadopt_2), dim = c(n, t, q))
#'
#' expo2 <- exposure(diffnet$graph, cumadopt = cumadopt_2)
#'
#' # With an attribute --
#'
#' X <- matrix(runif(n * t), nrow = n, ncol = t) # matrix n x T
#' ans3 <- exposure(diffnet$graph, cumadopt = cumadopt_2, attrs=X)
#'
#' X <- array(runif(n * t * q), dim = c(n, t, q)) # array n x T x q
#' ans4 <- exposure(diffnet$graph, cumadopt = cumadopt_2, attrs=X)
#'
#' # Exposure based on Structural Equivalence --
#'
#' diffnet_1 <- split_behaviors(diffnet)[[1]]
#' se <- struct_equiv(diffnet)
#' se <- lapply(se, function(x) {
#'   ans <- methods::as(x$SE, "dgCMatrix")
#'     ans@x <- 1/(ans@x + 1e-20)
#'       ans
#'       })
#' ans6 <- exposure(diffnet, cumadopt = cumadopt_2, alt.graph = se, valued=TRUE)
#'
#' @family statistics
#' @keywords univar
#' @return A matrix of size \eqn{n\times T}{n * T} with exposure for each node.
#' @export
#' @author George G. Vega Yon, Thomas W. Valente, and Aníbal Olivera M.
#' @name exposure
NULL

# Workhorse of exposure plotting
.exposure <- function(graph, cumadopt, attrs, outgoing, valued, normalized, self) {

  # Getting the parameters
  n <- nrow(graph)
  if (n!=ncol(graph))
    stop("-graph- is not squared.")

  # Checking values
  if (!valued)
    graph@x <- rep(1, length(graph@x))

  # Direction of the exposure
  if (!outgoing)
    graph <- t(graph)

  # Checking self
  if (!self) graph <- sp_diag(graph, rep(0, nnodes(graph)))

  norm <- graph %*% attrs + 1e-20

  if (!is.na(dim(cumadopt)[3])) {
    ans <- array(0, dim = c(dim(cumadopt)[1],dim(cumadopt)[3]))

    for (q in 1:dim(cumadopt)[3]) {
      if (normalized) {
        ans[,q] <- as.vector(graph %*% (attrs * cumadopt[,,q]) / norm)
      } else {
        ans[,q] <- as.vector(graph %*% (attrs * cumadopt[,,q]))
      }
    }
  } else {
    ans <- graph %*% (attrs * cumadopt)

    if (normalized) {
      ans <- ans/ norm
    }
  }

  return(as.vector(ans))
}

check_lags <- function(npers, lags) {

  # Checking length
  if (length(lags) != 1L)
    stop("-lags- should be a scalar (length 1). Right now it has lenght ",
         length(lags))

  # Checking class
  lags <- as.integer(lags)
  if (is.na(lags))
    stop("-lags- cannot be NA. It should be an integer scalar.")

  # Should fit the range of data. A lag cannot be greater than npers.
  # it has to be strictly smaller
  if (abs(lags) >= npers)
    stop("-abs(lags)- cannot be greater than ",npers-1L,". Right now lags=",lags,".")

  lags

}

#' @export
#' @rdname exposure
exposure <- function(
  graph,
  cumadopt,
  attrs      = NULL,
  alt.graph  = NULL,
  outgoing   = getOption("diffnet.outgoing", TRUE),
  valued     = getOption("diffnet.valued", FALSE),
  normalized = TRUE,
  groupvar   = NULL,
  self       = getOption("diffnet.self"),
  lags       = 0L,
  ...
  ) {

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
      cumadopt <- toa_mat(graph)$cumadopt
    }

  # Checking diffnet graph
  if (inherits(graph, "diffnet")) graph <- graph$graph

  # Checking attrs
  if (!length(attrs)) {
    if (!is.na(dim(cumadopt)[3])) {
    attrs <- array(1, dim = c(nrow(cumadopt), ncol(cumadopt), 1))}
    else {attrs <- matrix(1, ncol=ncol(cumadopt), nrow=nrow(cumadopt))}
  } else if (!is.na(dim(cumadopt)[3])) {
    attrs <- array(attrs, dim = c(nrow(attrs), ncol(attrs), 1))
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

    } else {
      # In the case of static nets
      if (inherits(alt.graph, "matrix"))
        alt.graph <- methods::as(alt.graph, "dgCMatrix")

      if (inherits(alt.graph, "dgCMatrix")) {
        warning("When -alt.graph- is static, will be repeated \"t\" times to fit the data.")
        alt.graph <- replicate(nslices(graph), alt.graph)
      }

      if (!valued)
        warning("The -alt.graph- will be treated as 0/1 graph (value=FALSE).")
      alt.graph
    }


  }

  # Checking lags
  lags <- check_lags(nslices(graph), lags)

  if ((is.array(graph) & !inherits(graph, "matrix")) | is.list(graph)) {
    exposure.list(as_spmat(graph), cumadopt, attrs, outgoing, valued, normalized,
                  self, lags)
  } else stopifnot_graph(graph)
}

# @rdname exposure
# @export
exposure.list <- function(
  graph, cumadopt, attrs,
  outgoing, valued, normalized, self, lags) {

  # attrs can be either
  #  degree, indegree, outdegree, or a user defined vector.
  #  by default is user equal to 1

  # dim(attrs) default n x T matrix of 1's
  if (!length(dim(attrs))) stop("-attrs- must be a matrix of size n by T.")

  if (!is.na(dim(cumadopt)[3])) {
    if (dim(cumadopt)[3]>1 && any(dim(attrs)[-3] != dim(cumadopt)[-3])) stop("Incorrect size for -attrs-. ",
                                              "Does not match n dim or t dim.")
  } else {
    if (any(dim(attrs) != dim(cumadopt))) stop("Incorrect size for -attrs-. ",
                                              "It must be of size that -cumadopt-.")
  }
  add_dimnames.mat(cumadopt)
  add_dimnames.mat(attrs)

  output <- exposure_for(graph, cumadopt, attrs, outgoing, valued, normalized,
                         self, lags)

  dimnames(output) <- dimnames(cumadopt)
  output

}

exposure_for <- function(
  graph,
  cumadopt,
  attrs,
  outgoing,
  valued,
  normalized,
  self,
  lags
  ) {

  if (!is.na(dim(cumadopt)[3])) {
    out <- array(NA, dim = c(dim(cumadopt)[1], dim(cumadopt)[2], dim(cumadopt)[3]))

    if (lags >= 0L) {
      for (i in 1:(nslices(graph) - lags)) {
        out[, i + lags, ] <- .exposure(graph[[i]],
                                       cumadopt[, i, , drop = FALSE],
                                       attrs[, i, , drop = FALSE],
                                       outgoing = outgoing,
                                       valued = valued,
                                       normalized = normalized,
                                       self = self)
      }
    } else {
      for (i in (1 - lags):nslices(graph)) {
        out[, i + lags, ] <- .exposure(graph[[i]],
                                       cumadopt[, i, , drop = FALSE],
                                       attrs[, i, , drop = FALSE],
                                       outgoing = outgoing,
                                       valued = valued,
                                       normalized = normalized,
                                       self = self)
      }
    }
  } else {
    out <- array(NA, dim = c(dim(cumadopt)[1], dim(cumadopt)[2]))

    if (lags >= 0L) {
      for (i in 1:(nslices(graph) - lags)) {
        out[, i + lags] <- .exposure(graph[[i]],
                                     cumadopt[, i, drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = outgoing,
                                     valued = valued,
                                     normalized = normalized,
                                     self = self)
      }
    } else {
      for (i in (1 - lags):nslices(graph)) {
        out[, i + lags] <- .exposure(graph[[i]],
                                     cumadopt[, i, drop = FALSE],
                                     attrs[, i, drop = FALSE],
                                     outgoing = outgoing,
                                     valued = valued,
                                     normalized = normalized,
                                     self = self)
      }
    }
  }

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
#' @author George G. Vega Yon & Thomas W. Valente
cumulative_adopt_count <- function(obj) {

  if (inherits(obj, "diffnet")) x <- obj$cumadopt
  else x <- obj

  # Checking colnames
  cn <- if (inherits(obj, "diffnet")) obj$meta$pers
  else colnames(obj)
  if (length(cn) == 0) cn <- as.character(1:ncol(obj))

  q <- colSums(x)
  t <- length(q)
  structure(
    rbind(
      q,
      q/nrow(x),
      c(0,(q[-1] - q[-t])/(q[-t] + 1e-15))
    ), dimnames = list(c("num", "prop", "rate"), cn)
  )
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
#' @param ... further arguments to be passed to the method.
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
#' @author George G. Vega Yon & Thomas W. Valente
hazard_rate <- function(obj, no.plot=FALSE, include.grid=TRUE, ...) {
  if (inherits(obj, "diffnet")) {
    dn  <- with(obj$meta, list(ids, pers))
    obj <- obj$cumadopt
    dimnames(obj) <- dn
  } else {
    if (!length(colnames(obj)))
      colnames(obj) <- seq_len(ncol(obj))
  }

  q <- colSums(obj)
  t <- length(q)

  x <- structure(
    rbind(c(0,(q[-1] - q[-t])/(nrow(obj) - q[-t] + 1e-15)))
    , dimnames = list("hazard", colnames(obj)),
    class=c("diffnet_hr", "matrix")
  )

  if (!no.plot) plot.diffnet_hr(x, include.grid=include.grid, ...)
  invisible(x)

}

#' @rdname hazard_rate
#' @export
plot_hazard <- function(x, ...) {
  hr <- hazard_rate(x, no.plot = TRUE)

  dots <- list(...)
  do.call(plot.diffnet_hr, c(list(x=hr), dots))
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
#' @param lags Integer scalar. Number of lags to consider when computing thresholds. \code{lags=1}
#'  defines threshold as exposure at \eqn{T-1}, where \code{T} is time of adoption.
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
#' @author George G. Vega Yon & Thomas W. Valente
threshold <- function(obj, toa, t0=min(toa, na.rm = TRUE), include_censored=FALSE,
                       lags=0L, ...) {

  if (inherits(obj, "diffnet")) {
    t0 <- min(obj$meta$pers)
    toa <- obj$toa
    obj <- exposure(obj, ...)
  } else {
    if (missing(toa))
      stop("-toa- should be provided when -obj- is not of class 'diffnet'")
  }

  toa <- toa - t0 + 1L

  # If lags are included
  toa <- toa - lags
  toa[(toa < 1) | (toa > ncol(obj))] <- NA

  ans <- structure(obj[cbind(1:length(toa),toa)], dim=c(length(toa),1),
            dimnames=list(rownames(obj), "threshold"))

  # Checking if whether to included censored or not
  if (include_censored) return(ans)
  else {
    ans[which(is.na(toa))] <- NA
    return(ans)
  }
}

#' Classify adopters accordingly to Time of Adoption and Threshold levels.
#'
#' Adopters are classified as in Valente (1995). In general, this is done
#' depending on the distance in terms of standard deviations from the mean of
#' Time of Adoption and Threshold.
#'
#' @param graph A dynamic graph.
#' @param toa Integer vector of length \eqn{n} with times of adoption.
#' @param t0 Integer scalar passed to \code{\link{threshold}} and \code{\link{toa_mat}}.
#' @param t1 Integer scalar passed to \code{\link{toa_mat}}.
#' @param expo Numeric matrix of size \eqn{n\times T}{n*T} with network exposures.
#' @param include_censored Logical scalar, passed to \code{\link{threshold}}.
#' @param x A \code{diffnet_adopters} class object.
#' @param ... Further arguments passed to the method.
#' @export
#' @details
#' Classifies (only) adopters according to time of adoption and threshold as
#' described in Valente (1995). In particular, the categories are defined as follow:
#'
#' For Time of Adoption, with \code{toa} as the vector of times of adoption:
#' \itemize{
#'  \item \emph{Early Adopters}: \code{toa[i] <= mean(toa) - sd(toa)},
#'  \item \emph{Early Majority}: \code{mean(toa) - sd(toa) < toa[i] <= mean(toa) },
#'  \item \emph{Late Majority}: \code{mean(toa) < toa[i] <= mean(toa) + sd(toa) }, and
#'  \item \emph{Laggards}: \code{mean(toa) + sd(toa) < toa[i] }.
#' }
#'
#' For Threshold levels, with \code{thr} as the vector of threshold levels:
#' \itemize{
#'  \item \emph{Very Low Thresh.}: \code{thr[i] <= mean(thr) - sd(thr)},
#'  \item \emph{Low Thresh.}: \code{mean(thr) - sd(thr) < thr[i] <= mean(thr) },
#'  \item \emph{High Thresh.}: \code{mean(thr) < thr[i] <= mean(thr) + sd(thr) }, and
#'  \item \emph{Very High. Thresh.}: \code{mean(thr) + sd(thr) < thr[i] }.
#' }
#'
#' By default threshold levels are not computed for left censored data. These
#' will have a \code{NA} value in the \code{thr} vector.
#'
#' The plot method, \code{plot.diffnet_adopters}, is a wrapper for the
#' \code{\link[graphics:plot.table]{plot.table}} method. This generates a
#' \code{\link[graphics:mosaicplot]{mosaicplot}} plot.
#'
#' @return A list of class \code{diffnet_adopters} with the following elements:
#' \item{toa}{A factor vector of length \eqn{n} with 4 levels:
#'  "Early Adopters", "Early Majority", "Late Majority", and "Laggards"}
#' \item{thr}{A factor vector of length \eqn{n} with 4 levels:
#'  "Very Low Thresh.", "Low Thresh.", "High Thresh.", and "Very High Thresh."}
#' @examples
#' # Classifying brfarmers -----------------------------------------------------
#'
#' x <- brfarmersDiffNet
#' diffnet.toa(x)[x$toa==max(x$toa, na.rm = TRUE)] <- NA
#' out <- classify_adopters(x)
#'
#' # This is one way
#' round(
#' with(out, ftable(toa, thr, dnn=c("Time of Adoption", "Threshold")))/
#'   nnodes(x[!is.na(x$toa)])*100, digits=2)
#'
#' # This is other
#' ftable(out)
#'
#' # Can be coerced into a data.frame, e.g. ------------------------------------
#'  str(classify(brfarmersDiffNet))
#'  ans <- cbind(
#'  as.data.frame(classify(brfarmersDiffNet)), brfarmersDiffNet$toa
#'  )
#'  head(ans)
#'
#' # Creating a mosaic plot with the medical innovations -----------------------
#' x <- classify(medInnovationsDiffNet)
#' plot(x)
#'
#' @family statistics
#' @references
#' Valente, T. W. (1995). "Network models of the diffusion of innovations"
#'  (2nd ed.). Cresskill N.J.: Hampton Press.
#' @author George G. Vega Yon
classify_adopters <- function(...) UseMethod("classify_adopters")

#' @export
#' @rdname classify_adopters
classify <- classify_adopters

#' @export
#' @rdname classify_adopters
classify_adopters.diffnet <- function(graph, include_censored=FALSE, ...) {
  classify_adopters.default(graph$graph, graph$toa,
                            t0=graph$meta$pers[1], t1=NULL,
                            expo=exposure(graph, ...),
                            include_censored=include_censored)
}

#' @export
#' @rdname classify_adopters
classify_adopters.default <- function(
  graph,
  toa,
  t0=NULL,
  t1=NULL,
  expo=NULL,
  include_censored=FALSE,
  ...
  ) {


  # Computing ranges
  ran_toa <- mean(toa, na.rm = TRUE) + sd(toa, na.rm = TRUE)*c(-1,0,1)

  # Getting threshold
  if (!length(expo)) expo <- exposure(graph, toa_mat(toa, t0=t0, t1=t1), ...)
  thr <- threshold(expo, toa, t0, include_censored)
  ran_thr <- mean(thr, na.rm = TRUE) + sd(thr, na.rm = TRUE)*c(-1,0,1)

  # Classifying
  c_toa <- c("Non-Adopters", "Early Adopters", "Early Majority", "Late Majority", "Laggards")
  c_toa <- factor(findInterval(toa, ran_toa), c(NA,0:3), c_toa, exclude=NULL)

  c_thr <- c("Non-Adopters", "Very Low Thresh.", "Low Thresh.", "High Thresh.", "Very High Thresh.")
  c_thr <- factor(findInterval(thr, ran_thr), c(NA,0:3), c_thr, exclude=NULL)

  structure(list(
    toa=c_toa,
    thr=c_thr,
    cutoffs=list(toa=ran_toa, thr=ran_thr)
    ), class="diffnet_adopters")
  #

}

#' @export
#' @param as.pcent Logical scalar. When \code{TRUE} returns a table with percentages
#' instead.
#' @param digits Integer scalar. Passed to \code{\link[base:round]{round}}.
#' @rdname classify_adopters
ftable.diffnet_adopters <- function(x, as.pcent=TRUE, digits=2, ...) {

  out <- with(x, stats::ftable(toa, thr, ...))

  if (as.pcent) round(out/sum(out)*100, digits)
  else out
}

#' @export
#' @param row.names Passed to \code{\link[base:as.data.frame]{as.data.frame}}.
#' @param optional Passed to \code{\link[base:as.data.frame]{as.data.frame}}.
#' @rdname classify_adopters
as.data.frame.diffnet_adopters <- function(x, row.names=NULL, optional=FALSE, ...) {
  as.data.frame(x[1:2], row.names, optional, ...)
}

#' @export
#' @param y Ignored.
#' @param ftable.args List of arguments passed to \code{\link[stats:ftable]{ftable}}.
#' @param table.args List of arguments passed to \code{\link{table}}.
#' @rdname classify_adopters
plot.diffnet_adopters <- function(x, y = NULL,
                                  ftable.args = list(),
                                  table.args=list(),...) {
  y <- do.call(ftable.diffnet_adopters, c(ftable.args, list(x=x)))
  y <- do.call(as.table, c(table.args, list(x=y)))
  plot(y, ...)
}

#' Computes covariate distance between connected vertices
#'
#' @param graph A square matrix of size \eqn{n} of class dgCMatrix.
#' @param X A numeric matrix of size \eqn{n \times K}{n * K}. Vertices attributes
#' @param p Numeric scalar. Norm to compute
#' @param S Square matrix of size \code{ncol(x)}. Usually the var-covar matrix.
#' @details
#'
#' Faster than \code{\link{dist}}, these functions compute distance metrics
#' between pairs of vertices that are connected (otherwise skip).
#'
#' The function \code{vertex_covariate_dist} is the simil of \code{\link{dist}}
#' and returns p-norms (Minkowski distance). It is implemented as follows (for
#' each pair of vertices):
#'
#' \deqn{%
#' D_{ij} = \left(\sum_{k=1}^K \left|X_{ik} - X_{jk}\right|^{p} \right)^{1/p}\mbox{ if }graph_{i,j}\neq 0
#' }{%
#' D(i,j) = [\sum_k abs(X(i,k) - X(j,k))^p]^(1/p)  if graph(i,j) != 0
#' }
#'
#' In the case of mahalanobis distance, for each pair of vertex \eqn{(i,j)}, the
#' distance is computed as follows:
#'
#' \deqn{%
#' D_{ij} = \left( (X_i - X_j)\times S \times (X_i - X_j)' \right)^{1/2}\mbox{ if }graph_{i,j}\neq 0
#' }{%
#' D(i,j) = sqrt[(X(i) - X(j)) \%*\% S  \%*\% t(X(i) - X(j))]  if graph(i,j) != 0
#' }
#'
#' @return A matrix of size \eqn{n\times n}{n * n} of class \code{dgCMatrix}. Will
#' be symmetric only if \code{graph} is symmetric.
#'
#' @export
#' @examples
#' # Distance (aka p norm) -----------------------------------------------------
#' set.seed(123)
#' G <- rgraph_ws(20, 4, .1)
#' X <- matrix(runif(40), ncol=2)
#'
#' vertex_covariate_dist(G, X)[1:5, 1:5]
#'
#' # Mahalanobis distance ------------------------------------------------------
#' S <- var(X)
#'
#' M <- vertex_mahalanobis_dist(G, X, S)
#'
#' # Example with diffnet objects ----------------------------------------------
#'
#' data(medInnovationsDiffNet)
#' X <- cbind(
#'   medInnovationsDiffNet[["proage"]],
#'   medInnovationsDiffNet[["attend"]]
#' )
#'
#' S <- var(X, na.rm=TRUE)
#' ans <- vertex_mahalanobis_dist(medInnovationsDiffNet, X, S)
#'
#' @name vertex_covariate_dist
#' @references
#' Mahalanobis distance. (2016, September 27). In Wikipedia, The Free Encyclopedia.
#' Retrieved 20:31, September 27, 2016, from
#' \url{https://en.wikipedia.org/w/index.php?title=Mahalanobis_distance&oldid=741488252}
#' @author George G. Vega Yon
#' @aliases p-norm mahalanobis minkowski
#' @family statistics
#' @family dyadic-level comparison functions
#' @seealso \code{\link[stats:mahalanobis]{mahalanobis}} in the stats package.
NULL

#' @export
#' @rdname vertex_covariate_dist
vertex_mahalanobis_dist <- function(graph, X, S) {

  # Analyzing
  cls <- class(graph)
  ans <- if ("matrix" %in% cls) {
    vertex_mahalanobis_dist_cpp(methods::as(graph, "dgCMatrix"), X, S)
  } else if ("dgCMatrix" %in% cls) {
    vertex_mahalanobis_dist_cpp(graph, X, S)
  } else if ("diffnet" %in% cls) {

    # Checking sizes
    if (!inherits(S, "list")) S <- lapply(1:nslices(graph), function(x) S)
    if (inherits(X, "character")) {
      X <- lapply(1:nslices(graph), function(x) cbind(as.matrix(graph[[X]])))
    }
    if (!inherits(X, c("list"))) {
      X <- lapply(1:nslices(graph), function(x) as.matrix(X))
    }

    with(graph, Map(function(g, x, s) vertex_mahalanobis_dist(g,x,s), g=graph, x=X, s=S))
  }

  return(ans)
}

#' Non-zero element-wise comparison between two sparse matrices
#'
#' Taking advantage of matrix sparseness, the function only evaluates
#' \code{fun} between pairs of elements of \code{A} and \code{B} where
#' either \code{A} or \code{B} have non-zero values. This can be helpful
#' to implement other binary operators between sparse matrices that may
#' not be implemented in the \pkg{Matrix} package.
#'
#' @param A A matrix of size \code{n*m} of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.
#' @param B A matrix of size \code{n*m} of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.
#' @param fun A function that receives 2 arguments and returns a scalar.
#'
#' @details Instead of comparing element by element, the function
#' loops through each matrix non-zero elements to make the comparisons, which
#' in the case of sparse matrices can be more efficient (faster). Algorithmically
#' it can be described as follows:
#'
#' \preformatted{
#' # Matrix initialization
#' init ans[n,m];
#'
#' # Looping through non-zero elements of A
#' for e_A in E_A:
#'   ans[e_A] = fun(A[e_A], B[e_A])
#'
#' # Looping through non-zero elements of B and applying the function
#' # in e_B only if it was not applied while looping in E_A.
#' for e_B in E_B:
#'   if (ans[e_B] == Empty)
#'     ans[e_B] = fun(A[e_B], B[e_B])
#'
#' }
#'
#' \code{compare_matrix} is just an alias for \code{matrix_compare}.
#'
#' @return An object of class \code{dgCMatrix} of size \code{n*m}.
#' @export
#' @examples
#' # These two should yield the same results -----------------------------------
#'
#' # Creating two random matrices
#' set.seed(89)
#' A <- rgraph_ba(t = 9, m = 4)
#' B <- rgraph_ba(t = 9, m = 4)
#' A;B
#'
#' # Comparing
#' ans0 <- matrix_compare(A,B, function(a,b) (a+b)/2)
#'
#' ans1 <- matrix(0, ncol=10, nrow=10)
#' for (i in 1:10)
#'   for (j in 1:10)
#'     ans1[i,j] <- mean(c(A[i,j], B[i,j]))
#'
#' # Are these equal?
#' all(ans0[] == ans1[]) # Should yield TRUE
#'
# # More elaborated example (speed) -------------------------------------------
#
# set.seed(123123123)
# A <- rgraph_ba(t = 5e3, m = 2)
# B <- rgraph_ba(t = 5e3, m = 2)
#
# Am <- as.matrix(A)
# Bm <- as.matrix(B)
#
# compfun <- function(a,b) {
#   ifelse(a > b, a, b)
# }
#
# t0 <- system.time(matrix_compare(A, B, compfun))
# t1 <- system.time(matrix(ifelse(Am > Bm, Am, Bm), ncol=ncol(Am)))
# t1/t0
#' @aliases binary-functions
#' @family dyadic-level comparison functions
matrix_compare <- function(A, B, fun) {

  # Checking objects class
  if (!inherits(A, c("dgCMatrix")))
    stop("-A- must be a dgCMatrix.")

  if (!inherits(B, c("dgCMatrix")))
    stop("-B- must be a dgCMatrix.")

  if (any(dim(A) != dim(B)))
    stop("-A- and -B- must have the same dimension.")

  matrix_compareCpp(A, B, fun)
}

#' @rdname matrix_compare
#' @export
compare_matrix <- matrix_compare

#' Degree and Time of Adoption Diagnostic
#'
#' Analyzes the correlation between in-degree, out-degree, and time of adoption
#' to identify whether opinion leaders were early adopters (supporters) or late
#' adopters (opposers).
#'
#' @param graph A \code{\link{diffnet}} object or a graph data structure (classes include 
#'   \code{array} (\eqn{n\times n \times T}{n*n*T}), \code{dgCMatrix} (sparse), 
#'   \code{igraph}, etc.; see \link{netdiffuseR-graphs}). 
#' @param toa Integer vector of length \eqn{n} (single behavior) or an \eqn{n\times Q}{n*Q}
#'   matrix (multi-behavior) with times of adoption. Required when \code{graph} is not a \code{diffnet}.
#' @param t0,t1 Optional integer scalars defining the first and last observed
#'   periods. If missing and \code{toa} is provided, \code{t0} defaults to 1
#'   and \code{t1} to \code{max(toa, na.rm=TRUE)}.
#' @param name Optional character scalars used only when coercing
#'   inputs into a \code{diffnet} object (passed to \code{new_diffnet}).
#' @param behavior Which behaviors to include when \code{toa} is a matrix (multi-diffusion).
#'   Can be \code{NULL} (all), a numeric index vector, or a character vector matching \code{colnames(toa)}.
#' @param combine Character scalar. How to combine multiple behaviors when \code{toa} is a matrix:
#'   \code{"none"} (analyze each behavior separately), \code{"pooled"} (stack rows across behaviors),
#'   \code{"average"} (per-actor mean of TOA across selected behaviors), or 
#'   \code{"earliest"} (per-actor minimum TOA). Ignored for single-behavior.
#' @param min_adopters Integer scalar. Minimum number of adopters required to compute correlations
#'   for any analysis cell (default 3).
#' @param degree_strategy Character scalar. How to aggregate degree measures across
#'   time periods: \code{"mean"} (default), \code{"first"}, or \code{"last"}.
#' @param bootstrap Logical scalar. Whether to compute bootstrap confidence intervals.
#' @param R Integer scalar. Number of bootstrap replicates (default 1000).
#' @param conf.level Numeric scalar. Confidence level for bootstrap intervals (default 0.95).
#' @param valued Logical scalar. Whether to use edge weights in degree calculations.
#' @param ... Additional arguments passed on when coercing to \code{diffnet}.
#'
#' @details
#' This diagnostic function computes correlations between degree centrality measures
#' (in-degree and out-degree) and time of adoption. Positive correlations suggest
#' that central actors (opinion leaders) adopted early, while negative correlations 
#' suggest they adopted late.
#'
#' The degree measures can be aggregated across time periods in different ways:
#' \itemize{
#'   \item \code{"mean"}: Average degree across all time periods
#'   \item \code{"first"}: Degree in the first time period
#'   \item \code{"last"}: Degree in the last time period
#' }
#' 
#' When \code{bootstrap = TRUE}, the function uses the \pkg{boot} package to
#' compute bootstrap confidence intervals for the correlations.
#' 
#' When \code{toa} is a matrix (multi-diffusion), degree vectors are computed once and
#' reused; the time of adoption is combined according to \code{combine}:
#' \itemize{
#'   \item \code{"none"}: computes separate results per behavior (see Value).
#'   \item \code{"pooled"}: stacks (actor, behavior) rows for adopters and runs a single analysis.
#'   \item \code{"average"}: one row per actor using the mean TOA of adopted behaviors.
#'   \item \code{"earliest"}: one row per actor using the minimum TOA of adopted behaviors.
#' }
#'
#' @return When analyzing a single behavior (or when \code{combine!="none"}), a list with:
#' \item{correlations}{Named numeric vector with correlations between in-degree/out-degree and time of adoption}
#' \item{bootstrap}{List with bootstrap results when \code{bootstrap = TRUE}, otherwise \code{NULL}}
#' \item{call}{The matched call}
#' \item{degree_strategy}{The degree aggregation strategy used}
#' \item{sample_size}{Number of rows included in the analysis (adopter rows)}
#' \item{combine}{\code{NULL} for single-behavior; otherwise the combination rule used.}
#'
#' When \code{combine="none"} with multiple behaviors, returns the same structure, except:
#' \itemize{
#'   \item \code{correlations} is a \eqn{2\times Q^*}{2 x Q*} matrix with rows \code{c("indegree_toa","outdegree_toa")}
#'     and one column per analyzed behavior.
#'   \item \code{bootstrap} is a named list with one entry per behavior (each like the single-behavior case), or \code{NULL} if \code{bootstrap=FALSE}.
#'   \item \code{sample_size} is an integer vector named by behavior.
#'   \item \code{combine} is \code{"none"}.
#' }
#'
#' @examples
#' # Basic usage with Korean Family Planning data
#' data(kfamilyDiffNet)
#' result_basics <- degree_adoption_diagnostic(kfamilyDiffNet, bootstrap = FALSE)
#' print(result_basics)
#' 
#' # With bootstrap confidence intervals
#' result_boot <- degree_adoption_diagnostic(kfamilyDiffNet)
#' print(result_boot)
#' 
#' # Different degree aggregation strategies
#' result_first <- degree_adoption_diagnostic(kfamilyDiffNet, degree_strategy = "first")
#' result_last  <- degree_adoption_diagnostic(kfamilyDiffNet, degree_strategy = "last")
#' 
#' # Multi-diffusion (toy) ----------------------------------------------------
#' \dontrun{
#' set.seed(999)
#' n <- 40; t <- 5; q <- 2
#' garr <- rgraph_ws(n, t, p=.3)
#' diffnet_multi <- rdiffnet(seed.graph = garr, t = t, seed.p.adopt = rep(list(0.1), q))
#'
#' # pooled (one combined analysis)
#' degree_adoption_diagnostic(diffnet_multi, combine = "pooled", bootstrap = FALSE)
#'
#' # per-behavior (matrix of correlations; one column per behavior)
#' degree_adoption_diagnostic(diffnet_multi, combine = "none", bootstrap = FALSE)
#' }
#'
#' @seealso \code{\link{dgr}}, \code{\link{diffreg}}, \code{\link{exposure}}
#' @family statistics
#' @export
degree_adoption_diagnostic <- function(
  graph,
  degree_strategy = c("mean", "first", "last"),
  bootstrap = TRUE,
  R = 1000,
  conf.level = 0.95,
  valued = getOption("diffnet.valued", FALSE),
  toa = NULL,
  t0 = NULL,
  t1 = NULL,
  name = "A diffusion network",
  behavior = NULL,
  combine = c("none","pooled","average","earliest"),
  min_adopters = 3L,
  ...
)
{
  # Coerce to diffnet if needed ----------------------------------------------
  if (!inherits(graph, "diffnet")) {
    if (is.null(toa))
      stop("When 'graph' is not a diffnet, argument 'toa' must be provided.", call. = FALSE)
    dn <- .ndd_coerce_to_diffnet(
      graph, toa = toa, t0 = t0, t1 = t1,
      name = name, behavior = behavior,
      valued = valued, ...
    )
    graph <- dn
  }

  degree_strategy <- match.arg(degree_strategy)
  combine <- match.arg(combine)

  if (!is.logical(bootstrap) || length(bootstrap) != 1)
    stop("'bootstrap' must be a logical scalar.", call. = FALSE)
  if (!is.numeric(R) || length(R) != 1 || R < 1)
    stop("'R' must be a positive integer.", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(min_adopters) || length(min_adopters) != 1 || min_adopters < 1)
    stop("'min_adopters' must be a positive integer.", call. = FALSE)

  # Extract TOA (vector or matrix) -------------------------------------------
  toa_obj <- graph$toa
  is_multi <- is.matrix(toa_obj)

  # Degree measures (aggregated once) ----------------------------------------
  indegree_full  <- dgr(graph, cmode = "indegree", valued = valued)
  outdegree_full <- dgr(graph, cmode = "outdegree", valued = valued)

  # Detect undirected: either meta flag or indegree==outdegree everywhere
  undirected_flag <- isTRUE(graph$meta$undirected)
  if (!undirected_flag) {
    if (identical(dim(indegree_full), dim(outdegree_full))) {
      undirected_flag <- all(indegree_full == outdegree_full, na.rm = TRUE)
    }
  }

  # Aggregate degrees according to strategy
  if (degree_strategy == "mean") {
    indegree_agg  <- rowMeans(indegree_full, na.rm = TRUE)
    outdegree_agg <- rowMeans(outdegree_full, na.rm = TRUE)
  } else if (degree_strategy == "first") {
    indegree_agg  <- indegree_full[, 1]
    outdegree_agg <- outdegree_full[, 1]
  } else { # "last"
    indegree_agg  <- indegree_full[, ncol(indegree_full)]
    outdegree_agg <- outdegree_full[, ncol(outdegree_full)]
  }

  # If undirected, collapse to a single degree vector (use "degree")
  if (undirected_flag) {
    degree_full <- dgr(graph, cmode = "degree", valued = valued)
    degree_agg <- switch(degree_strategy,
      mean  = rowMeans(degree_full, na.rm = TRUE),
      first = degree_full[, 1],
      last  = degree_full[, ncol(degree_full)]
    )
    indegree_agg  <- degree_agg
    outdegree_agg <- degree_agg
  }

  # Safe correlation (avoid warnings on zero variance)
  cor_safely <- function(x, y) {
    x <- as.numeric(x); y <- as.numeric(y)
    if (sum(stats::complete.cases(x, y)) < 2) return(NA_real_)
    sx <- stats::sd(x, na.rm = TRUE); sy <- stats::sd(y, na.rm = TRUE)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
    stats::cor(x, y, use = "complete.obs")
  }

  # Helper: bootstrap one dataset (3 cols: indeg, outdeg, toa)
  boot_one <- function(boot_data) {
    boot_cor <- function(data, indices) {
      d <- data[indices, , drop = FALSE]
      r1 <- cor_safely(d[, 1], d[, 3])
      r2 <- cor_safely(d[, 2], d[, 3])
      c(r1, r2)
    }
    boot_result <- boot::boot(boot_data, boot_cor, R = R)

    # Post-process allowing NAs in replicates
    mk <- function(idx) {
      t0 <- boot_result$t0[idx]
      tv <- boot_result$t[, idx]
      tv <- tv[is.finite(tv)]
      out <- list(correlation = t0)
      if (length(tv) >= 1) out$bias <- mean(tv) - t0
      if (length(tv) >= 2) out$std_error <- stats::sd(tv)
      if (length(tv) >= 3) {
        alpha <- (1 - conf.level)/2
        qs <- stats::quantile(tv, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE, type = 7)
        out$conf_int   <- as.numeric(qs)
        out$conf_level <- conf.level
      }
      out
    }

    list(
      R = R,
      boot_object = boot_result,
      indegree = mk(1),
      outdegree = mk(2)
    )
  }

  # ---------- SINGLE BEHAVIOR PATH ----------
  if (!is_multi) {
    adopters <- !is.na(toa_obj)
    if (sum(adopters) < min_adopters)
      stop("At least ", min_adopters, " adopters are required for correlation analysis.", call. = FALSE)

    toa_adopters       <- toa_obj[adopters]
    indegree_adopters  <- indegree_agg[adopters]
    outdegree_adopters <- outdegree_agg[adopters]

    cor_indegree_toa  <- cor_safely(indegree_adopters, toa_adopters)
    cor_outdegree_toa <- cor_safely(outdegree_adopters, toa_adopters)

    correlations <- c(indegree_toa = cor_indegree_toa,
                      outdegree_toa = cor_outdegree_toa)

    bootstrap_results <- NULL
    if (bootstrap) {
      boot_data <- cbind(indegree_adopters, outdegree_adopters, toa_adopters)
      boot_data <- boot_data[stats::complete.cases(boot_data), , drop = FALSE]
      if (nrow(boot_data) >= min_adopters) {
        bootstrap_results <- boot_one(boot_data)
      } else {
        warning("Insufficient complete cases for bootstrap analysis.", call. = FALSE)
      }
    }

    return(structure(list(
      correlations     = correlations,
      bootstrap        = bootstrap_results,
      call             = match.call(),
      degree_strategy  = degree_strategy,
      sample_size      = sum(adopters),
      combine          = NULL,
      undirected       = undirected_flag
    ), class = "degree_adoption_diagnostic"))
  }

  # ---------- MULTI-BEHAVIOR PATH ----------
  Q <- ncol(toa_obj)
  beh_names <- colnames(toa_obj)
  if (is.null(beh_names)) beh_names <- paste0("B", seq_len(Q))

  # Resolve behavior subset
  if (is.null(behavior)) {
    q_set <- seq_len(Q)
  } else if (is.numeric(behavior)) {
    q_set <- as.integer(behavior)
    if (any(q_set < 1 | q_set > Q)) stop("Some 'behavior' indices are out of range 1..", Q, ".")
  } else if (is.character(behavior)) {
    mm <- match(behavior, beh_names)
    if (anyNA(mm)) stop("Some 'behavior' names not found in TOA colnames.")
    q_set <- mm
  } else {
    stop("'behavior' must be NULL, numeric indices, or character names.")
  }
  beh_names <- beh_names[q_set]

  # Build analysis dataset depending on 'combine'
  make_boot_data <- switch(
    combine,
    pooled = function() {
      rows <- lapply(q_set, function(q) {
        a <- which(!is.na(toa_obj[, q]))
        if (length(a) == 0L) return(NULL)
        cbind(indegree_agg[a], outdegree_agg[a], toa_obj[a, q])
      })
      do.call(rbind, rows)
    },
    average = function() {
      toa_mean <- rowMeans(toa_obj[, q_set, drop = FALSE], na.rm = TRUE)
      adopters_any <- rowSums(!is.na(toa_obj[, q_set, drop = FALSE])) > 0
      a <- which(adopters_any & !is.nan(toa_mean))
      if (!length(a)) return(NULL)
      cbind(indegree_agg[a], outdegree_agg[a], toa_mean[a])
    },
    earliest = function() {
      toa_min <- apply(toa_obj[, q_set, drop = FALSE], 1, function(v) {
        v <- v[!is.na(v)]
        if (!length(v)) return(Inf)
        min(v)
      })
      adopters_any <- is.finite(toa_min)
      a <- which(adopters_any)
      if (!length(a)) return(NULL)
      cbind(indegree_agg[a], outdegree_agg[a], toa_min[a])
    },
    none = NULL
  )

  if (combine != "none") {
    boot_data <- make_boot_data()
    if (is.null(boot_data) || nrow(boot_data) < min_adopters)
      stop("Not enough adopter rows after combining behaviors (combine='", combine, "').", call. = FALSE)

    cor_indegree_toa  <- cor_safely(boot_data[,1], boot_data[,3])
    cor_outdegree_toa <- cor_safely(boot_data[,2], boot_data[,3])
    correlations <- c(indegree_toa = cor_indegree_toa,
                      outdegree_toa = cor_outdegree_toa)

    bootstrap_results <- NULL
    if (bootstrap) {
      boot_data2 <- boot_data[stats::complete.cases(boot_data), , drop = FALSE]
      if (nrow(boot_data2) >= min_adopters) {
        bootstrap_results <- boot_one(boot_data2)
      } else {
        warning("Insufficient complete cases for bootstrap analysis.", call. = FALSE)
      }
    }

    return(structure(list(
      correlations     = correlations,
      bootstrap        = bootstrap_results,
      call             = match.call(),
      degree_strategy  = degree_strategy,
      sample_size      = nrow(boot_data),
      combine          = combine,
      undirected       = undirected_flag
    ), class = "degree_adoption_diagnostic"))
  }

  # combine == "none": per-behavior results ---------------------------------
  K <- length(q_set)
  corr_mat <- matrix(NA_real_, nrow = 2, ncol = K,
                     dimnames = list(c("indegree_toa","outdegree_toa"), beh_names))
  boot_list <- if (bootstrap) vector("list", K) else NULL
  n_vec <- integer(K); names(n_vec) <- beh_names

  for (j in seq_along(q_set)) {
    q <- q_set[j]
    adopters <- !is.na(toa_obj[, q])
    n_vec[j] <- sum(adopters)
    if (n_vec[j] < min_adopters) next

    toa_q <- toa_obj[adopters, q]
    indeg_q <- indegree_agg[adopters]
    outdeg_q <- outdegree_agg[adopters]

    corr_mat["indegree_toa", j]  <- cor_safely(indeg_q, toa_q)
    corr_mat["outdegree_toa", j] <- cor_safely(outdeg_q, toa_q)

    if (bootstrap) {
      bd <- cbind(indeg_q, outdeg_q, toa_q)
      bd <- bd[stats::complete.cases(bd), , drop = FALSE]
      if (nrow(bd) >= min_adopters) boot_list[[j]] <- boot_one(bd)
    }
  }

  structure(list(
    correlations     = corr_mat,
    bootstrap        = boot_list,
    call             = match.call(),
    degree_strategy  = degree_strategy,
    sample_size      = n_vec,
    combine          = "none",
    undirected       = undirected_flag
  ), class = "degree_adoption_diagnostic")
}

# Internal: Coerce to diffnet from common graph inputs ------------------------
.ndd_coerce_to_diffnet <- function(graph, toa, t0 = NULL, t1 = NULL,
                                   name = "A diffusion network",
                                   behavior = NULL,
                                   valued = getOption("diffnet.valued", FALSE),
                                   ...) {
  # infer periods from TOA when not provided
  if (is.null(t0)) t0 <- 1L
  if (is.null(t1)) t1 <- max(toa, na.rm = TRUE)
  Tlen_target <- max(1L, as.integer(t1 - t0 + 1L))

  # If behavior not provided, derive a sensible default
  if (is.null(behavior)) {
    if (is.matrix(toa)) {
      behavior <- if (!is.null(colnames(toa))) colnames(toa) else paste0("behavior_", seq_len(ncol(toa)))
    } else {
      behavior <- "Innovation"
    }
  }

  # Helper: coerce a single graph slice to dgCMatrix
  as_slice_dgC <- function(g) {
    if (inherits(g, "dgCMatrix")) {
      G <- g
    } else if (inherits(g, "matrix")) {
      G <- methods::as(g, "dgCMatrix")
    } else if (inherits(g, "igraph")) {
      G <- methods::as(igraph::as_adjacency_matrix(g, sparse = TRUE), "dgCMatrix")
    } else if (inherits(g, "network")) {
      tmp <- as_generic_graph.network(g)
      G <- tmp$graph[[1]]
    } else {
      stop("Unsupported graph slice class for coercion.")
    }
    if (!valued) G@x <- rep(1, length(G@x))
    G
  }

  # Build list of time slices (each a dgCMatrix)
  if (inherits(graph, "list")) {
    Glist <- lapply(graph, as_slice_dgC)
  } else if (inherits(graph, "array")) {
    if (length(dim(graph)) != 3L || dim(graph)[1] != dim(graph)[2])
      stop("When 'graph' is an array it must be n x n x T.")
    Glist <- lapply(seq_len(dim(graph)[3]), function(i) methods::as(graph[, , i], "dgCMatrix"))
    if (!valued) Glist <- lapply(Glist, function(G) { G@x <- rep(1, length(G@x)); G })
  } else {
    # Single static graph: make one slice
    Glist <- list(as_slice_dgC(graph))
  }

  # If only one slice provided but multiple periods implied by TOA, repeat it
  if (length(Glist) == 1L && Tlen_target > 1L) {
    Glist <- rep(Glist, Tlen_target)
    names(Glist) <- as.character(seq.int(t0, t1))
  }

  # Sanity check
  if (!all(vapply(Glist, function(x) inherits(x, "dgCMatrix"), logical(1))))
    stop("All graph slices must be coercible to 'dgCMatrix'.")

  # Align vertex names with TOA (vector or matrix)
  rn <- rownames(Glist[[1]])
  if (length(rn)) {
    toa_names <- if (is.matrix(toa)) rownames(toa) else names(toa)
    if (!is.null(toa_names) && !identical(toa_names, rn)) {
      ord <- match(rn, toa_names)
      if (anyNA(ord)) stop("Names in 'toa' do not match graph vertex names.")
      if (is.matrix(toa)) toa <- toa[ord, , drop = FALSE] else toa <- toa[ord]
    }
  }

  new_diffnet(graph = Glist, toa = toa, t0 = t0, t1 = t1,
              name = name, behavior = behavior)
}

#' @export
print.degree_adoption_diagnostic <- function(x, ...) {
  cat("Degree and Time of Adoption Diagnostic\n")
  cat("======================================\n\n")

  cat("Degree aggregation strategy:", x$degree_strategy, "\n")

  # Sample size can be a scalar or a vector (per behavior when combine='none')
  if (length(x$sample_size) == 1L) {
    cat("Sample size (adopters only):", x$sample_size, "\n\n")
  } else {
    cat("Sample size (adopters only): total =", sum(x$sample_size, na.rm = TRUE), "\n")
    cat("  By behavior:\n")
    beh_names <- names(x$sample_size)
    for (j in seq_along(x$sample_size)) {
      cat(sprintf("    - %s: %d\n", if (length(beh_names)) beh_names[j] else paste0("B", j), x$sample_size[j]))
    }
    cat("\n")
  }

  undirected <- isTRUE(x$undirected)

  # Helper used in both modes ------------------------------------------------
  explain <- function(label, r, ci, lvl_arg = NA_real_, thr = 0.10) {
    # Handle NA correlations gracefully
    if (is.na(r)) {
      if (is.null(ci)) {
        cat(sprintf(
          "  %s: Weak relationship between centrality and adoption timing:\n             r is NA; no CI.\n",
          label
        ))
      } else {
        lvl_local <- if (!is.na(lvl_arg)) lvl_arg else 95
        ci_includes_zero <- (length(ci) >= 2) && is.finite(ci[1]) && is.finite(ci[2]) && (ci[1] <= 0 && ci[2] >= 0)
        cat(sprintf(
          "  %s: Weak relationship between centrality and adoption timing; %s statistically supported:\n             r is NA; CI (%.1f%%) %s 0.\n",
          label,
          if (ci_includes_zero) "NOT" else "",
          lvl_local,
          if (ci_includes_zero) "includes" else "excludes"
        ))
      }
      return(invisible())
    }

    abs_big <- abs(r) > thr

    if (is.null(ci)) {
      if (!abs_big) {
        cat(sprintf("  %s: Weak relationship between %s and adoption timing:\n             |r| \u2264 %.1f; no CI.\n",
                    label,
                    if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                    thr))
      } else if (r > 0) {
        cat(sprintf("  %s: Central actors (high %s) tended to adopt early (supporters):\n             |r| > %.1f; no CI.\n",
                    label,
                    if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                    thr))
      } else {
        cat(sprintf("  %s: Central actors (high %s) tended to adopt late (opposers):\n             |r| > %.1f; no CI.\n",
                    label,
                    if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                    thr))
      }
      return(invisible())
    }

    # CI available
    ci_includes_zero <- (length(ci) >= 2) && is.finite(ci[1]) && is.finite(ci[2]) && (ci[1] <= 0 && ci[2] >= 0)
    lvl_local <- if (!is.na(lvl_arg)) lvl_arg else 95

    if (!abs_big && !ci_includes_zero) {
      cat(sprintf("  %s: Weak relationship between centrality and adoption timing; statistically supported:\n             |r| \u2264 %.1f; CI (%.1f%%) excludes 0.\n",
                  label, 0.10, lvl_local))
    } else if (!abs_big && ci_includes_zero) {
      cat(sprintf("  %s: Weak relationship between centrality and adoption timing; NOT statistically supported:\n             |r| \u2264 %.1f; CI (%.1f%%) includes 0.\n",
                  label, 0.10, lvl_local))
    } else if (abs_big && !ci_includes_zero && r > 0) {
      cat(sprintf("  %s: Central actors (high %s) tended to adopt early (supporters); statistically supported:\n             |r| > %.1f; CI (%.1f%%) excludes 0.\n",
                  label,
                  if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                  0.10, lvl_local))
    } else if (abs_big && !ci_includes_zero && r < 0) {
      cat(sprintf("  %s: Central actors (high %s) tended to adopt late (opposers); statistically supported:\n             |r| > %.1f; CI (%.1f%%) excludes 0.\n",
                  label,
                  if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                  0.10, lvl_local))
    } else if (abs_big && ci_includes_zero && r > 0) {
      cat(sprintf("  %s: Central actors (high %s) tended to adopt early (supporters); NOT statistically supported:\n             |r| > %.1f; CI (%.1f%%) includes 0.\n",
                  label,
                  if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                  0.10, lvl_local))
    } else if (abs_big && ci_includes_zero && r < 0) {
      cat(sprintf("  %s: Central actors (high %s) tended to adopt late (opposers); NOT statistically supported:\n             |r| > %.1f; CI (%.1f%%) includes 0.\n",
                  label,
                  if (label == "In-degree") "in-degree" else if (label == "Out-degree") "out-degree" else "degree",
                  0.10, lvl_local))
    }
  }

  # Branch on combine mode ---------------------------------------------------
  if (is.null(x$combine) || x$combine != "none") {
    # -------- single behavior OR combined pooled/average/earliest ----------
    cat("\nCorrelations:\n")
    if (undirected) {
      rdeg <- x$correlations[["indegree_toa"]] # same as outdegree_toa
      cat(sprintf("  Degree     - Time of Adoption: %6.3f\n", rdeg))
    } else {
      cat(sprintf("  In-degree  - Time of Adoption: %6.3f\n", x$correlations[["indegree_toa"]]))
      cat(sprintf("  Out-degree - Time of Adoption: %6.3f\n", x$correlations[["outdegree_toa"]]))
    }

    if (!is.null(x$bootstrap)) {
      cat("\nBootstrap Results:\n")
      cat(sprintf("  Replicates: %d\n", x$bootstrap$R))
      if (!is.null(x$bootstrap$indegree$conf_level))
        cat(sprintf("  Confidence level: %.1f%%\n\n", x$bootstrap$indegree$conf_level * 100))

      if (!is.null(x$bootstrap$indegree$conf_int)) {
        if (undirected) {
          cat("  Degree correlation:\n")
          cat(sprintf("    Estimate: %6.3f\n", x$bootstrap$indegree$correlation))
          if (!is.null(x$bootstrap$indegree$bias))      cat(sprintf("    Bias: %6.3f\n", x$bootstrap$indegree$bias))
          if (!is.null(x$bootstrap$indegree$std_error)) cat(sprintf("    Std. Error: %6.3f\n", x$bootstrap$indegree$std_error))
          cat(sprintf("    %.1f%% CI: [%6.3f, %6.3f]\n\n",
                      x$bootstrap$indegree$conf_level * 100,
                      x$bootstrap$indegree$conf_int[1], x$bootstrap$indegree$conf_int[2]))
        } else {
          cat("  In-degree correlation:\n")
          cat(sprintf("    Estimate: %6.3f\n", x$bootstrap$indegree$correlation))
          if (!is.null(x$bootstrap$indegree$bias))      cat(sprintf("    Bias: %6.3f\n", x$bootstrap$indegree$bias))
          if (!is.null(x$bootstrap$indegree$std_error)) cat(sprintf("    Std. Error: %6.3f\n", x$bootstrap$indegree$std_error))
          cat(sprintf("    %.1f%% CI: [%6.3f, %6.3f]\n\n",
                      x$bootstrap$indegree$conf_level * 100,
                      x$bootstrap$indegree$conf_int[1], x$bootstrap$indegree$conf_int[2]))

          if (!is.null(x$bootstrap$outdegree$conf_int)) {
            cat("  Out-degree correlation:\n")
            cat(sprintf("    Estimate: %6.3f\n", x$bootstrap$outdegree$correlation))
            if (!is.null(x$bootstrap$outdegree$bias))      cat(sprintf("    Bias: %6.3f\n", x$bootstrap$outdegree$bias))
            if (!is.null(x$bootstrap$outdegree$std_error)) cat(sprintf("    Std. Error: %6.3f\n", x$bootstrap$outdegree$std_error))
            cat(sprintf("    %.1f%% CI: [%6.3f, %6.3f]\n\n",
                        x$bootstrap$outdegree$conf_level * 100,
                        x$bootstrap$outdegree$conf_int[1], x$bootstrap$outdegree$conf_int[2]))
          }
        }
      } else {
        if (undirected) {
          cat("  Bootstrap estimates (CI unavailable):\n")
          cat(sprintf("    Degree: Est=%6.3f%s%s\n",
                      x$bootstrap$indegree$correlation,
                      if (!is.null(x$bootstrap$indegree$bias)) paste0(", Bias=", sprintf("%6.3f", x$bootstrap$indegree$bias)) else "",
                      if (!is.null(x$bootstrap$indegree$std_error)) paste0(", SE=", sprintf("%6.3f", x$bootstrap$indegree$std_error)) else ""))
          cat("\n")
        } else {
          cat("  Bootstrap estimates (CI unavailable):\n")
          cat(sprintf("    In-degree: Est=%6.3f%s%s\n",
                      x$bootstrap$indegree$correlation,
                      if (!is.null(x$bootstrap$indegree$bias)) paste0(", Bias=", sprintf("%6.3f", x$bootstrap$indegree$bias)) else "",
                      if (!is.null(x$bootstrap$indegree$std_error)) paste0(", SE=", sprintf("%6.3f", x$bootstrap$indegree$std_error)) else ""))
          cat(sprintf("    Out-degree: Est=%6.3f%s%s\n\n",
                      x$bootstrap$outdegree$correlation,
                      if (!is.null(x$bootstrap$outdegree$bias)) paste0(", Bias=", sprintf("%6.3f", x$bootstrap$outdegree$bias)) else "",
                      if (!is.null(x$bootstrap$outdegree$std_error)) paste0(", SE=", sprintf("%6.3f", x$bootstrap$outdegree$std_error)) else ""))
        }
      }
    }

    cat("Interpretation:\n")
    if (undirected) {
      rdeg <- x$correlations[["indegree_toa"]]
      has_boot <- !is.null(x$bootstrap)
      deg_ci <- if (has_boot && !is.null(x$bootstrap$indegree$conf_int)) x$bootstrap$indegree$conf_int else NULL
      lvl <- if (has_boot && !is.null(x$bootstrap$indegree$conf_level)) x$bootstrap$indegree$conf_level * 100 else NA_real_
      explain("Degree", rdeg, deg_ci, lvl_arg = lvl)
    } else {
      indeg <- x$correlations[["indegree_toa"]]
      outdeg <- x$correlations[["outdegree_toa"]]
      has_boot <- !is.null(x$bootstrap)
      indeg_ci <- if (has_boot && !is.null(x$bootstrap$indegree$conf_int)) x$bootstrap$indegree$conf_int else NULL
      outdeg_ci <- if (has_boot && !is.null(x$bootstrap$outdegree$conf_int)) x$bootstrap$outdegree$conf_int else NULL
      lvl <- if (has_boot && !is.null(x$bootstrap$indegree$conf_level)) x$bootstrap$indegree$conf_level * 100 else NA_real_
      explain("In-degree",  indeg, indeg_ci, lvl_arg = lvl)
      explain("Out-degree", outdeg, outdeg_ci, lvl_arg = lvl)
    }
    invisible(x)
    return(invisible(x))
  }

  # -------------------------- combine == "none" -----------------------------
  cat("Combination mode:", "none (per behavior)\n\n")

  # Correlations section
  cat("Correlations (per behavior):\n")
  beh_names <- colnames(x$correlations)
  for (j in seq_len(ncol(x$correlations))) {
    bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
    cat(sprintf("  [%s]\n", bname))
    if (undirected) {
      rdeg <- x$correlations["indegree_toa", j]
      cat(sprintf("    Degree     - Time of Adoption: %6.3f\n", rdeg))
    } else {
      r_in  <- x$correlations["indegree_toa", j]
      r_out <- x$correlations["outdegree_toa", j]
      cat(sprintf("    In-degree  - Time of Adoption: %6.3f\n", r_in))
      cat(sprintf("    Out-degree - Time of Adoption: %6.3f\n", r_out))
    }
  }
  cat("\n")

  # Bootstrap section (if any)
  if (!is.null(x$bootstrap)) {
    any_boot <- any(vapply(x$bootstrap, function(b) !is.null(b), logical(1)))
    if (any_boot) {
      cat("Bootstrap Results (per behavior):\n")
      for (j in seq_len(ncol(x$correlations))) {
        bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
        b <- x$bootstrap[[j]]
        if (is.null(b)) next
        cat(sprintf("  [%s]\n", bname))
        if (!is.null(b$R)) cat(sprintf("    Replicates: %d\n", b$R))
        if (!is.null(b$indegree$conf_level)) cat(sprintf("    Confidence level: %.1f%%\n", b$indegree$conf_level * 100))
        if (!is.null(b$indegree$conf_int)) {
          if (undirected) {
            cat("    Degree correlation:\n")
            cat(sprintf("      Estimate: %6.3f\n", b$indegree$correlation))
            if (!is.null(b$indegree$bias))      cat(sprintf("      Bias: %6.3f\n", b$indegree$bias))
            if (!is.null(b$indegree$std_error)) cat(sprintf("      Std. Error: %6.3f\n", b$indegree$std_error))
            cat(sprintf("      %.1f%% CI: [%6.3f, %6.3f]\n",
                        b$indegree$conf_level * 100, b$indegree$conf_int[1], b$indegree$conf_int[2]))
          } else {
            cat("    In-degree correlation:\n")
            cat(sprintf("      Estimate: %6.3f\n", b$indegree$correlation))
            if (!is.null(b$indegree$bias))      cat(sprintf("      Bias: %6.3f\n", b$indegree$bias))
            if (!is.null(b$indegree$std_error)) cat(sprintf("      Std. Error: %6.3f\n", b$indegree$std_error))
            cat(sprintf("      %.1f%% CI: [%6.3f, %6.3f]\n",
                        b$indegree$conf_level * 100, b$indegree$conf_int[1], b$indegree$conf_int[2]))
          }
        }
        if (!undirected && !is.null(b$outdegree$conf_int)) {
          cat("    Out-degree correlation:\n")
          cat(sprintf("      Estimate: %6.3f\n", b$outdegree$correlation))
          if (!is.null(b$outdegree$bias))      cat(sprintf("      Bias: %6.3f\n", b$outdegree$bias))
          if (!is.null(b$outdegree$std_error)) cat(sprintf("      Std. Error: %6.3f\n", b$outdegree$std_error))
          cat(sprintf("      %.1f%% CI: [%6.3f, %6.3f]\n",
                      b$outdegree$conf_level * 100, b$outdegree$conf_int[1], b$outdegree$conf_int[2]))
        }
      }
      cat("\n")
    }
  }

  # Interpretation (per behavior)
  cat("Interpretation (per behavior):\n")
  for (j in seq_len(ncol(x$correlations))) {
    bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
    if (undirected) {
      rdeg <- x$correlations["indegree_toa", j]
      b <- if (!is.null(x$bootstrap)) x$bootstrap[[j]] else NULL
      deg_ci <- if (!is.null(b) && !is.null(b$indegree$conf_int)) b$indegree$conf_int else NULL
      lvl <- if (!is.null(b) && !is.null(b$indegree$conf_level)) b$indegree$conf_level * 100 else NA_real_
      cat(sprintf("  [%s]\n", bname))
      explain("Degree",  rdeg,  deg_ci, lvl_arg = lvl)
    } else {
      r_in  <- x$correlations["indegree_toa", j]
      r_out <- x$correlations["outdegree_toa", j]
      b <- if (!is.null(x$bootstrap)) x$bootstrap[[j]] else NULL
      indeg_ci <- if (!is.null(b) && !is.null(b$indegree$conf_int)) b$indegree$conf_int else NULL
      outdeg_ci <- if (!is.null(b) && !is.null(b$outdegree$conf_int)) b$outdegree$conf_int else NULL
      lvl <- if (!is.null(b) && !is.null(b$indegree$conf_level)) b$indegree$conf_level * 100 else NA_real_
      cat(sprintf("  [%s]\n", bname))
      explain("In-degree",  r_in,  indeg_ci, lvl_arg = lvl)
      explain("Out-degree", r_out, outdeg_ci, lvl_arg = lvl)
    }
  }

  invisible(x)
}

