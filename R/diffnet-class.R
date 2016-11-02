# Creates n column names
make_col_names <- function(n, is.dynamic, prefix="v.") {
  expression <- paste0("%s%s%0", nchar(n), "d")
  if (is.dynamic) sprintf(expression,prefix,"dyn.", seq_len(n))
  else sprintf(expression,prefix,"static.", seq_len(n))
}

# Checks attributes to be added to a diffnet object
check_as_diffnet_attrs <- function(attrs, meta, is.dynamic, id.and.per.vars=NULL,
                                   attr.class="vertex") {
  # Getting meta
  n    <- meta$n
  nper <- meta$nper

  # In case of merging
  iddf  <- data.frame(id  = meta$ids, `_original_sort`=seq_len(n), check.names = FALSE )
  perdf <- data.frame(per = meta$pers)

  # Vertex Attributes  ---------------------------------------------------------
  if (attr.class == "vertex") {
    if (is.dynamic) {
      # If null
      if (!length(attrs)) {
        attrs <-
          lapply(meta$pers, function(x) {
            as.data.frame(
              matrix(ncol=0, nrow=n, dimnames = list(meta$ids, NULL))
            )
          })

        # Labeling
        names(attrs) <- meta$pers
        return(attrs)
      }


      if (inherits(attrs, "list")) {
        # Must have the same number of elements as time periods
        if (length(attrs) != nper)
          stop("The list -vertex.dyn.attrs-'s length must be equal to the number of slices.",
               " Has ",length(attrs), " and should have ", nper,".")

        # Checking that all are data.frames/matrices or vectors
        isdf  <- sapply(attrs, inherits, what="data.frame")
        ismat <- sapply(attrs, inherits, what="matrix")
        isvec <- sapply(attrs, is.vector)

        test <- which(!isdf & !ismat & !isvec)
        if (length(test))
          stop("Some in the list -vertex.dyn.attrs- are not supported:\n\t:",
               paste0(test, collapse=", "), ".")

        # Further, all should be of the same class
        if ( (any(isdf) + any(ismat) + any(isvec))>1 )
          stop("All elements in the list -vertex.dyn.attrs- should be of the same class.")

        if (all(isdf) | all(ismat)) {
          # Dynamic: List of matrices/data.frames --------------------------------
          # Checking dims
          test <- which(sapply(attrs, nrow) != n)
          if (length(test)) {
            stop("Some matrices/data frames in -vertex.dyn.attrs- don't have n (",
                 n, ") rows:\n\t", paste0(test, collapse = ", "), ".")
          }

          # Checking number of columns
          test <- which(ncol(attrs[[1]]) != sapply(attrs, ncol))
          if (length(test))
            stop("All elements in -vertex.dyn.attrs- must have the same number of columns.")

          # Checking columnnames
          test <- colnames(attrs[[1]])
          test <- if (length(test)) which(sapply(attrs, function(x) all(colnames(x) != test)))
          else which(!!sapply(attrs, function(x) length(colnames(x))))

          if (length(test))
            stop("Not all the matrices/data frames in -vertex.dyn.attrs- have the same colname(.):\n\t",
                 paste0(test, collapse = ", "))

          # Making up names
          cnames <- colnames(attrs[[1]])
          if (!length(cnames))
            cnames <- make_col_names(ncol(attrs[[1]]), TRUE)

          # Returning as data.frame
          attrs <- lapply(attrs, function(x) {
            x <- as.data.frame(x)

            # Checking if must be ordered
            if (length(id.and.per.vars)) {
              x <- merge(iddf, x, by.x="id", by.y=id.and.per.vars[1],
                         sort=FALSE, all.x=TRUE, all.y=FALSE)

              # Sorting names back
              colnames(x)[1] <- id.and.per.vars[1]
              x <- x[,c(cnames, "_original_sort")]

              # Sorting back rows and removing sort column
              x <- x[order(x[["_original_sort"]]),]
              x[["_original_sort"]] <- NULL

              x
            }

            # Naming
            dimnames(x) <- list(meta$ids, cnames)
            x
            })

          # Adding names
          names(attrs) <- meta$pers

          return(attrs)

        } else { # Dynamic: List of vectors --------------------------------------
          # Checking dim
          test <- which(sapply(attrs, length) != n)
          if (length(test))
            stop("Some vectors in -vertex.dyn.attrs- have a different lengths than expected:\n\t",
                 paste0(test, collapse=", "), ".")

          cnames <- make_col_names(1, TRUE)

          # Returning as data.frame
          attrs <- lapply(attrs, function(x) {
            x <- as.data.frame(x)

            dimnames(x) <- list(meta$ids, cnames)
            x
            })

          # Adding names
          names(attrs) <- meta$pers

          return(attrs)
        }

      } else if (inherits(attrs, "data.frame") | inherits(attrs, "matrix")) {
        # Dynamic: Matrix or data.frame n*nper -----------------------------------
        if (nrow(attrs) != n*nper)
          stop("The matrix/data.frame -vertex.dyn.attrs- has incorrect number of rows.",
               " Has ", nrow(attrs), " and should have ", n, "*", nper," = ",
               n*nper, ".")

        # Checking colnames
        cnames <- colnames(attrs)
        if (!length(cnames))
          cnames <- make_col_names(ncol(attrs), TRUE)

        # Coercing into a list of data.frames
        attrs <- lapply(meta$pers, function(x) {

          # Subseting
          x <- if (!length(id.and.per.vars[2])) attrs[((x-1)*n + 1):(x*n),,drop=FALSE]
          else attrs[attrs[,id.and.per.vars[2]] == x,,drop=FALSE]

          x <- as.data.frame(x)

          # Checking if must be ordered
          if (length(id.and.per.vars)) {
            x <- merge(iddf, x, by.x="id", by.y=id.and.per.vars[1],
                       sort=FALSE, all.x=TRUE, all.y=FALSE)

            # Sorting names (back)
            colnames(x)[1] <- id.and.per.vars[1]
            x <- x[,c(cnames, "_original_sort")]

            # Sorting back rows and removing sort column
            x <- x[order(x[["_original_sort"]]),]
            x[["_original_sort"]] <- NULL

            x
          }

          dimnames(x) <- list(meta$ids, cnames)
          x
          }
        )

        # Adding names
        names(attrs) <- meta$pers

        # Returning
        return(attrs)

      } else if (is.vector(attrs)) {
        # Dynamic: vector of size n*nper -----------------------------------------
        if (length(attrs) != n*nper)
          stop("The vector -vertex.dyn.attrs- has incorrect length.",
               " Has ", length(attrs), " and should have ", n, "*", nper," = ",
               n*nper, ".")

        # Colname
        cnames <- make_col_names(1, TRUE)

        # Coercing into a list of data.frames
        attrs <- lapply(meta$pers, function(x) {
          x <- as.data.frame(attrs[((x-1)*n + 1):(x*n)])
          dimnames(x) <- list(meta$ids, cnames)
          x}
        )

        # Adding names
        names(attrs) <- meta$pers

        # Returning
        return(attrs)

      } else
        stop("The class of -vertex.dyn.attrs-, \'", class(attrs),
             "\', is not supported.")

    } else { # Static attributes
      # If null
      if (!length(attrs)) {
        return(as.data.frame(
          matrix(ncol=0, nrow=n, dimnames = list(meta$ids, NULL))
        ))
      }

      if (inherits(attrs, "data.frame") | inherits(attrs, "matrix")) {
        # Static: Data frame or matrix -------------------------------------------
        if (nrow(attrs) != n)
          stop("The matrix/data.frame -vertex.static.attrs- has incorrect number of rows.",
               " Has ", nrow(attrs), " and should have ", n, ".")

        # If we need to sort using id
        attrs <- as.data.frame(attrs)
        cnames <- colnames(attrs)
        if (length(id.and.per.vars)) {
          attrs <- merge(iddf, attrs, by.x="id", by.y=id.and.per.vars[1],
                     sort=FALSE, all.x=TRUE, all.y=FALSE)

          # Sorting names back
          colnames(attrs)[1] <- id.and.per.vars[1]
          attrs <- attrs[,c(cnames, "_original_sort")]

          # Sorting back rows and removing sort column
          attrs <- attrs[order(attrs[["_original_sort"]]),]
          attrs[["_original_sort"]] <- NULL
        }

        # Checking colnames
        if (!length(cnames))
          cnames <- make_col_names(ncol(attrs), FALSE)

        # Returning
        dimnames(attrs) <- list(meta$ids, cnames)
        return(attrs)

      } else if (is.vector(attrs)) {
        # Static: Vector -------------------------------------------------------
        if (length(attrs) != n)
          stop("The vector -vertex.static.attrs- has incorrect length.",
               " Has ", length(attrs), " and should have ", n, ".")

        # Returning
        attrs <- as.data.frame(attrs)
        dimnames(attrs) <- list(meta$ids, make_col_names(1, FALSE))

        return(attrs)

      } else
        stop("The class of -vertex.static.attrs-, \'", class(attrs),
             "\', is not supported.")

    }
  } else if (attr.class == "graph") {
    # if Null
    if (!length(attrs)) {
      return(as.data.frame(
          matrix(ncol=0, nrow=nper, dimnames = list(meta$pers, NULL))
        ))
    }

    # Graph attributes ---------------------------------------------------------
    if (inherits(attrs, "data.frame") | inherits(attrs, "matrix")) {
      # Checking number of rows
      if (nrow(attrs) != nper)
        stop("The matrix/data.frame -graph.attrs- has incorrect number of rows.",
             " Has ", nrow(attrs), " and should have ", nper, ".")

      # # Checking the order of the data
      # attrs <- as.data.frame(attrs)
      # cnames <- colnames(attrs)
      # if (length(id.and.per.vars[2])) {
      #   attrs <- merge(perdf, attrs, by.x="per", by.y=id.and.per.vars[2],
      #                  all.x=TRUE, all.y=FALSE, sort=TRUE)
      #   # Sorting names back
      #   colnames(attrs)[1] <- id.and.per.vars[1]
      #   attrs <- attrs[,cnames]
      # }

      # Checking colnames
      if (!length(cnames))
        cnames <- make_col_names(ncol(attrs), TRUE, prefix = "graph")

      # Returning
      dimnames(attrs) <- list(meta$pers, cnames)

      return(attrs)

    } else if (is.vector(attrs) && inherits(attrs, c("numeric", "character", "raw", "integer"))) {

      # Checking the length
      if (length(attrs) != nper)
        stop("The vector -graph.attrs- has incorrect length.",
             " Has ", length(attrs), " and should have ", nper, ".")

      # Returning
      attrs <- as.data.frame(attrs)
      dimnames(attrs) <- list(meta$pers, make_col_names(1, TRUE, prefix = "graph"))

      return(attrs)

    } else
      stop("The class of -graph.attrs-, \'", class(attrs),
           "\', is not supported.")
  }
}

#' Creates a \code{diffnet} class object
#'
#' \code{diffnet} objects contain difussion networks. With adjacency
#' matrices and time of adoption (toa) vector as its main components, most of the
#' package's functions have methods for this class of objects.
#'
#' @templateVar dynamic TRUE
#' @templateVar undirected TRUE
#' @templateVar self TRUE
#' @templateVar valued TRUE
#' @templateVar multiple TRUE
#' @template graph_template
#' @param gmode Character scalar. Passed to \code{\link[sna:gplot]{sna::gplot}.}
#' @param toa Numeric vector of size \eqn{n}. Times of adoption.
#' @param t0 Integer scalar. Passed to \code{\link{toa_mat}}.
#' @param t1 Integer scalar. Passed to \code{\link{toa_mat}}.
#' @param ... In the case of \code{plot}, further arguments passed to \code{gplot}, for
#' \code{summary}, further arguments to be passed to
#' \code{\link[igraph:distances]{igraph::distances}}, otherwise is ignored.
#' @param x A \code{diffnet} object.
#' @param object A \code{diffnet} object.
#' @param y Ignored.
#' @param t Integer scalar indicating the time slice to plot.
#' @param displaylabels Logical scalar. When TRUE, \code{plot} shows vertex labels.
#' @param vertex.col Character scalar/vector. Color of the vertices.
#' @param vertex.cex Numeric scalar/vector. Size of the vertices.
#' @param edge.col Character scalar/vector. Color of the edges.
#' @param mode Character scalar. In the case of \code{plot}, passed to
#' \code{\link[sna:gplot]{gplot}}. Otherwise, in the case of summary, passed to
#' \code{\link[igraph:distances]{igraph::distances}}.
#' @param layout.par Layout parameters (see details).
#' @param main Character. A title template to be passed to sprintf.
#' @param i Indices specifying elements to replace. See \code{\link[base:Extract]{Extract}}.
#' @param value In the case of \code{diffnet.toa}, replacement, otherwise see below.
#' @param vertex.dyn.attrs Vertices dynamic attributes (see details).
#' @param vertex.static.attrs Vertices static attributes (see details).
#' @param graph.attrs Graph dynamic attributes (not supported yet).
#' @param id.and.per.vars A character vector of length 2. Optionally specified to check the
#' order of the rows in the attribute data.
#' @param slices Either an integer or character vector. While integer vectors are used as
#' indexes, character vectors are used jointly with the time period labels.
#' @param element Character vector/scalar. Indicates what to retrieve/alter.
#' @param attr.class Character vector/scalar. Indicates the class of the attribute, either dynamic (\code{"dyn"}),
#' or static (\code{"static"}).
#' @param as.df Logical scalar. When TRUE returns a data.frame.
#' @param no.print Logical scalar. When TRUE suppress screen messages.
#' @param skip.moran Logical scalar. When TRUE Moran's I is not reported (see details).
#' @param name Character scalar. Name of the diffusion network (descriptive).
#' @param behavior Character scalar. Name of the behavior been analyzed (innovation).
#' @export
#' @seealso Default options are listed at \code{\link{netdiffuseR-options}}
#' @details
#'
#' \code{diffnet} objects hold both, static and dynamic vertex attributes. When
#' creating \code{diffnet} objects, these can be specified using the arguments
#' \code{vertex.static.attrs} and \code{vertex.dyn.attrs}; depending on whether
#' the attributes to specify are static or dynamic, \pkg{netdiffuseR} currently
#' supports the following objects:
#'
#' \tabular{llr}{
#' \strong{Class}    \tab \strong{Dimension}               \tab \strong{Check sorting}\cr
#' \emph{Static attributes} \cr
#' \code{matrix}     \tab with \eqn{n} rows                \tab \code{id} \cr
#' \code{data.frame} \tab with \eqn{n} rows                \tab \code{id} \cr
#' \code{vector}     \tab of length \eqn{n}                \tab - \cr\cr
#' \emph{Dynamic attributes} \cr
#' \code{matrix}     \tab with \eqn{n\times T}{n * T} rows \tab \code{id}, \code{per} \cr
#' \code{data.frame} \tab with \eqn{n\times T}{n * T} rows \tab \code{id}, \code{per}\cr
#' \code{vector}     \tab of length \eqn{n\times T}{n*T}   \tab - \cr
#' \code{list}       \tab of length \eqn{T} with matrices or data.frames of \eqn{n} rows\tab \code{id}, \code{per}\cr
#' }
#'
#' The last column, \strong{Check sorting}, lists the variables that
#' the user should specify if he wants the function to check the order of the rows
#' of the attributes (notice that this is not possible for the case of vectors).
#' By providing the name of the vertex id variable, \code{id}, and the time period
#' id variable, \code{per}, the function makes sure that the attribute data is
#' presented in the right order. See the example below. If the user does not
#' provide the names of the vertex id and time period variables then the function
#' does not check the way the rows are sorted, further it assumes that the data
#' is in the correct order.
#'
#' Plotting is done via the function \code{\link[sna:gplot]{gplot}},
#' and its layout via \code{\link[sna:gplot.layout]{gplot.layout}}, both from
#' the (\pkg{sna}) package.
#'
#' \code{vertex.cex} can either be a numeric scalar, a numeric vector or a character
#' scalar taking any of the following values \code{"degree"}, \code{"indegree"}, or
#' \code{"outdegree"}. The later will be passed to \code{\link{dgr}} to calculate
#' degree of the selected slice and will be normalized as
#'
#' \deqn{vertex.cex = d/[max(d) - min(d)]\times 2 + .5}{vertex.cex = d/[max(d) - min(d)]* 2 + .5}
#'
#' where \code{d=sqrt(dgr(graph))}.
#'
#' In the case of the \code{summary} method, Moran's I is calculated over the
#' cumulative adoption matrix using as weighting matrix the inverse of the geodesic
#' distance matrix. All this via \code{\link{moran}}. For each time period \code{t},
#' this is calculated as:
#'
#' \preformatted{
#'  m = moran(C[,t], G^(-1))
#' }
#'
#' Where \code{C[,t]} is the t-th column of the cumulative adoption matrix,
#' \code{G^(-1)} is the element-wise inverse of the geodesic matrix at time \code{t},
#' and \code{moran} is \pkg{netdiffuseR}'s moran's I routine. When \code{skip.moran=TRUE}
#' Moran's I is not reported. This can be useful when the graph is particuarly
#' large (tens of thousands of vertices) as when doing so geodesic distances are
#' not calculated, which avoids allocating a square matrix of size \eqn{n} on
#' the memory. As a difference from the adjacency matrices, the matrix with the
#' geodesic distances can't be stored as a sparse matrix (saving space).
#'
#' @section Auxiliary functions:
#'
#' \code{diffnet.attrs} Allows retriving network attributes. In particular, by default
#' returns a list of length \eqn{T} with data frames with the following columns:
#'
#' \enumerate{
#'  \item \code{per} Indicating the time period to which the observation corresponds.
#'  \item \code{toa} Indicating the time of adoption of the vertex.
#'  \item Further columns depending on the vertex and graph attributes.
#' }
#'
#' Each vertex static attributes' are repeated \eqn{T} times in total so that these
#' can be binded (\code{rbind}) to dynamic attributes.
#'
#' When \code{as.df=TRUE}, this convenience function is useful as it can be used
#' to create event history (panel data) datasets used for model fitting.
#'
#' Conversely, the replacement method allows including new vertex or graph
#' attributes either dynamic or static (see examples below).
#'
#' \code{diffnet.toa(graph)} works as an alias of \code{graph$toa}.
#' The replacement method, \code{diffnet.toa<-} used as \code{diffnet.toa(graph)<-...},
#' is the right way of modifying times of adoption as when doing so it
#'  performs several checks on the time ranges, and
#' recalculates adoption and cumulative adoption matrices using \code{toa_mat}.
#'
#' \code{nodes(graph)} is an alias for \code{graph$meta$ids}.
#'
#'
#' @family diffnet methods
#' @family data management functions
#' @aliases diffnet diffnet-class
#' @examples
#'
#' # Creating a random graph
#' set.seed(123)
#' graph <- rgraph_ba(t=9)
#' graph <- lapply(1:5, function(x) graph)
#'
#' # Pretty TOA
#' names(graph) <- 2001L:2005L
#' toa <- sample(c(2001L:2005L,NA), 10, TRUE)
#'
#' # Creating diffnet object
#' diffnet <- as_diffnet(graph, toa)
#' diffnet
#' summary(diffnet)
#'
#' # Plotting slice 4
#' plot(diffnet, t=4)
#'
#' # ATTRIBUTES ----------------------------------------------------------------
#'
#' # Retrieving attributes
#' diffnet.attrs(diffnet, "vertex", "static")
#'
#' # Now as a data.frame (only static)
#' diffnet.attrs(diffnet, "vertex", "static", as.df = TRUE)
#'
#' # Now as a data.frame (all of them)
#' diffnet.attrs(diffnet, as.df = TRUE)
#'
#' # Unsorted data -------------------------------------------------------------
#' # Loading example data
#' data(fakesurveyDyn)
#'
#' # Creating a diffnet object
#' fs_diffnet <- survey_to_diffnet(
#'    fakesurveyDyn, "id", c("net1", "net2", "net3"), "toa", "group",
#'    timevar = "time", keep.isolates=TRUE, warn.coercion=FALSE)
#'
#' # Now, we extract the graph data and create a diffnet object from scratch
#' graph <- fs_diffnet$graph
#' ids <- fs_diffnet$meta$ids
#' graph <- Map(function(g) {
#'   dimnames(g) <- list(ids,ids)
#'   g
#'   }, g=graph)
#' attrs <- diffnet.attrs(fs_diffnet, as.df=TRUE)
#' toa   <- diffnet.toa(fs_diffnet)
#'
#' # Lets apply a different sorting to the data to see if it works
#' n <- nrow(attrs)
#' attrs <- attrs[order(runif(n)),]
#'
#' # Now, recreating the old diffnet object (notice -id.and.per.vars- arg)
#' fs_diffnet_new <- as_diffnet(graph, toa=toa, vertex.dyn.attrs=attrs,
#'    id.and.per.vars = c("id", "per"))
#'
#' # Now, retrieving attributes. The 'new one' will have more (repeated)
#' attrs_new <- diffnet.attrs(fs_diffnet_new, as.df=TRUE)
#' attrs_old <- diffnet.attrs(fs_diffnet, as.df=TRUE)
#'
#' # Comparing elements!
#' tocompare <- intersect(colnames(attrs_new), colnames(attrs_old))
#' all(attrs_new[,tocompare] == attrs_old[,tocompare], na.rm = TRUE) # TRUE!
#'
#' # diffnetLapply -------------------------------------------------------------
#'
#' data(medInnovationsDiffNet)
#' diffnetLapply(medInnovationsDiffNet, function(x, cumadopt, ...) {sum(cumadopt)})
#'
#' @return
#' A list of class \code{diffnet} with the following elements:
#' \item{graph}{A list of length \eqn{T}. Containing sparse square matrices of size \eqn{n}
#' and class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.}
#' \item{toa}{An integer vector of size \eqn{T} with times of adoption.}
#' \item{adopt, cumadopt}{Numeric matrices of size \eqn{n\times T}{n*T} as those returned
#' by \code{\link{toa_mat}}.}
#' \item{vertex.static.attrs}{If not NULL, a data frame with \eqn{n} rows with vertex static
#' attributes.}
#' \item{vertex.dyn.attrs}{A list of length \eqn{T} with data frames containing vertex attributes
#' throught time (dynamic).}
#' \item{graph.attrs}{A data frame with \eqn{T} rows.}
#' \item{meta}{A list of length 9 with the following elements:
#' \itemize{
#'  \item \code{type}: Character scalar equal to \code{"dynamic"}.
#'  \item \code{class}: Character scalar equal to \code{"list"}.
#'  \item \code{ids}: Character vector of size \eqn{n} with vertices' labels.
#'  \item \code{pers}: Integer vector of size \eqn{T}.
#'  \item \code{nper}: Integer scalar equal to \eqn{T}.
#'  \item \code{n}: Integer scalar equal to \eqn{n}.
#'  \item \code{self}: Logical scalar.
#'  \item \code{undirected}: Logical scalar.
#'  \item \code{multiple}: Logical scalar.
#'  \item \code{name}: Character scalar.
#'  \item \code{behavior}: Character scalar.
#' }
#' }
#' @author George G. Vega Yon
as_diffnet <- function(graph, toa, t0=min(toa, na.rm = TRUE), t1=max(toa, na.rm = TRUE),
                       vertex.dyn.attrs = NULL, vertex.static.attrs= NULL,
                       id.and.per.vars = NULL,
                       graph.attrs = NULL,
                       undirected=getOption("diffnet.undirected"),
                       self=getOption("diffnet.self"),
                       multiple=getOption("diffnet.multiple"),
                       name = "Diffusion Network",
                       behavior = "Unspecified") {

  # Step 0.0: Check if its diffnet! --------------------------------------------
  if (inherits(graph, "diffnet")) {
    message("Nothing to do, the graph is already of class \"diffnet\".")
    return(graph)
  }

  # Step 1.1: Check graph ------------------------------------------------------
  meta <- classify_graph(graph)
  if (meta$type=="static") stop("-graph- should be dynamic (see ?`netdiffuseR-graphs`).")

  # Step 1.2: Checking that lengths fit
  if (length(toa)!=meta$n) stop("-graph- and -toa- have different lengths (",
                                meta$n, " and ", length(toa), " respectively). ",
                                "-toa- should be of length n (number of vertices).")

  # Step 2.1: Checking class of TOA and coercing if necesary -------------------
  if (!inherits(toa, "integer")) {
    warning("Coercing -toa- into integer.")
    toa <- as.integer(toa)
  }

  # Step 2.2: Checking names of toa
  if (!length(names(toa)))
    names(toa) <- meta$ids

  # Step 3.1: Creating Time of adoption matrix ---------------------------------
  mat <- toa_mat(toa, labels = meta$ids, t0=t0, t1=t1)

  # Step 3.2: Verifying dimensions and fixing meta$pers
  tdiff <- meta$nper - ncol(mat[[1]])
  if (tdiff < 0)
    stop("Range of -toa- is bigger than the number of slices in -graph- (",
         ncol(mat[[1]]), " and ", length(graph) ," respectively). ",
         "There must be at least as many slices as range of toa.")
  else if (tdiff > 0)
    stop("Range of -toa- is smaller than the number of slices in -graph- (",
         ncol(mat[[1]]), " and ", length(graph) ," respectively). ",
         "Please provide lower and upper boundaries for the values in -toa- ",
         "using -t0- and -t- (see ?toa_mat).")

  meta$pers <- as.integer(colnames(mat$adopt))

  # Step 4.0: Checking the attributes ------------------------------------------

  # Vertex dyn attrs
  vertex.dyn.attrs <-
    check_as_diffnet_attrs(vertex.dyn.attrs, meta, TRUE, id.and.per.vars)

  # Vertex static attrs
  vertex.static.attrs <-
    check_as_diffnet_attrs(vertex.static.attrs, meta, FALSE, id.and.per.vars)

  # Vertex static attrs
  graph.attrs <-
    check_as_diffnet_attrs(graph.attrs, meta, FALSE, id.and.per.vars, "graph")


  # Step 4.1: Change the class (or set the names) of the graph -----------------
  if (meta$class=="array") {
    graph <- apply(graph, 3, methods::as, Class="dgCMatrix")
  } else { # Setting names (if not before)
    if (!length(names(graph))) names(graph) <- meta$pers
    else if (any(names(graph) != meta$pers)) names(graph) <- meta$pers
  }

  # Step 5: Compleating attributes and building the object and returning
  meta$self       <- self
  meta$undirected <- undirected
  meta$multiple   <- multiple
  meta$name       <- ifelse(!length(name), "", ifelse(is.na(name), "",
                                                      as.character(name)))
  meta$behavior   <- ifelse(!length(behavior), "", ifelse(is.na(behavior), "",
                                                          as.character(behavior)))

  # Removing dimnames
  graph                  <- Map(function(x) Matrix::unname(x), x=graph)
  dimnames(toa)          <- NULL
  dimnames(mat$adopt)    <- NULL
  dimnames(mat$cumadopt) <- NULL

  return(structure(list(
    graph = graph,
    toa   = toa,
    adopt = mat$adopt,
    cumadopt = mat$cumadopt,
    # Attributes
    vertex.static.attrs = vertex.static.attrs,
    vertex.dyn.attrs    = vertex.dyn.attrs,
    graph.attrs         = graph.attrs,
    meta = meta
  ), class="diffnet"))
}

#' @export
#' @rdname as_diffnet
diffnet.attrs <- function(graph, element=c("vertex","graph"), attr.class=c("dyn","static"),as.df=FALSE) {
  nper <- graph$meta$nper
  pers <- graph$meta$pers
  n    <- graph$meta$n

  # Only for diffnet objects
  if (!inherits(graph, "diffnet")) stopifnot_graph(graph)

  # Checking elements
  if (any(!(element %in% c("vertex", "graph"))))
    stop("-element- should only have 'vertex', and/or 'graph'.")

  # Checking classes
  if (any(!(attr.class %in% c("dyn", "static"))))
    stop("-attr.class- should only have 'dyn', and/or 'static'.")

  # Expanding graph static attr
  g.static <- NULL
  v.dyn    <- NULL
  v.static <- NULL
  if ("graph" %in% element) g.static <- graph$graph.attrs
  if ("vertex" %in% element) {
    if ("dyn"    %in% attr.class) v.dyn    <- graph$vertex.dyn.attrs
    if ("static" %in% attr.class) v.static <- graph$vertex.static.attrs
  }

  # Parsing attributes
  if (!length(g.static)) g.static <- as.data.frame(matrix(ncol=0, nrow=n))
  if (!length(v.static)) v.static <- as.data.frame(matrix(ncol=0, nrow=n))
  if (!length(v.dyn[[1]])) v.dyn <- lapply(1:nper, function(y) as.data.frame(matrix(ncol=0, nrow=n)))

  attrs <- cbind(toa=graph$toa, v.static)
  out <- lapply(1:nper, function(y) {
    cbind(per=rep(pers[y], n),attrs, v.dyn[[y]])
  })

  if (as.df) {
    out <- do.call(rbind, out)
    return(data.frame(out, id=rep(graph$meta$ids, nper), row.names = NULL))
  }

  names(out) <- graph$meta$pers
  out
}

#' @rdname as_diffnet
#' @export
`diffnet.attrs<-` <- function(graph, element="vertex", attr.class="static", value) {
  .Defunct("[[<-.diffnet")
}

#' @rdname as_diffnet
#' @export
diffnet.toa <- function(graph) {
  if (!inherits(graph, "diffnet")) stop("-graph- must be a 'diffnet' object")
  graph$toa
}

#' @rdname as_diffnet
#' @export
`diffnet.toa<-` <- function(graph, i, value) {
  if (!inherits(graph, "diffnet")) stop("-graph- must be a 'diffnet' object")
  if (missing(i)) i <- graph$meta$ids

  # Checking values of the data: normalizing
  test <- !(value %in% c(graph$meta$pers, NA))
  if (any(test)) stop("Some elements of -value- (",
                      paste0(head(value[test], 20), collapse=", "),
                      ifelse(length(value[test]) > 20,", ...", "")
                      ,") are not within the range of the original graph.")

  # Changing the value of toa
  names(graph$toa) <- graph$meta$ids
  graph$toa[i] <- value
  graph$toa <- unname(graph$toa)

  # Recalculating adopt and cumadopt
  mat <- toa_mat(graph$toa, t0=graph$meta$pers[1], t1=graph$meta$pers[graph$meta$nper])

  # checking stack
  nper <- ncol(mat[[1]])
  graph$adopt    <- mat$adopt
  graph$cumadopt <- mat$cumadopt

  graph

}

