#' Read foreign graph formats
#'
#' Reading pajek and Ucinet files, this function returns weighted edgelists in the form of
#' data frames including a data frame of the vertices. (function on development)
#'
#' @param x Character scalar. Path to the file to be imported.
#' @return In the case of \code{read_pajek}, a list with three elements
#' \item{vertices}{A data frame with \eqn{n} rows and two columns: id and label}
#' \item{edges}{If not null, a list of data frames with three columns: ego, alter, w (weight)}
#' \item{edgelist}{If not null, a list of data frame with three columns: ego, alter, w (weight)}
#'
#' For \code{read_ml}, a list with two elements:
#' \item{adjmat}{An array with the graph}
#' \item{meta}{A list with metadata}
#'
#' @details
#' Since .net files allow working with multi-relational networks (more than one
#' class of edge), the function returns lists of edges and edgeslist with the corresponding
#' tag on the .net file. For example, if the .net file contains
#'
#' \preformatted{
#'  *Arcslist :9 "SAMPPR"
#'  ...
#'  *Arcslist :10 "SAMNPR"
#' }
#'
#' The output will include data frames of edgelists with those tags.
#'
#' @examples
#' # From .net: Sampson monastery data from UCINET dataset ---------------------
#'
#' # Reading the arcs/edges format
#' path <- system.file("extdata", "SAMPSON.NET", package = "netdiffuseR")
#' SAMPSON <- read_pajek(path)
#'
#' # Reading the arcslist/edgelist format
#' path <- system.file("extdata", "SAMPSONL.NET", package = "netdiffuseR")
#' SAMPSONL <- read_pajek(path)
#'
#' # From DL (UCINET): Sampson monastery data (again) --------------------------
#' path <- system.file("extdata", "SAMPSON.DAT", package = "netdiffuseR")
#' SAMPSONL <- read_ml(path)
#'
#' @author George G. Vega Yon
#' @export
#' @aliases read_net read_dl
#' @family graph formats
#' @source From the pajek manual \url{http://mrvar.fdv.uni-lj.si/pajek/pajekman.pdf}
read_pajek <- function(x) {
  # (1) Reading the file and finding the tags ----------------------------------
  lines <- readLines(x)

  # Removing empty lines (white)
  lines <- lines[!grepl("^\\s*$",lines)]

  nlines <- length(lines)
  tags <- which(grepl("^[*]",lines))
  tagnames <- tolower(gsub("^[*]", "", lines[tags]))

  # Finding start/end of the tags
  tags <- data.frame(start=tags+1, end=c(tags[-1] - 1, nlines), lab=tags)
  rownames(tags) <- tagnames

  # (2) Checking each tag ------------------------------------------------------

  # Reading vertices
  test <- which(grepl("^vertices", rownames(tags)))
  if (length(test)) {
    n <- as.integer(gsub("[a-zA-Z* ]", "", lines[tags$lab[test[1]]]))
    vertices <- lines[with(tags[test[1],,drop=FALSE], start:end)]

    # Parsing
    vertices <- as.data.frame(do.call(rbind, sapply(vertices, strsplit, split="\\s+")))
    vertices[,2] <- gsub("\"", "", vertices[,2])

    rownames(vertices) <- NULL
    colnames(vertices) <- c("id", "label")

  } else vertices <- NULL

  # Reading edges
  test <- which(grepl("^(arcslist|edgeslist)", rownames(tags)))
  edgelist <- NULL
  if (length(test)) {
    # Creating empty object
    edgelist <- vector("list", length(test))
    names(edgelist) <- lines[tags$lab[test]]

    for (i in test) {
      # Getting the lines of the class of arc
      subarc <- lines[tags$start[i]:tags$end[i]]

      # Coercing into edgelist
      subarc <- do.call(rbind, lapply(subarc, function(x) {
        x <- as.numeric(strsplit(x, split="\\s+")[[1]])
        if (length(x)-1) data.frame(ego=x[1], alter=x[-1], w=1)
        else data.frame(ego=x, alter=x, w=1)
      }))

      # Including it into the edgelist
      edgelist[[ lines[tags$lab[i]] ]] <- subarc
    }

  }

  # Reading edgelist
  test <- which(grepl("^(arcs|edges)\\s+", rownames(tags)))
  edges <- NULL
  if (length(test)) {

    # Creating empty object
    edges <- vector("list", length(test))
    names(edges) <- lines[tags$lab[test]]

    for (i in test) {
      # Getting the lines of the class of arc
      subedge <- lines[tags$star[i]:tags$end[i]]

      # Coercing into edgelist
      subedge <- do.call(rbind, lapply(subedge, function(x) {
        x <- as.numeric(strsplit(x, split="\\s+")[[1]])

        if (length(x) > 2) {
          data.frame(ego=x[1], alter=x[2], w=x[3])
        } else data.frame(ego=x[1], alter=x[2], w=1)
      }))

      # Including it into the edgelist
      edges[[ lines[tags$lab[i]] ]] <- subedge
    }
  }

#   # Reading Matrices
#   test <- which(grepl("^matrix", names(tags)))
#   matrices <- NULL
#   if (length(test)) {
#     for (i in test)
#       matrices <- c(matrices, lines[tags[i]])
#   }

  return(list(vertices=vertices, edges=edges, edgelist=edgelist))
}


# regmatches
gregexec <- function(s, pattern) {
  lapply(regmatches(s, gregexpr(pattern, s)),
         function(e) regmatches(e, regexec(pattern, e)))
}

# gregexec("dl n=4 format=fullmatrix", "n[= ,]+([0-9]+)")[[1]]
# gregexec("dl n,4 format=fullmatrix", "n[= ,]+([0-9]+)")[[1]]
# gregexec("dl n 4 format=fullmatrix", "n[= ,]+([0-9]+)")[[1]]
# http://mrvar.fdv.uni-lj.si/sola/info4/multinet/multinet.htm
# x <- c("n=4 patt:abc ajajaj c:76 as")


#' Read UCINET graph files
#' Other datasets http://moreno.ss.uci.edu/data.html
#' @rdname read_pajek
#' @export
read_ml <- function(x) {
  # (1) Reading the file and finding the tags ----------------------------------
  lines <- readLines(x)

  # Collapsing/splitting the first lines: In order to read DL file easily, we
  #   will assume that data labels in the form of [a-zA-Z ]+: start immediatly
  #   from the second element in the vector. If not we will force it to do it so.
  #
  #   To work this around, we also have a function that splits an element into two
  #   or more lines until all [a...]: start in different lines.
  lines <- unlist(strsplit(lines, split="\\s+(?=\\s+[a-zA-Z]+[:])", perl = TRUE))
  lines <- unlist(strsplit(lines, split="(?<=[:])", perl = TRUE))

  # Correcting badsplits (COLUMN LAB, ROW LAB, etc.)
  tocheck <- which(
    grepl("^labels[:]$", lines[-1], ignore.case = TRUE) &
      grepl("(COLUMN|COL|ROW|LEVEL)$", lines[-length(lines)], ignore.case = TRUE)
    )

  if (length(tocheck)) {
    lines[tocheck+1] <- paste(lines[tocheck], lines[tocheck + 1])
    lines[tocheck]   <- gsub("(COLUMN|COL|ROW|LEVEL)$", "", lines[tocheck], ignore.case = TRUE)
  }

  # Identifying which should be merged
  start <- which(grepl("[a-zA-Z ]+[:]", lines))[1]

  if (start == 1) { # split
    lines <- c(
      gsub("[a-zA-Z ]+[:].+", "", lines[1]),
      gsub(".+(?=[a-zA-Z ]+[:])", "" ,lines[1]),
      lines[-1])
  } else if (start > 2) { # Merging
    subsec <- 1:(start-1)
    lines <- c(
      paste0(lines[subsec], collapse = " "),
      lines[-subsec]
      )
  }

  # Trimming data
  lines <- trimws(lines, "both")

  # Taking the first row and all the : into lowercase
  lines[1] <- tolower(lines[1])
  test <- which(grepl("^[a-zA-Z ]+[:]$", lines))
  lines[test] <- tolower(lines[test])

  nlines <- length(lines)
  n      <- as.integer(gregexec(lines[1], "n[= ,]+([0-9]+)")[[1]][[1]][2])

  # By default format is fullmatrix
  if (grepl("format[= ,]+[a-z]+", lines[1])) format <- gregexec(lines[1], "format[= ,]+([a-z]+)")[[1]][[1]][2]
  else format <- "fullmatrix"

  # By default the number of matrices is 1
  if (grepl("nm[= ,]+[0-9]+", lines[1])) nm <- as.integer(gregexec(lines[1], "nm[= ,]+([0-9]+)")[[1]][[1]][2])
  else nm <- 1L

  # Others: row labels embedded, col labels embedded, labels embedded, diagonal = presen
  if (grepl("row\\s+labels\\s+embedded", lines[1])) row_labels_embedded <- TRUE
  else row_labels_embedded <- FALSE

  if (grepl("col\\s+labels\\s+embedded", lines[1])) col_labels_embedded <- TRUE
  else col_labels_embedded <- FALSE

  if (grepl("(?<!(col|row))\\s+labels\\s+embedded", lines[1], perl = TRUE)) col_labels_embedded <- TRUE
  else col_labels_embedded <- FALSE

  if (grepl("diagonal[= ,]+[a-z]+", lines[1])) diagonal <- gregexec(lines[1], "diagonal[= ,]+([a-z]+)")[[1]][[1]][2]
  else diagonal <- "present"

  # Creating meta
  meta <- list(n=n, format=format, nm=nm, row_labels_embedded=row_labels_embedded,
               col_labels_embedded=col_labels_embedded, diagonal=diagonal)

  # ----------------------------------------------------------------------------
  # Finding start/end of the tags
  tags <- which(grepl("[a-zA-Z]+:", lines))
  tagnames <- tolower(gsub("[:]$", "", lines[tags]))
  tags <- data.frame(start=tags+1, end=c(tags[-1] - 1, nlines), lab=tags)
  rownames(tags) <- tagnames

  # ----------------------------------------------------------------------------
  # Reading the data
  dat <- tags["data",,drop=FALSE]
  dat <- with(dat, lines[start:end])
  dat <- as.numeric(unlist(strsplit(dat, "\\s+")))

  # R reads in column fashion, so we need to transpose later
  dat <- array(dat, dim=c(n,n,nm))
  dat <- array(apply(dat, 3, t), dim=c(n,n,nm))

  # ----------------------------------------------------------------------------
  # Checking colnames
  test <- which(grepl("col(umn)?\\s+labels", rownames(tags)))
  if (length(test)) {
    cnames <- tags[test,,drop=FALSE]
    cnames <- with(cnames, lines[start:end])
    cnames <- strsplit(cnames, split = ",")
    colnames(dat) <- cnames
  } else colnames(dat) <- 1:ncol(dat)

  # Checking Rownames
  test <- which(grepl("row\\s+labels", rownames(tags)))
  if (length(test)) {
    rnames <- tags[test,,drop=FALSE]
    rnames <- with(rnames, lines[start:end])
    rnames <- strsplit(rnames, split = ",")
    rownames(dat) <- rnames
  } else rownames(dat) <- 1:row(dat)

  # Checking levels
  test <- which(grepl("level\\s+labels", rownames(tags)))
  if (length(test)) {
    lnames <- tags[test,,drop=FALSE]
    lnames <- with(lnames, lines[start:end])
    lnames <- strsplit(lnames, split = ",")
    dimnames(dat)[[3]] <- lnames
  } else dimnames(dat)[[3]] <- 1:dim(dat)[3]

  # ----------------------------------------------------------------------------
  # Out for data
  return(
    list(adjmat=dat, meta=meta)
    )
}


UCINET_datatype <- c(
  "nodt","bytedt","booleandt","shortintdt","worddt","smallintdt","longintdt",
  "singledt","realdt","doubledt","compdt","extendeddt","labeldt","setdt",
  "stringdt","pointerdt","chardt","integerdt","nodelistdt","sparsedt","int64dt"
)
# UCINET_datatype <- data.frame(
#   des   = UCINET_datatype,
#   bytes = c(0, 1, 1, 2)
# )

#' Reads UCINET files
#' @param f Character scalar. Name of the header file. e.g. \code{mydata.##h}.
#' @export
#' @return An array including dimnames (if there are) and the following attributes:
#' \item{headerversion}{Character scalar}
#' \item{year}{Integer. Year the file was created}
#' \item{month}{Integer. Month of the year the file was created.}
#' \item{day}{Integer. Day of the month the file was created.}
#' \item{dow}{Integer. Day of the week the file was created.}
#' \item{labtype}{}
#' \item{infile.dt}{Character scalar. Type of data of the array.}
#' \item{dim}{Integer vector. Dimmensions of the array.}
#' \item{tit}{Character scalar. Title of the file.}
#' \item{haslab}{Logical vector. Whether  each dim has a label.}
#' @aliases ucinet UCINET
#' @family graph formats
read_ucinet_head <- function(f) {
  con <- file(f,"rb")
  on.exit(close(con))
  readBin(con, raw(), n=2L)

  # Common meta
  meta <- list(
    headerversion = rawToChar(readBin(con, raw(), n=4L)),
    year          = readBin(con, integer(), size=2L, signed = FALSE),
    month         =readBin(con, integer(), size=2L, signed = FALSE),
    day           = readBin(con, integer(), size=2L, signed = FALSE),
    dow           = readBin(con, integer(), size=2L, signed = FALSE),
    labtype       = readBin(con, integer(), size=2L, signed = FALSE),
    infile.dt     = UCINET_datatype[readBin(con, integer(), size=1L, signed = FALSE)]
  )

  # Dims
  ndim <- readBin(con, integer(), size=2L, signed = FALSE)
  meta$dim <- readBin(con, integer(), size=4L, n=ndim)

  # Title
  ntit      <- readBin(con, integer(), size=1L, signed = FALSE)
  meta$tit  <- readChar(con, nchars = ntit)

  # Labels
  meta$haslab <- readBin(con, logical(), size = 1L, n=3L)

  # Col labels
  if (meta$haslab[1]) {
    meta$clab <- vector("character", meta$dim[1])
    for (i in 1:meta$dim[1]) {
      lablen       <- readBin(con, integer(), size=2L)
      meta$clab[i] <- paste(readBin(con, character(), n=lablen/2, size=1L), collapse="")
    }
  } else {
    meta$clab <- NULL
  }

  # Row labels
  if (meta$haslab[2]) {
    meta$rlab <- vector("character", meta$dim[2])
    for (i in 1:meta$dim[2]) {
      lablen       <- readBin(con, integer(), size=2L)
      meta$rlab[i] <- paste(readBin(con, character(), n=lablen/2, size=1L), collapse="")
    }
  } else {
    meta$rlab <- NULL
  }

  # Level labels
  if (meta$haslab[3]) {
    meta$llab <- vector("character", meta$dim[3])
    for (i in 1:meta$dim[3]) {
      lablen       <- readBin(con, integer(), size=2L)
      meta$llab[i] <- paste(readBin(con, character(), n=lablen/2, size=1L), collapse="")
    }
  } else {
    meta$llab <- NULL
  }


  return(meta)
}

#' Read UCINET files (binary)
#' @param echo Logical scalar. When \code{TRUE} shows a message.
#' @export
#' @rdname read_ucinet_head
read_ucinet <- function(f, echo=FALSE) {

  # Recursion
  if (length(f) > 1)
    return(lapply(f, read_ucinet, echo=echo))

  # Checking extension
  f <- gsub("\\.[#]{2}(h|d)$", "", f)
  if (echo) message(sprintf("Reading file %s", f))

  # File names (do they exits)
  fhead <- paste0(f,".##h")
  fdata <- paste0(f,".##d")

  if (!file.exists(fhead))
    stop(sprintf("File %s not found"))

  if (!file.exists(fdata))
    stop(sprintf("File %s not found"))

  # Getting Metadata
  meta <- tryCatch(read_ucinet_head(fhead), error=function(e) e)
  if (inherits(meta, "error")) {
    stop(sprintf("An error ocurred while processing the header %s\n", fhead),
         meta)
  }

  # Reading data and coercing into an array
  con <- file(fdata, "rb")
  on.exit(close(con))

  s <- prod(meta$dim)
  dat <- readBin(con, numeric(), s, size=4L)
  dat <- do.call(structure, c(list(.Data=dat), meta))

  # Transposing (since the data is stored by row in the .##d file)
  for (i in 1:dim(dat)[3])
    dat[,,i] <- t(dat[,,i])

  dimnames(dat) <- list(meta$rlab, meta$clab, meta$llab)

  return(dat)
}
