# Auxiliar function to warn and coerce classes
check_var_class_and_coerce <- function(var, dat, class.ok, class.target, warn.coercion) {

  # Has the right class?
  cl <- class(dat[[var]])
  if (!any(cl %in% class.ok))
    stop("The variable -",var,"- is of class ", cl, " which is ",
         "not supported. Supported class(es): ", paste0(class.ok, collapse=", "),
         ".")

  # Should be coherced?
  if (warn.coercion && !(class.target %in% cl))
    warning("Coercing -",var, "- into ",class.target,".")

  # Returning
  return(methods::as(dat[[var]], class.target))
}

#' Convert survey-like data and edgelists to a \code{diffnet} object
#'
#' These convenient functions turn network nomination datasets and edgelists with
#' vertex attributes datasets into diffnet objects. Both work as wrappers of
#' \code{\link{edgelist_to_adjmat}} and \code{\link{as_diffnet}}.
#'
#' @inheritParams edgelist_to_adjmat
#' @param dat A data frame.
#' @param idvar Character scalar. Name of the id variable.
#' @param netvars Character vector. Names of the network nomination variables.
#' @param toavar Character scalar. Name of the time of adoption variable.
#' @param timevar Character sacalar. In the case of longitudinal data, name of the time var.
#' @param groupvar Character scalar. Name of cohort variable (e.g. city).
#' @param no.unsurveyed Logical scalar. When \code{TRUE} the nominated individuals
#' that do not show in \code{idvar} are set to \code{NA} (see details).
#' @param warn.coercion Logical scalar. When \code{TRUE} warns coercion from numeric to integer.
#' @param ... Further arguments to be passed to \code{\link{as_diffnet}}.
#' @details
#'
#' All of \code{netvars}, \code{toavar} and \code{groupvar}
#' must be integers. Were these numeric they are coerced into integers, otherwise,
#' when neither of both, the function returns with error. \code{idvar}, on the
#' other hand, should only be integer when calling \code{survey_to_diffnet},
#' on the contrary, for \code{edgelist_to_diffnet}, \code{idvar} may be character.
#'
#' In field work it is not unusual that some respondents nominate unsurveyed
#' individuals. In such case, in order to exclude them from the analysis,
#' the user can set \code{no.unsurveyed=TRUE} (the default), telling the
#' function to exclude such individuals from the adjacency matrix. This is
#' done by setting variables in \code{netvars} equal to \code{NA} when the
#' nominated id can't be found in \code{idvar}.
#'
#' If the network nomination process was done in different groups (location
#' for example) the survey id numbers may be define uniquely within each group
#' but not across groups (there may be many individuals with \code{id=1},
#' for example). To encompass this issue, the user can tell the function what
#' variable can be used to distinguish between groups through the \code{groupvar}
#' argument. When \code{groupvar} is provided, function redifines \code{idvar}
#' and the variables in \code{netvars} as follows:
#'
#' \preformatted{
#'    dat[[idvar]] <- dat[[idvar]] + dat[[groupvar]]*z
#' }
#'
#' Where \code{z = 10^nchar(max(dat[[idvar]]))}.
#'
#' For longitudinal data, it is assumed that the \code{toavar} holds the same
#' information through time, this is, time-invariable. This as the package does
#' not yet support variable times of adoption.
#'
#' The \code{fill.missing} option can take any of these three values: \code{"edgelist"},
#' \code{"dat"}, or \code{"both"}. This argument works as follows:
#' \enumerate{
#'    \item When \code{fill.missing="edgelist"} (or \code{"both"}) the function
#'    will check which vertices show in \code{dat} but do not show in \code{edgelist}.
#'    If there is any, the function will include these in \code{edgelist} as ego to
#'    \code{NA} (so they have no link to anyone), and, if specified, will fill
#'    the \code{t0}, \code{t1} vectors with \code{NA}s for those cases. If
#'    \code{w} is also specified, the new vertices will be set to
#'    \code{min(w, na.rm=TRUE)}.
#'    \item When \code{fill.missing="dat"} (or \code{"both"}) the function
#'    checks which vertices show in \code{edgelist} but not in \code{dat}. If
#'    there is any, the function will include these in \code{dat} by adding
#'    one row per individual.
#' }
#'
#' @export
#' @return A \code{\link{diffnet}} object.
#' @seealso \code{\link{fakesurvey}}, \code{\link{fakesurveyDyn}}
#' @family data management functions
#' @author
#' Vega Yon
#' @examples
#' # Loading a fake survey (data frame)
#' data(fakesurvey)
#'
#' # Diffnet object keeping isolated vertices ----------------------------------
#' dn1 <- survey_to_diffnet(fakesurvey, "id", c("net1", "net2", "net3"), "toa",
#'    "group", keep.isolates=TRUE)
#'
#' # Diffnet object NOT keeping isolated vertices
#' dn2 <- survey_to_diffnet(fakesurvey, "id", c("net1", "net2", "net3"), "toa",
#'    "group", keep.isolates=FALSE)
#'
#' # dn1 has an extra vertex than dn2
#' dn1
#' dn2
#'
#' # Loading a longitudinal survey data (two waves) ----------------------------
#' data(fakesurveyDyn)
#'
#' groupvar <- "group"
#' x <- survey_to_diffnet(
#'    fakesurveyDyn, "id", c("net1", "net2", "net3"), "toa", "group" ,
#'    timevar = "time", keep.isolates = TRUE, warn.coercion=FALSE)
#'
#' plot_diffnet(x, vertex.cex = 1.5, displaylabels = TRUE)
#'
#' # Reproducing medInnovationsDiffNet object ----------------------------------
#' data(medInnovations)
#'
#' # What are the netvars
#' netvars <- names(medInnovations)[grepl("^net", names(medInnovations))]
#'
#' medInnovationsDiffNet2 <- survey_to_diffnet(
#'    medInnovations,
#'    "id", netvars, "toa", "city",
#'    warn.coercion=FALSE)
#'
#' medInnovationsDiffNet2
#'
#' # Comparing with the package's version
#' all(diffnet.toa(medInnovationsDiffNet2) == diffnet.toa(medInnovationsDiffNet)) #TRUE
#' all(
#'    diffnet.attrs(medInnovationsDiffNet2, as.df = TRUE) ==
#'    diffnet.attrs(medInnovationsDiffNet, as.df = TRUE),
#'    na.rm=TRUE) #TRUE
#'
#'
survey_to_diffnet <- function(
  dat, idvar, netvars, toavar,
  groupvar=NULL,
  no.unsurveyed=TRUE,
  timevar=NULL,
  t = NULL,
  undirected = getOption("diffnet.undirected", FALSE),
  self = getOption("diffnet.self", FALSE),
  multiple=getOption("diffnet.multiple", FALSE),
  keep.isolates=TRUE, recode.ids=TRUE,
  warn.coercion=TRUE,
  ...) {

  # Creating a varlist
  varlist <- c(idvar, groupvar, netvars, toavar, timevar)

  # Are all in the dataset??
  test <- varlist %in% colnames(dat)
  if (any(!test))
    stop("Variables -", paste(varlist[!test], collapse = "-, -"),"- can't be found on -dat-.")

  # Coercing data into numeric variables
  for (x in varlist) {
    dat[[x]] <- check_var_class_and_coerce(
      x, dat, c("numeric", "integer"), "integer", warn.coercion)
  }

  # Changing ids
  if (length(groupvar)) {
    idord <- nchar(max(dat[[idvar]]))
    for (x in c(netvars, idvar))
      dat[[x]] <- dat[[x]] + dat[[groupvar]]*(10^(idord))
  }

  # Removing unsurveyed
  if (no.unsurveyed) {
    surveyed <- unique(dat[[idvar]])
    for (x in netvars)
      dat[[x]][which(!(dat[[x]] %in% surveyed))] <- NA
  }

  # Analyzing time data
  if (length(timevar)) {
    # Checking if data is complete
    test <- complete.cases(dat[[timevar]])
    if (any(!test)) {
      test <- which(!test)
      test <- paste0(test, collapse=", ")
      stop("Some elements of -timevar- have missing data:\n\t",
           ifelse(nchar(test) > 80, paste(strtrim(test,80), "..."), test),".")
    }

    tvar <- dat[[timevar]]
    tran <- range(tvar, na.rm=TRUE)
    tran <- tran[1]:tran[2]
  } else {
    tvar <- rep(1, nrow(dat))
    tran    <- 1
  }

  # Reshaping data (so we have an edgelist)
  dat.long     <- NULL
  t0           <- NULL
  vertex.attrs <- vector("list", length(tran))
  colstoexport <- which(!(colnames(dat) %in% c(toavar)))
  for (i in 1:length(tran)) {
    subs            <- dat[tvar == tran[i],]
    vertex.attrs[[i]] <- subs[,colstoexport]

    # Reshaping
    subs <- reshape(
      subs[,c(idvar, netvars)], v.names= "net",
      varying = netvars,
      idvar="id", direction="long")[,c(idvar, "net")]

    # Creating edgelist
    dat.long <- rbind(dat.long, subs)

    # Times for dyn networks
    t0 <- c(t0, rep(tran[i], nrow(subs)))
  }
  t1 <- t0

  # If the time range equals 1, then it implies that the graph data is static
  if (length(tran) == 1) {
    # Computing the times of adoption
    rtoa <- range(dat[[toavar]], na.rm = TRUE)
    t    <- rtoa[2] - rtoa[1] + 1

    # Creating the adjacency matrix
    graph <- edgelist_to_adjmat(edgelist = dat.long, t = t,
                         undirected=undirected, self=self, multiple = multiple,
                         keep.isolates = keep.isolates, recode.ids = recode.ids)

  } else {
    # Creating the adjacency matrix
    graph <- edgelist_to_adjmat(edgelist = dat.long, t0 = t0, t1=t1,
                         undirected=undirected, self=self, multiple = multiple,
                         keep.isolates = keep.isolates, recode.ids = recode.ids)

  }

  # Used vertices
  used.vertex <- data.frame(rownames(graph[[1]]),
                            `_original_sort` = seq_len(nnodes(graph)),
                            check.names = FALSE)

  colnames(used.vertex)[1] <- idvar
  for (i in 1:length(tran)) {
    vertex.attrs[[i]] <- merge(
      used.vertex, vertex.attrs[[i]],
      by = idvar,
      all.x=TRUE, sort=FALSE)

    # Sorting rows
    vertex.attrs[[i]] <- vertex.attrs[[i]][
      order(vertex.attrs[[i]][["_original_sort"]]),]

    # Removing the idvar
    test <- colnames(vertex.attrs[[i]]) %in% c(idvar, "_original_sort")
    vertex.attrs[[i]] <- vertex.attrs[[i]][,which(!test)]
  }

  # Times of adoption
  dat <- unique(dat[,c(idvar, toavar)])
  toa <- merge(used.vertex, dat, by=idvar, all.x=TRUE, sort=FALSE)

  # Sorting rows
  toa <- toa[order(toa[["_original_sort"]]),][[toavar]]

  if (length(toa) != nrow(used.vertex))
    stop("It seems that -toavar- is not time-invariant.")

  if (length(tran) == 1) {
    as_diffnet(
      graph=graph, toa=toa,
      vertex.static.attrs = vertex.attrs[[1]],
      ...
    )
  } else {
    as_diffnet(
      graph=graph, toa=toa,
      vertex.dyn.attrs = vertex.attrs,
      ...
    )
  }
}

#' @rdname survey_to_diffnet
#' @param fill.missing Character scalar. In the case of having unmatching ids
#' between \code{dat} and \code{edgelist}, fills the data (see details).
#' @export
edgelist_to_diffnet <- function(edgelist, w=NULL,
                                t0=NULL, t1=NULL ,
                                dat, idvar, toavar, timevar=NULL,
                                undirected = getOption("diffnet.undirected", FALSE),
                                self = getOption("diffnet.self", FALSE),
                                multiple=getOption("diffnet.multiple", FALSE),
                                fill.missing=NULL,
                                keep.isolates=TRUE, recode.ids=TRUE,
                                warn.coercion=TRUE) {

  # Step 0.1: Checking dat -----------------------------------------------------
  # Creating a varlist
  varlist <- c(idvar, toavar, timevar)

  # Are all in the dataset??
  test <- varlist %in% colnames(dat)
  if (any(!test))
    stop("Variables -", paste(varlist[!test], collapse = "-, -"),"- can't be found on -dat-.")

  # Is it complete? toa may be empty
  test <- which(!complete.cases(dat[,varlist[-2]]))
  if (length(test))
    stop("Incomplete cases in -dat-. All observations must have -idvar- and,",
         " if specified, -timevar-. The following rows are incomplete:\n\t",
         paste0(test, collapse=", "), ".")

  # Coercing data into numeric variables. idvar can be names
  for (x in varlist[-1])
    dat[[x]] <- check_var_class_and_coerce(
      x, dat, c("numeric", "integer"), "integer", warn.coercion)

  # Step 0.2: Checking FILL data -----------------------------------------------
  ids.edgelist <- unique(c(edgelist[,1,drop=TRUE], edgelist[,2,drop=TRUE]))
  ids.dat      <- unique(dat[[idvar]])
  if (length(fill.missing)) {
    # Checking argument
    if (!inherits(fill.missing, "character")) {
      stop("-fill.missing- should be a character scalar")
    }
    # If something goes wrong
    if (!(fill.missing %in% c("edgelist", "dat", "both"))) {
      stop("The only values currently supported for -fill.missing- are:\n\t",
           "'edgelist', 'dat', or 'both'.")
    }

    # Filling missing pieces
    if (fill.missing %in% c("edgelist", "both")) {
      test <- ids.dat[which(!( ids.dat %in% ids.edgelist))]

      # If some missing, then filling with more edges
      if (length(test)) {
        warning("The following ids will be added to -edgelist-:\n\t",
                paste0(test, collapse=", "),".")

        nedgelist <- nrow(edgelist)
        edgelist  <- rbind(
          edgelist,
          edgelist[1:length(test),, drop=FALSE]
          )

        edgelist[(nedgelist + 1):nrow(edgelist),]   <- NA
        edgelist[(nedgelist + 1):nrow(edgelist),1]  <- test

        if (length(t0)) t0 <- c(t0, rep(NA, length(test)))
        if (length(t1)) t1 <- c(t1, rep(NA, length(test)))
        if (length(w))   w <- c(w, rep(min(w, na.rm = TRUE), length(test)))
      }
    }
    if (fill.missing %in% c("dat", "both")) {
      test <- ids.edgelist[which(!(ids.edgelist %in% ids.dat))]

      # If some missing, then filling with more edges
      if (length(test)) {
        warning("The following ids will be added to -dat-:\n\t",
                paste0(test, collapse=", "),".")

        ndat <- nrow(dat)
        dat  <- rbind(
          dat,
          dat[1:length(test),, drop=FALSE])

        dat[(ndat + 1):nrow(dat),] <- NA
        dat[[idvar]][(ndat + 1):nrow(dat)]  <- test
      }
    }
  }

  # Step 1.1: Converting edgelist into adjmat ------------------------------------
  adjmat <- edgelist_to_adjmat(
    edgelist, w=w, t0=t0, t1=t1,
    undirected = undirected, self=self, multiple=multiple,
    keep.isolates = keep.isolates,
    recode.ids = recode.ids, simplify = FALSE)

  # Step 1.2: Checking times in edgelist and in dat (if any) -------------------
  suppressWarnings(dat.ran.toavar <- range(dat[[toavar]], na.rm = TRUE))

  # If the range turns out to be infinite, then error
  if (any(!is.finite(dat.ran.toavar)))
    stop("Invalid Times of Adoption. When computing its is undefine.")

  if (length(timevar)) {
    dat.ran.timevar <- range(dat[[timevar]])

    # range(toa) %within% range(timevar)
    if (dat.ran.toavar[1] < dat.ran.timevar[1] ||
        dat.ran.toavar[2] > dat.ran.timevar[2])
      stop("Invalid range in -toavar- (",dat.ran.toavar[1], " to ",
           dat.ran.toavar[2],"). It should be within the range of -timevar-",
           " (",dat.ran.timevar[1]," to ",dat.ran.timevar[2],").")

    # Setting the range from the data to be the timevar
    dat.ran <- dat.ran.timevar

  } else dat.ran <- dat.ran.toavar

  # Auxiliary time range
  tran <- dat.ran[1]:dat.ran[2]

  # Number of observations in adjmat
  if (length(t0) | length(t1)) { # If dynamic, we have to check
    edge.ran <- as.integer(names(adjmat))
    edge.ran <- c(edge.ran[1], edge.ran[length(adjmat)])

    # Range of dat and edgelist should be equal
    if (any(dat.ran != edge.ran))
      stop("Time ranges in -edgelist- and -dat- should be the same. Currently ",
           "they are ",paste0(edge.ran, collapse = " to "), " and ",
           paste0(dat.ran, collapse=" to "), " respectively.")

  } else { # If no dynamic, then simply replicate it
    adjmat        <- lapply(tran, function(x) adjmat[[1]])
    names(adjmat) <- tran
  }

  # Step 2: Getting the ids and checking everything is in order ----------------
  used.vertex           <- data.frame(rownames(adjmat[[1]]),
                                      `_original_sort` = seq_len(nnodes(adjmat)),
                                      check.names = FALSE)
  colnames(used.vertex)[1] <- idvar

  # All in the edgelist?
  dat.idvar <- unique(dat[[idvar]])
  test  <- which(!(dat.idvar %in% used.vertex[[idvar]]))
  if (length(test))
    warning("Some -ids- not present on the adjacency matrix:\n\t",
            paste0(dat.idvar[test],collapse = ", "),".")

  # Step 3: Checking attributes ------------------------------------------------
  # Creating the attributes (this depends on whether these are dynamic or not)
  vertex.attrs        <- vector("list", length(tran))
  names(vertex.attrs) <- tran

  if (length(timevar)) {
    for (i in tran) {
      # In order to access the slices, i do so using a character
      ichar <- as.character(i)

      vertex.attrs[[ichar]] <- merge(
        used.vertex,
        dat[is.na(dat[[timevar]]) | (dat[[timevar]] == i),],
        by = idvar,
        all.x=TRUE, sort=FALSE)

      # Sorting back
      vertex.attrs[[ichar]] <- vertex.attrs[[ichar]][
        order(vertex.attrs[[ichar]][["_original_sort"]]),]

      # Removing the id var, the per var and the toa var
      test <- colnames(vertex.attrs[[ichar]]) %in% c(varlist, "_original_sort")
      vertex.attrs[[ichar]] <- vertex.attrs[[ichar]][,which(!test),drop=FALSE]
    }
  } else {
    # Creating data.frame
    vertex.attrs <- merge(used.vertex, dat, by=idvar, all.x=TRUE, sort=FALSE)

    # Sorting back
    vertex.attrs <- vertex.attrs[
      order(vertex.attrs[["_original_sort"]]),]

    # Removing the idvar
    test <- colnames(vertex.attrs) %in% c(varlist, "_original_sort")
    vertex.attrs <- vertex.attrs[,which(!test),drop=FALSE]
  }

  # Times of Adoption vector
  toa <- unique(dat[,c(idvar, toavar)])
  toa <- merge(used.vertex, toa, by=idvar, all.x=TRUE,
               all.y=FALSE, sort=FALSE)

  # Sorting
  toa <- toa[order(toa[["_original_sort"]]),][[toavar]]

  # It should be of the same length as the used vertex
  if (length(toa) != nrow(used.vertex))
    stop("Multiple -toavar- by individual. Multiple adoption times are not ",
         "supported yet by the package.")

  # Step 4: Wrapping all together, creating the diffnet object -----------------
  if (length(timevar)) {
    as_diffnet(adjmat, toa=toa, t0 = dat.ran[1], t1=dat.ran[2],
               vertex.dyn.attrs = vertex.attrs,
               undirected=undirected, self=self,
               multiple=multiple)
  } else {
    as_diffnet(adjmat, toa=toa, t0 = dat.ran[1], t1=dat.ran[2],
               vertex.static.attrs = vertex.attrs,
               undirected=undirected, self=self,
               multiple=multiple)
  }

}
