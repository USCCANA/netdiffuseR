#' Create a diffusion graph
#' @param graph An array
#' @param cumadopt nxT matrix
#' @param vcols A vector of size 2 with colors
#' @param mode Layout
#' @param layout.par Layout parameters
#' @param ... Further arguments to be passed to gplot
#' @return NULL
#' @export
plot_diffnet <- function(graph, cumadopt, vcols=c("blue","grey"), mode="fruchtermanreingold", layout.par=NULL,
                         mfcol.par=NULL, ...) {
  t <- dim(graph)[3]
  n <- dim(graph)[1]

  cols <- matrix(ncol=t, nrow=n)
  cumgraph <- matrix(0, n, n)
  for(i in 1:t) {
    cols[,i] <- ifelse(cumadopt[,i], vcols[1], vcols[2])
    cumgraph <- cumgraph + graph[,,i]
  }

  # Getting the coords
  fun <- getFromNamespace(paste0("gplot.layout.",mode), "sna")
  coords <- fun(cumgraph, layout.par)

  # Figuring out the dimension
  if (!length(mfcol.par)) {
    if (t<4) mfcol.par <- c(1,t)
    else if (t==4) mfcol.par <- c(2,2)
    else if (t==5) mfcol.par <- c(2,3)
    else if (t==6) mfcol.par <- c(2,3)
    else if (t==7) mfcol.par <- c(2,4)
    else if (t==8) mfcol.par <- c(2,4)
    else if (t==9) mfcol.par <- c(3,4)
    else if (t==10) mfcol.par <- c(3,4)
    else if (t==11) mfcol.par <- c(3,4)
    else if (t==12) mfcol.par <- c(3,4)
    else mfcol.par <- c(ceiling(t/4),4)
  }

  # Plotting
  curseed <- .Random.seed
  par(mfcol=c(1,3))
  for(i in 1:t)  {
    set.seed(curseed)
    sna::gplot(graph[,,i],displaylabels =  TRUE, vertex.col = cols[,i], coord=coords, ...)
  }
  par(mfcol=c(1,1))

}
