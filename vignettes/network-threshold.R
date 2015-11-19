## ----Setup, echo=FALSE---------------------------------------------------
doctype <- knitr::opts_knit$get("rmarkdown.pandoc.to")

# Whether is html or not
ishtml <- !is.null(doctype) && doctype=='html'
knitr::opts_chunk$set(fig.height=5, fig.width=8, fig.align="center", 
                      out.width=ifelse(ishtml, 600, ".8\\linewidth"))

## ------------------------------------------------------------------------
library(netdiffuseR)

## ----DGP-----------------------------------------------------------------
# Generating random data
set.seed(123)
graph <- rand_graph(n=1000,t=5)
toa <- sample(2000:(2000+dim(graph)[3]-1), dim(graph)[1], TRUE)

# Threshold
adopt <- toa_mat(toa)
expo <- exposure(graph, adopt$cumadopt)
sunflowerplot(threshold(expo, toa))

