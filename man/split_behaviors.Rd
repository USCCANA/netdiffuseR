% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdiffnet.r
\name{split_behaviors}
\alias{split_behaviors}
\title{Splitting behaviors}
\usage{
split_behaviors(diffnet_obj)
}
\arguments{
\item{diffnet_obj}{A multi-diffusion diffnet object.}
}
\value{
A list of diffnet objects. Each element represent a unique behavior.
}
\description{
Split each behavior within multi-diffusion diffnet object. The function gets
\code{toa}, \code{adopt}, \code{cumadopt}, and the \code{behavior} name from
each behavior, and returns a list where each element is a single behavior.
All the rest of the structure remains the same for each element in the list.
}
\examples{
# Running a multi-diffusion simulation
set.seed(1231)
diffnet_multi <- rdiffnet(50, 5, seed.p.adopt = list(0.1,0.1))

diffnet_multi_list <- split_behaviors(diffnet_multi)
diffnet_single <- diffnet_multi_list[[1]]

# You can now run standard functions for a single behavior
# Plotting single behavior
plot_diffnet(diffnet_single, slices = c(1, 3, 5))

}
\author{
George G. Vega Yon & Aníbal Olivera M.
}
