FROM rocker/tidyverse:4.5.0

# RUN \
#   echo 'options(repos=c(CRAN="https://packagemanager.posit.co/cran/__linux__/bookworm/latest"))' >> ~/.Rprofile && \
#   Rscript --vanilla -e 'getOption("repos")'

RUN install2.r --error languageserver && installGithub.r nx10/httpgd

RUN apt-get update && apt-get install --no-install-recommends -y valgrind gdb \
  libglpk-dev

# Adding R packages
RUN install2.r --error Rcpp sna network networkDynamic Matrix MASS MatchIt SparseM igraph \
  viridisLite covr testthat knitr rmarkdown ape RSiena survival RcppArmadillo

CMD ["bash"]
