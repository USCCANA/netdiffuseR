rm(list=ls())
library(foreign)

# Preparing the data -----------------------------------------------------------
brfarmers <- read.dta("data-raw/brfarmers.dta")
save(brfarmers,file =  "data/brfarmers.rdata",
     compress = "xz")
