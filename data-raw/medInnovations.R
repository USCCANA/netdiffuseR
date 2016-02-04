rm(list=ls())
library(foreign)

# Preparing the data -----------------------------------------------------------
medInnovations <- read.dta("data-raw/mi_v2.dta")

save(medInnovations, file="data/medInnovations.rdata")
