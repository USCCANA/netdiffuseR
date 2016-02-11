rm(list=ls())

library(foreign)

kfamily <- read.dta("data-raw/kfp_v3.dta")
save(kfamily,file =  "data/kfamily.rdata")
