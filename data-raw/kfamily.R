rm(list=ls())

library(foreign)

kfamily <- read.dta("data-raw/kfp_v3_labels_fixed.dta")
save(kfamily,file =  "data/kfamily.rdata")
