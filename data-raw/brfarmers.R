rm(list=ls())
library(foreign)
source("data-raw/listing_duplicated_columns.r")
# Preparing the data -----------------------------------------------------------
brfarmers <- read.dta("data-raw/bf_v2.dta")

listing_duplicated_columns(brfarmers)
# brfarmers <- subset(brfarmers, select = c(-ado))

save(brfarmers,file =  "data/brfarmers.rdata",
     compress = "xz")
