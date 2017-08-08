rm(list=ls())
library(foreign)

source("data-raw/listing_duplicated_columns.r")

# Preparing the data -----------------------------------------------------------
medInnovations <- read.dta("data-raw/mi_v2.dta")
# > listing_duplicated_columns(medInnovations)
# [,1]     [,2]
# [1,] "city"   "commun"
# [2,] "detail" "detail2"
# [3,] "ado"    "adopt"
# [4,] "ado"    "toa"
# [5,] "adopt"  "toa"
# medInnovations <- medInnovations[
#   ,
#   setdiff(
#     colnames(medInnovations),
#     c("commun", "detail2", "adopt", "ado"))
#   ]

save(medInnovations, file="data/medInnovations.rdata",
     compress = "xz")
