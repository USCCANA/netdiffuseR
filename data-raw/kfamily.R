rm(list=ls())

library(foreign)
source("data-raw/listing_duplicated_columns.r")
kfamily <- read.dta("data-raw/kfp_v3_labels_fixed.dta", convert.factors = FALSE)

# [,1]      [,2]
# "village" "area1"
# "village" "area2"
# "village" "area3"
# "village" "area4"
# "village" "area5"
# "village" "area6"
# "village" "area7"
# "village" "commun"
# "id"      "id1"
# "id"      "id2"
# "id"      "id3"
# "id"      "id4"
# "id"      "id5"
# "id"      "id6"
# "id"      "id7"
# "recno1"  "studno1"
# "recno1"  "studno2"
# "recno1"  "studno4"
# "recno1"  "studno5"
# "recno1"  "studno6"
# "recno1"  "studno7"
# "studno1" "studno2"
# "studno1" "studno4"
# "studno1" "studno5"
# "studno1" "studno6"
# "studno1" "studno7"
# "area1"   "area2"
# "area1"   "area3"
# "area1"   "area4"
# "area1"   "area5"
# "area1"   "area6"
# "area1"   "area7"
# "area1"   "commun"
# "id1"     "id2"
# "id1"     "id3"
# "id1"     "id4"
# "id1"     "id5"
# "id1"     "id6"
# "id1"     "id7"
# "studno2" "studno4"
# "studno2" "studno5"
# "studno2" "studno6"
# "studno2" "studno7"
# "area2"   "area3"
# "area2"   "area4"
# "area2"   "area5"
# "area2"   "area6"
# "area2"   "area7"
# "area2"   "commun"
# "id2"     "id3"
# "id2"     "id4"
# "id2"     "id5"
# "id2"     "id6"
# "id2"     "id7"
# "area3"   "area4"
# "area3"   "area5"
# "area3"   "area6"
# "area3"   "area7"
# "area3"   "commun"
# "id3"     "id4"
# "id3"     "id5"
# "id3"     "id6"
# "id3"     "id7"
# "studno4" "studno5"
# "studno4" "studno6"
# "studno4" "studno7"
# "area4"   "area5"
# "area4"   "area6"
# "area4"   "area7"
# "area4"   "commun"
# "id4"     "id5"
# "id4"     "id6"
# "id4"     "id7"
# "studno5" "studno6"
# "studno5" "studno7"
# "area5"   "area6"
# "area5"   "area7"
# "area5"   "commun"
# "id5"     "id6"
# "id5"     "id7"
# "studno6" "studno7"
# "area6"   "area7"
# "area6"   "commun"
# "id6"     "id7"
# "area7"   "commun"
# "awe2t9"  "awe2t10"
# "awe2t9"  "awe2t12"
# "awe3t9"  "awe3t10"
# "awe3t9"  "awe3t11"
# "awe3t9"  "awe3t12"
# "awe2t10" "awe2t12"
# "awe3t10" "awe3t11"
# "awe3t10" "awe3t12"
# "awe3t11" "awe3t12"
# "ado"     "toa"

# All these are duplicated
toremove <- c("area1",
"area2"  ,
"area3"  ,
"area4"  ,
"area5"  ,
"area6"  ,
"area7"  ,
"commun" ,
"id1"    ,
"id2"    ,
"id3"    ,
"id4"    ,
"id5"    ,
"id6"    ,
"id7"    ,
"studno1",
"studno2",
"studno4",
"studno5",
"studno6",
"studno7",
"studno2",
"studno4",
"studno5",
"studno6",
"studno7",
"area2"  ,
"area3"  ,
"area4"  ,
"area5"  ,
"area6"  ,
"area7"  ,
"commun" ,
"id2"    ,
"id3"    ,
"id4"    ,
"id5"    ,
"id6"    ,
"id7"    ,
"studno4",
"studno5",
"studno6",
"studno7",
"area3"  ,
"area4"  ,
"area5"  ,
"area6"  ,
"area7"  ,
"commun" ,
"id3"    ,
"id4"    ,
"id5"    ,
"id6"    ,
"id7"    ,
"area4"  ,
"area5"  ,
"area6"  ,
"area7"  ,
"commun" ,
"id4"    ,
"id5"    ,
"id6"    ,
"id7"    ,
"studno5",
"studno6",
"studno7",
"area5"  ,
"area6"  ,
"area7"  ,
"commun" ,
"id5"    ,
"id6"    ,
"id7"    ,
"studno6",
"studno7",
"area6"  ,
"area7"  ,
"commun" ,
"id6"    ,
"id7"    ,
"studno7",
"area7"  ,
"commun" ,
"id7"    ,
"commun" ,
"awe2t10",
"awe2t12",
"awe3t10",
"awe3t11",
"awe3t12",
"awe2t12",
"awe3t11",
"awe3t12",
"awe3t12", "ado")

# kfamily <- kfamily[,setdiff(colnames(kfamily), toremove)]
save(kfamily, file =  "data/kfamily.rdata",
     compress = "xz")
