library(foreign)

x <- read.dta("data-raw/brfarmers.dta")

#               storage  display     value
# variable name   type   format      label      variable label
# ----------------------------------------------------------------------------
# net31           int    %8.0g                  nomination friend 1
# net32           int    %8.0g                  nomination friend 2
# net33           int    %8.0g                  nomination friend 3
# net21           int    %8.0g                  nomination influential 1
# net22           int    %8.0g                  nomination influential 2
# net23           int    %8.0g                  nomination influential 3
# net11           int    %8.0g                  nomination practice A
# net12           int    %8.0g                  nomination practice B
# net13           int    %8.0g                  nomination practice C
# net41           int    %8.0g                  nomination coop comm proj
# yr : Year of adoption (46' -> 66')

# Influential graph
# Rogers et al. (1970)
brfarmers <- subset(x, select=c(idold, net31, net32, net33, yr, village))
brfarmers$yr <- as.integer(brfarmers$yr) + 1900L

# Creating an ID
brfarmers$id <- with(brfarmers, idold + village*100L)
brfarmers$net31 <- with(brfarmers, net31 + village*100L)
brfarmers$net32 <- with(brfarmers, net32 + village*100L)
brfarmers$net33 <- with(brfarmers, net33 + village*100L)

# Reshaping and droping lost
brfarmers.long <- reshape(
  brfarmers, v.names= "net", varying = c("net31", "net32", "net33"),
  timevar = "level", idvar="id", direction="long", drop = c("idold"))

length(unique(unlist(subset(brfarmers.long, select=c(id,net)))))

library(diffusiontest)

# Creating the graph object
graph <- with(brfarmers.long, edgelist_to_adjmat(cbind(id, net), times=yr))
