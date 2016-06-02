# This R script Demonstrates how to USE NetdifuseR to analyze the
# 3 classic network diffusion datasets and estimate cohesion vs
# Structural equivalence influences on adoption

rm(list = ls())
library(foreign)
library(netdiffuseR)

# Read know data
know16_att <- read.dta("../../clases/2016_PM542/netdiffuseR/know16.dta")
know16diffnet <- survey_to_diffnet(know16_att, idvar="id", netvars=c("nom1", "nom2", "nom3", "nom4",
                                                             "nom5", "nom6", "nom7"),
                               toavar="toa", warn.coercion = FALSE)
# Make sure it makes sense
summary(know16diffnet)
# Calculate exposures
know16diffnet[["cohexp"]] <- exposure(know16diffnet)
know16diffnet[["seexp2"]] <- exposure(know16diffnet, alt.graph = "se")

know16difnet <- diffnet.attrs(know16diffnet, as.df = TRUE)                  # Convert to dataframe
know16difnet$adopted <- as.integer(with(know16difnet, toa == per))          # Set adoption variable
know16difnet <- know16difnet[know16difnet$per <=  know16diffnet$toa, ]                   # Keep pre-adoption time only
mod_all <- as.formula(paste("adopted ~ factor(per) + cohexp + seexp2  "))
out_all <- glm(mod_all, data=know16difnet, family = binomial(link="logit"))
summary(out_all)

plot_diffnet(know16diffnet)
plot_threshold(know16diffnet,
               main="TOA and Threshold for Know16", vertex.label=NULL)

struct_test(know16diffnet, function(g) mean(threshold(g)), R=1000)

########################################################################
#                       The End                                        #
########################################################################


