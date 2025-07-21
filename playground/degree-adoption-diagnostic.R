#' This function adds degree attributes to the diffnet object, creates a panel
#' dataset, and models the likelihood of adoption in a given period as a
#' function of the node's centrality.
#'
#' @param diffnet_obj An object of class 'diffnet'.
#' @return A summary of the fitted glm object.
#'
run_longitudinal_model <- function(diffnet_obj) {

  # adding degree attributes to the diffnet object
  diffnet_obj[["indegree"]] <- dgr(diffnet_obj, "indegree")
  diffnet_obj[["outdegree"]] <- dgr(diffnet_obj, "outdegree")

  # create a panel dataset.
  panel_data <- as.data.frame(diffnet_obj)

  # create a binary variable: 1 for adoption in the period, 0 otherwise
  panel_data$adopt_in_period <- ifelse(panel_data$toa == panel_data$per, 1, 0)

  # This keeps the individual in the dataset only while they have not adopted.
  at_risk_data <- subset(panel_data, is.na(toa) | toa >= per)

  # run the logit model
  model <- glm( adopt_in_period ~ indegree + outdegree,
                data = at_risk_data,
                family = binomial(link = "logit")
  )

  # 6. DEVOLVER EL RESUMEN DEL MODELO
  return(summary(model))
}

library(netdiffuseR)

data(medInnovationsDiffNet)
data(kfamilyDiffNet)

run_longitudinal_model(medInnovationsDiffNet) # Only in-degree is ""significant""
run_longitudinal_model(kfamilyDiffNet)

# Some plots ----

# Rate of adoption plot
# Shows the porcentual growth in the base of adopters, a value > 1.0 indicates explosive growth in a period.
plot_adoption_rate <- function(diffnet_obj, main_title, leader_metric = "outdegree") {
  avg_degree <- rowMeans(dgr(diffnet_obj, leader_metric), na.rm = TRUE)
  leader_ids <- diffnet_obj$meta$ids[avg_degree >= quantile(avg_degree, 0.9, na.rm = TRUE)]
  follower_ids <- setdiff(diffnet_obj$meta$ids, leader_ids)

  rate_leaders <- cumulative_adopt_count(diffnet_obj[leader_ids, , ])[3,]
  rate_followers <- cumulative_adopt_count(diffnet_obj[follower_ids, , ])[3,]

  plot(rate_leaders, type = "l", col = "blue", lwd = 2,
       ylim = c(0, max(c(rate_leaders, rate_followers), na.rm=T, finite=T)),
       xlab = "Time", ylab = "Adoption Rate", main = main_title)
  lines(rate_followers, type = "l", col = "red", lwd = 2, lty = 2)
  legend("bottomright", legend = c("Opinion Leaders (Top 10%)", "Followers"),
         col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")
}

# in-degree seems to keep higher adoption rates for longer
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot_adoption_rate(medInnovationsDiffNet, "Adoption Rate (Leaders: in-degree) Medical Innovations", "indegree")
plot_adoption_rate(medInnovationsDiffNet, "Adoption Rate (Leaders: out-degree) Medical Innovations", "outdegree")
par(mfrow = c(1, 1))

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot_adoption_rate(kfamilyDiffNet, "Adoption Rate (Leaders: in-degree) KFP", "indegree")
plot_adoption_rate(kfamilyDiffNet, "Adoption Rate (Leaders: out-degree) KFP", "outdegree")
par(mfrow = c(1, 1))

# New adopters plot
# The solid line shows an adjusted trend
plot_new_adopters <- function(diffnet_obj, main_title, leader_metric = "outdegree") {
  # identify leaders and followers
  avg_degree <- rowMeans(dgr(diffnet_obj, leader_metric), na.rm = TRUE)
  leader_ids <- diffnet_obj$meta$ids[avg_degree >= quantile(avg_degree, 0.9, na.rm = TRUE)]
  follower_ids <- setdiff(diffnet_obj$meta$ids, leader_ids)

  # calculate the proportion of new adopters
  cum_adopt_leaders <- cumulative_adopt_count(diffnet_obj[leader_ids, , ])[1, ]
  cum_adopt_followers <- cumulative_adopt_count(diffnet_obj[follower_ids, , ])[1, ]
  prop_new_leaders <- c(cum_adopt_leaders[1], diff(cum_adopt_leaders)) / length(leader_ids)
  prop_new_followers <- c(cum_adopt_followers[1], diff(cum_adopt_followers)) / length(follower_ids)

  # plot the proportions
  time_periods <- 1:length(prop_new_leaders)
  plot(time_periods, prop_new_leaders, type = "l", lty = 3, col = "blue",
       ylim = c(0, max(c(prop_new_leaders, prop_new_followers), na.rm = TRUE) * 1.1),
       xlab = "Time", ylab = "Proportion of New Adopters", main = main_title)
  lines(time_periods, prop_new_followers, type = "l", lty = 3, col = "red")

  abline(lm(prop_new_leaders ~ time_periods), col = "blue", lwd = 2.5)
  abline(lm(prop_new_followers ~ time_periods), col = "red", lwd = 2.5)

  legend("topright", legend = c("Opinion Leaders (Top 10%)", "Followers"),
         col = c("blue", "red"), lwd = 2, lty = 1, bty = "n")
}

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot_new_adopters(medInnovationsDiffNet, "New Adopters (Leaders: in-degree) Medical Innovations", "indegree")
plot_new_adopters(medInnovationsDiffNet, "New Adopters (Leaders: out-degree) Medical Innovations", "outdegree")
par(mfrow = c(1, 1))

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot_new_adopters(kfamilyDiffNet, "New Adopters (Leaders: in-degree) KFP", "indegree")
plot_new_adopters(kfamilyDiffNet, "New Adopters (Leaders: out-degree) KFP", "outdegree")
par(mfrow = c(1, 1))
