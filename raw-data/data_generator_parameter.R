library(imuGAP)
library(ggplot2)
base_dir <- "."
setwd(base_dir)
source("raw-data/simulate_imuGAP_data.R")

# ── Default (matches the package fixture exactly) ──────────────
sim <- simulate_imuGAP_data()
print_sim_summary(sim)

# ── Try higher school-level noise ──────────────────────────────
sim2 <- simulate_imuGAP_data(sigma_sch = 1.5, sigma_cnty = 0.6)
print_sim_summary(sim2)

# ── Fewer schools, different seed ──────────────────────────────
sim3 <- simulate_imuGAP_data(seed = 42, n_schools = c(4, 4, 4))
print_sim_summary(sim3)

# ── Lower overall coverage (shift phi_st down) ─────────────────
low_phi <- sim$params  # just to see the default, then modify
sim4 <- simulate_imuGAP_data(phi_st = rep(0.5, 30))

# ── Access the tables directly ─────────────────────────────────
nc_observations<- sim$observations
nc_populations <- sim$populations
nc_locations <- sim$locations






# Define base directory (change as needed; here assumes script run from project root, or set to your data root)
base_dir <- file.path("data")
derived_dir <- file.path(base_dir, "derived")


# Canonicalize and validate
canonical_locations <- canonicalize_locations(nc_locations)
head(canonical_locations)


canonical_locations <-canonicalize_locations(nc_locations)

canonical_observations <- canonicalize_observations(nc_observations)

canonical_populations <- canonicalize_populations(populations = nc_populations, observations = nc_observations, locations = nc_locations)




### fit the model


fit_sim <- sampling(
  
  nc_observations,
  nc_populations,
  nc_locations,
  stan_opts = stan_options(
    iter = 2000, chains = 4, refresh = 0, seed = 1L
  )
)

## Save the fitted model to an RDS file

saveRDS(fit_sim, "fit_sim1.rds")

##target_simulation

target_sim <- create_target(
  fit = fit_sim, location = unique(nc_locations$loc_id), age = 1:18,
  cohort = max(nc_populations$cohort) - 18, dose = c(1, 2), mode = "snapshot"
)
head(target_sim)

nc_observations

##prediction 

predict_sim <- predict(object = fit_sim, target = target_sim, posterior_size = 100)



# Calculate the posterior mean coverage probability for each location and dose at age 5
summary_predict <- summary(predict_sim)
head(summary_predict)


summary_predict |>
  subset(loc_id == "State" & dose == 2 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_line(aes(y = q50)) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  labs(x = "Age", y = "Coverage", title = "State-level two dose coverage")

summary_predict |>
  subset(loc_id %in% c("Scruggs", "Simone", "Watson") & dose == 2 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_line(aes(y = q50, color = loc_id)) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5, fill = loc_id), alpha = 0.2) +
  theme_bw() +
  theme(
    legend.position = "inside", legend.position.inside = c(.2, 0.05),
    legend.justification.inside = c(0, 0)
  ) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_color_discrete(NULL, aesthetics = c("color", "fill")) +
  labs(
    x = "Age", y = "Coverage", title = "County-level two dose coverage"
  )

  scruggs_schools <- locations_sim[parent_id == "Scruggs", loc_id]
summary_predict |>
  subset(
    loc_id %in% scruggs_schools & dose == 2 & age > 4
  ) |>
  ggplot() +
  geom_boxplot(aes(x = factor(age), y = q50)) +
  theme_bw() +
  labs(
    x = "Age", y = "Coverage",
    title = "Distribution of School-Level Coverage For Scruggs County"
  )


  schools <- c(
  "Towhee Children's Academy", # ~380 per grade
  "Flycatcher Elementary", # ~110 per grade
  "Sparrow School" # ~60 per grade
)

# Subset to targets of interest (all retained posterior draws)
predict_sub <- predict_sim |>
  subset(loc_id %in% schools & dose == 2 & age > 4)

# Get the pre-computed background coverage matching the subsetted target
latent_ref <- copy(predict_sub$target)
latent_ref$coverage <- latent_params_sim$coverage[predict_sub$target$obs_id]


# Convert predictions to a long-format data.frame
draws_df <- as.data.frame(predict_sub)

# Now plot it all
ggplot() +
  aes(age, coverage, color = loc_id) +
  geom_point(
    data = draws_df,
    alpha = 1 / 256, shape = 16,
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.15
    )
  ) +
  geom_point(
    data = latent_ref,
    mapping = aes(shape = "True value"),
    fill = NA
  ) +
  theme_bw() +
  scale_shape_manual(
    name = "",
    values = c("True value" = 24)
  ) +
  scale_color_discrete(NULL, aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  theme(legend.position = "bottom") +
  labs(color = "School", x = "Age", y = "Coverage")


sim1_fit <- fit_sim
predict_sim1 <- predict_sim
saveRDS(predict_sim, "predict_sim1.rds")

saveRDS(sim1_fit, "fit_sim1.rds")


nc_observations <- sim4$observations
nc_populations <- sim4$populations
nc_locations <- sim4$locations






# Define base directory (change as needed; here assumes script run from project root, or set to your data root)
base_dir <- file.path("data")
derived_dir <- file.path(base_dir, "derived")


# Canonicalize and validate
canonical_locations <- canonicalize_locations(nc_locations)
head(canonical_locations)


canonical_locations <-canonicalize_locations(nc_locations)

canonical_observations <- canonicalize_observations(nc_observations)

canonical_populations <- canonicalize_populations(populations = nc_populations, observations = nc_observations, locations = nc_locations)




### fit the model


fit_sim <- sampling(
  
  nc_observations,
  nc_populations,
  nc_locations,
  stan_opts = stan_options(
    iter = 2000, chains = 4, refresh = 0, seed = 1L
  )
)

## Save the fitted model to an RDS file

saveRDS(fit_sim, "fit_sim4.rds")

fit_sim4 <-fit_sim

##target_simulation

target_sim <- create_target(
  fit = fit_sim, location = unique(nc_locations$loc_id), age = 1:18,
  cohort = max(nc_populations$cohort) - 18, dose = c(1, 2), mode = "snapshot"
)
head(target_sim)

nc_observations

##prediction 

predict_sim <- predict(object = fit_sim, target = target_sim, posterior_size = 100)

predict_sim4 <-predict_sim

# Calculate the posterior mean coverage probability for each location and dose at age 5
summary_predict <- summary(predict_sim)
head(summary_predict)


summary_predict |>
  subset(loc_id == "State" & dose == 2 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_line(aes(y = q50)) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  labs(x = "Age", y = "Coverage", title = "State-level two dose coverage")

summary_predict |>
  subset(loc_id %in% c("Scruggs", "Simone", "Watson") & dose == 2 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_line(aes(y = q50, color = loc_id)) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5, fill = loc_id), alpha = 0.2) +
  theme_bw() +
  theme(
    legend.position = "inside", legend.position.inside = c(.2, 0.05),
    legend.justification.inside = c(0, 0)
  ) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_color_discrete(NULL, aesthetics = c("color", "fill")) +
  labs(
    x = "Age", y = "Coverage", title = "County-level two dose coverage"
  )

  scruggs_schools <- locations_sim[parent_id == "Scruggs", loc_id]
summary_predict |>
  subset(
    loc_id %in% scruggs_schools & dose == 2 & age > 4
  ) |>
  ggplot() +
  geom_boxplot(aes(x = factor(age), y = q50)) +
  theme_bw() +
  labs(
    x = "Age", y = "Coverage",
    title = "Distribution of School-Level Coverage For Scruggs County"
  )


  schools <- c(
  "Towhee Children's Academy", # ~380 per grade
  "Flycatcher Elementary", # ~110 per grade
  "Sparrow School" # ~60 per grade
)

# Subset to targets of interest (all retained posterior draws)
predict_sub <- predict_sim |>
  subset(loc_id %in% schools & dose == 2 & age > 4)

# Get the pre-computed background coverage matching the subsetted target
latent_ref <- copy(predict_sub$target)
latent_ref$coverage <- latent_params_sim$coverage[predict_sub$target$obs_id]


# Convert predictions to a long-format data.frame
draws_df <- as.data.frame(predict_sub)

# Now plot it all
ggplot() +
  aes(age, coverage, color = loc_id) +
  geom_point(
    data = draws_df,
    alpha = 1 / 256, shape = 16,
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.15
    )
  ) +
  geom_point(
    data = latent_ref,
    mapping = aes(shape = "True value"),
    fill = NA
  ) +
  theme_bw() +
  scale_shape_manual(
    name = "",
    values = c("True value" = 24)
  ) +
  scale_color_discrete(NULL, aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  theme(legend.position = "bottom") +
  labs(color = "School", x = "Age", y = "Coverage")


saveRDS(predict_sim, "predict_sim4.rds")

#average for simulation 1

summary_dt_sim1 <- as.data.table(summary(predict_sim1))

avg_by_age_sim1 <- summary_dt_sim1[loc_id == "State" & dose == 2, .(
  mean_coverage = round(mean(mean),  3),
  lower_ci      = round(mean(q2_5),  3),
  upper_ci      = round(mean(q97_5), 3)
), by = age][order(age)]

# Non-zero only
average_1<-avg_by_age_sim1[mean_coverage > 0]
average_1

#avergae for simulation 4
summary_dt_sim4 <- as.data.table(summary(predict_sim4))

avg_by_age_sim4 <- summary_dt_sim4[loc_id == "State" & dose == 2, .(
  mean_coverage = round(mean(mean),  3),
  lower_ci      = round(mean(q2_5),  3),
  upper_ci      = round(mean(q97_5), 3)
), by = age][order(age)]

# Non-zero only
average_4<-avg_by_age_sim4[mean_coverage > 0]


