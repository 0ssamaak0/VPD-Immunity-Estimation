library(imuGAP)
library(data.table)
library(plyr)
library(ggplot2)

#step 1: load the data

# Start from imuGAP's built-in synthetic fixture to learn the expected schemas.
data("locations_sim", package = "imuGAP")
head(locations_sim)

# Canonicalize and validate
canonical_locations <- canonicalize_locations(locations_sim)
head(canonical_locations)

data("observations_sim", package = "imuGAP")
head(observations_sim[, .(obs_id, loc_id, positive, sample_n, censored)])

# Canonicalize and validate
canonical_observations <- canonicalize_observations(observations_sim)
head(canonical_observations)

data("populations_sim", package = "imuGAP")
head(populations_sim)

# Canonicalize and validate
canonical_populations <- canonicalize_populations(
  populations_sim, observations_sim, locations_sim
)
head(canonical_populations)
# Define base directory (change as needed; here assumes script run from project root, or set to your data root)
base_dir <- file.path("data")
derived_dir <- file.path(base_dir, "derived")




data("locations_sim", package = "imuGAP")
head(locations_sim)

# Canonicalize and validate
canonical_locations <- canonicalize_locations(locations_sim)
head(canonical_locations)

# Reuse the package fixture as the working "nc_*" inputs for this baseline fit.
nc_locations <- locations_sim
nc_observations <- observations_sim
nc_populations <- populations_sim

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

saveRDS(fit_sim, "fit_sim2.rds")

##target_simulation

# Build prediction targets for all hierarchy levels, ages, and doses.
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

mcv1_baseline <- subset(summary_predict, dose == 1)
# Save the one-dose baseline for downstream comparison with reduced-coverage scenarios.
write.csv(mcv1_baseline, "mcv1_coverage_baseline.csv", row.names = FALSE)

summary_predict |>
  subset(loc_id == "State" & dose == 1 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_line(aes(y = q50)) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  labs(x = "Age", y = "Coverage", title = "State-level MCV1 (one-dose) coverage")

summary_predict |>
  subset(loc_id %in% c("Scruggs", "Simone", "Watson") & dose == 1 & age > 4) |>
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
    x = "Age", y = "Coverage", title = "County-level MCV1 (one-dose) coverage"
  )

  scruggs_schools <- locations_sim[parent_id == "Scruggs", loc_id]
summary_predict |>
  subset(
    loc_id %in% scruggs_schools & dose == 1 & age > 4
  ) |>
  ggplot() +
  geom_boxplot(aes(x = factor(age), y = q50)) +
  theme_bw() +
  labs(
    x = "Age", y = "Coverage",
    title = "Scruggs County school-level MCV1 (one-dose) coverage"
  )


  schools <- c(
  "Towhee Children's Academy", # ~380 per grade
  "Flycatcher Elementary", # ~110 per grade
  "Sparrow School" # ~60 per grade
)

# Subset to targets of interest (all retained posterior draws)
predict_sub <- predict_sim |>
  subset(loc_id %in% schools & dose == 1 & age > 4)

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