library(imuGAP)
library(ggplot2)
base_dir <- "."
setwd(base_dir)

# Consistent report styling for coverage plots.
coverage_theme <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 3),
      plot.subtitle    = element_text(color = "gray35", size = base_size - 1),
      axis.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

county_colors <- c(
  "Scruggs" = "#1b9e77",
  "Simone"  = "#d95f02",
  "Watson"  = "#7570b3"
)

# Shared y-axis range for the two state-level coverage plots (baseline below,
# low-coverage scenario further down) so the drop is directly comparable.
# Spans both scenarios' 95% credible intervals (~0.78-0.95) with a little padding.
state_ylim <- c(0.77, 0.97)

# This driver compares the default synthetic coverage scenario with a low-coverage scenario.
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
# Use the default scenario first; the same fitting block is repeated below for sim4.
nc_observations<- sim$observations
nc_populations <- sim$populations
nc_locations <- sim$locations






# Define base directory (change as needed; here assumes script run from project root, or set to your data root)
base_dir <- file.path("data")
derived_dir <- file.path(base_dir, "derived")


# Canonicalize and validate
# The canonicalizers are the schema gate before running the Stan fit.
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

# Predict a snapshot for every location, age 1-18, and both dose levels.
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
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5, fill = "95% credible interval"),
              alpha = 0.22, color = NA) +
  geom_line(aes(y = q50, color = "Median"), linewidth = 1.3) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  coord_cartesian(ylim = state_ylim) +   # shared y-axis with low-coverage plot
  scale_color_manual(NULL, values = c("Median" = "#1b9e77")) +
  scale_fill_manual(NULL, values = c("95% credible interval" = "#1b9e77")) +
  coverage_theme() +
  labs(
    x = "Age", y = "Coverage",
    title = "State-level two-dose coverage",
    subtitle = "Median posterior estimate with 95% credible interval"
  )
ggsave("report/images/to_git/state_level_coverage_before.png",
       width = 1121 / 120, height = 835 / 120, dpi = 120)

summary_predict |>
  subset(loc_id %in% c("Scruggs", "Simone", "Watson") & dose == 2 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5, fill = loc_id), alpha = 0.16, color = NA) +
  geom_line(aes(y = q50, color = loc_id), linewidth = 1.25) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_color_manual("County", values = county_colors, aesthetics = c("color", "fill")) +
  coverage_theme() +
  labs(
    x = "Age", y = "Coverage",
    title = "County-level two-dose coverage",
    subtitle = "Median posterior estimate with uncertainty ribbons"
  )

  scruggs_schools <- locations_sim[parent_id == "Scruggs", loc_id]
summary_predict |>
  subset(
    loc_id %in% scruggs_schools & dose == 2 & age > 4
  ) |>
  ggplot() +
  geom_boxplot(aes(x = factor(age), y = q50),
               fill = "#9ecae1", color = "#08519c",
               linewidth = 0.6, outlier.alpha = 0.6) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  coverage_theme() +
  labs(
    x = "Age", y = "Median coverage",
    title = "School-level two-dose coverage in Scruggs County",
    subtitle = "Distribution across schools by age"
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
    alpha = 0.08, shape = 16, size = 1.3,
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.15
    )
  ) +
  geom_point(
    data = latent_ref,
    mapping = aes(shape = "Latent true value"),
    fill = "white",
    color = "black",
    size = 3,
    stroke = 1.1
  ) +
  scale_shape_manual(
    name = "",
    values = c("Latent true value" = 24)
  ) +
  scale_color_discrete(NULL, aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  coverage_theme() +
  labs(
    color = "School", x = "Age", y = "Coverage",
    title = "Posterior draws for selected schools",
    subtitle = "Jittered draws with latent true values highlighted"
  )


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
# Repeat validation for the low-coverage scenario before fitting it.
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
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5, fill = "95% credible interval"),
              alpha = 0.22, color = NA) +
  geom_line(aes(y = q50, color = "Median"), linewidth = 1.3) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  coord_cartesian(ylim = state_ylim) +   # shared y-axis with baseline plot
  scale_color_manual(NULL, values = c("Median" = "#1b9e77")) +
  scale_fill_manual(NULL, values = c("95% credible interval" = "#1b9e77")) +
  coverage_theme() +
  labs(
    x = "Age", y = "Coverage",
    title = "Low-coverage state-level two-dose coverage",
    subtitle = "Median posterior estimate with 95% credible interval"
  )
ggsave("report/images/to_git/after_line.png",
       width = 1121 / 120, height = 835 / 120, dpi = 120)

summary_predict |>
  subset(loc_id %in% c("Scruggs", "Simone", "Watson") & dose == 2 & age > 4) |>
  ggplot() +
  aes(x = age) +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5, fill = loc_id), alpha = 0.16, color = NA) +
  geom_line(aes(y = q50, color = loc_id), linewidth = 1.25) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_color_manual("County", values = county_colors, aesthetics = c("color", "fill")) +
  coverage_theme() +
  labs(
    x = "Age", y = "Coverage",
    title = "Low-coverage county-level two-dose coverage",
    subtitle = "Median posterior estimate with uncertainty ribbons"
  )

  scruggs_schools <- locations_sim[parent_id == "Scruggs", loc_id]
summary_predict |>
  subset(
    loc_id %in% scruggs_schools & dose == 2 & age > 4
  ) |>
  ggplot() +
  geom_boxplot(aes(x = factor(age), y = q50),
               fill = "#9ecae1", color = "#08519c",
               linewidth = 0.6, outlier.alpha = 0.6) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  coverage_theme() +
  labs(
    x = "Age", y = "Median coverage",
    title = "Low-coverage school-level two-dose coverage in Scruggs County",
    subtitle = "Distribution across schools by age"
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
    alpha = 0.08, shape = 16, size = 1.3,
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.15
    )
  ) +
  geom_point(
    data = latent_ref,
    mapping = aes(shape = "Latent true value"),
    fill = "white",
    color = "black",
    size = 3,
    stroke = 1.1
  ) +
  scale_shape_manual(
    name = "",
    values = c("Latent true value" = 24)
  ) +
  scale_color_discrete(NULL, aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = 5:18, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  coverage_theme() +
  labs(
    color = "School", x = "Age", y = "Coverage",
    title = "Low-coverage posterior draws for selected schools",
    subtitle = "Jittered draws with latent true values highlighted"
  )


saveRDS(predict_sim, "predict_sim4.rds")

#average for simulation 1

summary_dt_sim1 <- as.data.table(summary(predict_sim1))

avg_by_age_sim1 <- summary_dt_sim1[loc_id == "State" & dose == 2, .(
  mean_coverage = round(mean(mean),  3),
  lower_ci      = round(mean(q2_5),  3),
  upper_ci      = round(mean(q97_5), 3)
), by = age][order(age)]

# Non-zero only
# These state-level age averages feed the Reed-Frost phase as background immunity.
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


