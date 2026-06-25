library(imuGAP)
library(data.table)
library(plyr)

#step 1: load the data


# Define base directory (change as needed; here assumes script run from project root, or set to your data root)
base_dir <- file.path("data")
derived_dir <- file.path(base_dir, "derived")

# These CSVs are the Python-prepared NC tables in imuGAP's three-table format.
school <- read.csv(file.path(base_dir, "all-schools.csv"))

nc_locations <- read.csv(file.path(derived_dir, "nc_locations.csv"),
                         na.strings = c("", "NA"))

nc_observations <- read.csv(file.path(derived_dir, "nc_observations.csv"))
setDT(nc_observations)

nc_populations <- read.csv(file.path(derived_dir, "nc_populations.csv"))

# Canonicalization checks hierarchy layers and observation/population consistency.
canonical_locations <-canonicalize_locations(nc_locations)

canonical_observations <- canonicalize_observations(nc_observations)

canonical_populations <- canonicalize_populations(populations = nc_populations, observations = nc_observations, locations = nc_locations)


head(canonical_observations)

head(canonical_locations)

head(canonical_populations)

max_layer <- max(canonical_locations$layer)

head(nc_locations)

layer_4_list <- canonical_locations[layer == max_layer, loc_id]



### fit the model


# Fit on the prepared NC inputs after they pass the canonicalizer checks.
fit_sim <- sampling(
  
  nc_observations,
  nc_populations,
  nc_locations,
  stan_opts = stan_options(
    iter = 2000, chains = 4, refresh = 0, seed = 1L
  )
)
canonical_observations <- canonicalize_observations(observations_sim)
head(canonical_observations)

nc_observations
