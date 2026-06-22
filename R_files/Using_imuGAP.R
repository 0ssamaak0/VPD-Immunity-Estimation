library(imuGAP)
library(data.table)
library(plyr)

#step 1: load the data
school <- read.csv("~/Downloads/MMEDGit/VPD-Immunity-Estimation/data/all-schools.csv")

nc_locations <- read.csv("~/Downloads/MMEDGit/VPD-Immunity-Estimation/data/derived/nc_locations.csv")

nc_observations <- read.csv("~/Downloads/MMEDGit/VPD-Immunity-Estimation/data/derived/nc_observations.csv")
setDT(nc_observations)

nc_populations <- read.csv("~/Downloads/MMEDGit/VPD-Immunity-Estimation/data/derived/nc_populations.csv")

canonical_locations <-canonicalize_locations(nc_locations)

canonical_observations <- canonicalize_observations(nc_observations)

canonical_populations <- canonicalize_populations(nc_populations, nc_populations, nc_observations, nc_locations)

head(canonical_observations)

head(canonical_locations)

head(canonical_population)


#step 2

invalid_obs <- copy(nc_observations[, .(obs_id, loc_id, positive, sample_n, censored)])
invalid_obs[1, positive := sample_n + 10]

# This will fail validation and throw an error:
tryCatch(
  canonicalize_observations(invalid_obs),
  error = function(e) message("Caught expected error: ", e$message)
)

data("observations_sim", package = "imuGAP")
head(observations_sim[, .(obs_id, loc_id, positive, sample_n, censored)])

# Canonicalize and validate
canonical_observations <- canonicalize_observations(observations_sim)
head(canonical_observations)

nc_observations
