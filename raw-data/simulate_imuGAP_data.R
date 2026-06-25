## simulate_imuGAP_data.R
##
## Wraps the data-generation logic into a single callable function so you can
## experiment with parameters interactively before committing to the full
## Stan fit.  No Stan, no file I/O, no usethis required.
##
## Usage:
##   source("simulate_imuGAP_data.R")
##   sim <- simulate_imuGAP_data()          # defaults match the package fixture
##   sim <- simulate_imuGAP_data(seed = 42, sigma_sch = 1.2, n_schools = c(5,5,4))
##   print_sim_summary(sim)

library(data.table)
library(dplyr)

# Base R replacement for EnvStats::rlnormTrunc() — no extra package needed.
# Draws n values from a truncated log-normal via rejection sampling.
.rlnorm_trunc <- function(n, meanlog, sdlog, min = -Inf, max = Inf) {
  out <- numeric(n)
  i   <- 0L
  while (i < n) {
    x    <- rlnorm(n - i, meanlog, sdlog)
    x    <- x[x >= min & x <= max]
    take <- min(length(x), n - i)
    out[(i + 1L):(i + take)] <- x[seq_len(take)]
    i    <- i + take
  }
  out
}

# =============================================================================
#' Simulate imuGAP input data
# =============================================================================
simulate_imuGAP_data <- function(
  seed              = 93254L,
  n_yr              = 33L,
  n_cohort          = 30L,
  phi_st            = c(
    0.8401733, 0.8458791, 0.8515769, 0.8572586, 0.8629160,
    0.8685411, 0.8741259, 0.8796623, 0.8851422, 0.8905575,
    0.8958959, 0.9011275, 0.9062182, 0.9111339, 0.9158404,
    0.9203035, 0.9244892, 0.9283632, 0.9318916, 0.9350400,
    0.9377351, 0.9397467, 0.9408054, 0.9407130, 0.9395576,
    0.9375024, 0.9347246, 0.9314054, 0.9277256, 0.9298663
  ),
  lambda            = c(2.8, 3.0),
  sigma_sch         = 0.8,
  sigma_cnty        = 0.4,
  other_vax_reduction = 0.95,
  n_schools         = c(9L, 8L, 7L),
  dose_schedule     = c(1L, 4L),
  county_names      = c("Scruggs", "Simone", "Watson"),
  school_names      = c(
    "Chickadee Elementary",    "Nuthatch Academy",
    "Blue Heron School",       "Flycatcher Elementary",
    "Bluebird Learning Center","Catbird Academy",
    "Finch Elementary",        "Sparrow School",
    "Towhee Children's Academy","Warbler Elementary",
    "Egret Elementary",        "Cardinal Academy",
    "Bunting School",          "Tanager Academy",
    "Oriole Youth Academy",    "Grosbeak Learning Center",
    "Junco Elementary",        "Meadowlark School",
    "Goldfinch Elementary",    "Mockingbird Academy",
    "Kinglet Learning Center", "Vireo School",
    "Kingfisher Academy",      "Cormorant Elementary"
  )
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    "phi_st must have length == n_cohort"       = length(phi_st) == n_cohort,
    "phi_st values must be in (0, 1)"           = all(phi_st > 0 & phi_st < 1),
    "lambda must have length 2"                  = length(lambda) == 2L,
    "lambda values must be positive"             = all(lambda > 0),
    "n_schools must have length 3"              = length(n_schools) == 3L,
    "n_schools must sum to <= length(school_names)" =
      sum(n_schools) <= length(school_names),
    "county_names must have length 3"           = length(county_names) == 3L,
    "sigma_sch must be positive"                = sigma_sch > 0,
    "sigma_cnty must be positive"               = sigma_cnty > 0,
    "other_vax_reduction must be in (0, 1]"     =
      other_vax_reduction > 0 & other_vax_reduction <= 1
  )

  set.seed(seed)

  n_doses <- length(lambda)

  # ── Dose uptake matrix ──────────────────────────────────────────────────────
  doses <- matrix(0, ncol = n_doses, nrow = n_yr)
  for (i in seq_along(dose_schedule)) {
    doses[(dose_schedule[i] + 1):nrow(doses), i] <- 1
  }

  # cov stores the expected cumulative uptake by age/year for each dose.
  cov <- matrix(nrow = n_yr, ncol = n_doses)
  cov[1, ] <- 0
  for (d in seq_len(n_doses)) {
    ref      <- if (d == 1L) rep(1, n_yr) else cov[, d - 1L]
    survival <- 1 - exp(-lambda[d] * doses[, d])
    for (i in 2:n_yr) {
      cov[i, d] <- cov[i - 1, d] + (ref[i] - cov[i - 1, d]) * survival[i]
    }
  }

  # ── Random effects ──────────────────────────────────────────────────────────
  n_total_schools <- sum(n_schools)
  cnty_offset     <- rnorm(3L,              0, sigma_cnty)
  sch_offset      <- rnorm(n_total_schools, 0, sigma_sch)
  cnty_ids        <- rep(seq_along(n_schools), times = n_schools)

  # ── Child vax view ──────────────────────────────────────────────────────────
  n24 <- round(runif(n_cohort, 250, 450))
  n36 <- round(runif(n_cohort, 250, 450))

  sim_child <- dplyr::bind_rows(
    data.frame(
      vaxview_type = "child", Year = 1:n_cohort, age = "24 months",
      y_obs = rbinom(n_cohort, n24, phi_st * cov[2, 1] * other_vax_reduction),
      y_smp = n24, censored = 1L
    ),
    data.frame(
      vaxview_type = "child", Year = 1:n_cohort, age = "36 months",
      y_obs = rbinom(n_cohort, n36, phi_st * cov[3, 1] * other_vax_reduction),
      y_smp = n36, censored = 1L
    )
  )

  # ── Teen vax view ───────────────────────────────────────────────────────────
  teen_yrs <- 18:30
  sim_teen <- data.frame(
    vaxview_type = "teen", Year = teen_yrs,
    y_obs = numeric(length(teen_yrs)), y_smp = numeric(length(teen_yrs))
  )
  for (i in seq_len(nrow(sim_teen))) {
    samp_size          <- round(runif(1, 40, 70))
    sim_teen$y_smp[i]  <- samp_size * 5
    sim_teen$y_obs[i]  <- sum(
      rbinom(5, samp_size,
             phi_st[(teen_yrs[i] - 17):(teen_yrs[i] - 13)] * cov[18:14, 2])
    )
  }

  # ── School-level data ───────────────────────────────────────────────────────
  sch_yrs   <- 6:30
  nsch_base <- round(
    .rlnorm_trunc(n_total_schools, log(75), log(2.5), min = 10, max = 450)
  )

  kg_sim_full <- vector("list", n_total_schools)
  for (s in seq_len(n_total_schools)) {
    nsch    <- numeric(length(sch_yrs))
    nsch[1] <- nsch_base[s]
    for (y in 2:length(sch_yrs)) {
      nsch[y] <- max(4, round(runif(1, nsch[y - 1] - 5, nsch[y - 1] + 5)))
    }
    offset   <- sch_offset[s] + cnty_offset[cnty_ids[s]]
    # Random effects shift each school's coverage on the logit scale.
    cov_temp <- plogis(qlogis(phi_st[sch_yrs - 5]) + offset) * cov[5, 2]

    kg_sim_full[[s]] <- data.frame(
      year        = sch_yrs,
      enc_unit_id = cnty_ids[s] + 1L,
      unit_id     = s,
      y_obs       = rbinom(length(sch_yrs), nsch, cov_temp),
      y_smp       = nsch
    )
  }
  kg_sim_full <- dplyr::bind_rows(kg_sim_full)

  # ── School vax view (aggregate) ─────────────────────────────────────────────
  annual_tots <- kg_sim_full |>
    dplyr::group_by(year) |>
    dplyr::summarize(tot_enr = sum(y_smp), tot_vax = sum(y_obs), .groups = "drop")

  # The school survey is a noisy aggregate of the simulated school records.
  if (nrow(annual_tots) > 0) {
    probs <- rep_len(phi_st[sch_yrs - 5] * cov[5, 2], nrow(annual_tots))
    sim_school <- data.frame(
      vaxview_type = rep("school", nrow(annual_tots)),
      Year         = annual_tots$year,
      y_smp        = round(annual_tots$tot_enr * 0.9),
      y_obs        = rbinom(
        nrow(annual_tots),
        round(annual_tots$tot_enr * 0.9),
        probs
      )
    )
  } else {
    sim_school <- data.frame(
      vaxview_type = character(0),
      Year = integer(0),
      y_smp = integer(0),
      y_obs = integer(0)
    )
  }

  # ── Assign names ────────────────────────────────────────────────────────────
  kg_sim <- kg_sim_full
  kg_sim$unit_id <- kg_sim$unit_id + 4L

  kg_sim$county <- county_names[kg_sim$enc_unit_id - 1L]
  kg_sim$school <- school_names[kg_sim$unit_id - 4L]

  kg_sim <- kg_sim |>
    dplyr::select(loc_id = school, parent_id = county, year, enc_unit_id, unit_id,
                  y_obs, y_smp)

  vv_sim <- dplyr::bind_rows(sim_child, sim_school, sim_teen) |>
    dplyr::rename(year = Year) |>
    dplyr::mutate(loc_id = "State")

  # ── Calendar years ──────────────────────────────────────────────────────────
  kg_sim$year <- kg_sim$year + 1995L
  vv_sim$year <- vv_sim$year + 1995L

  # ── Weight / dose / lag-year info ───────────────────────────────────────────
  kg_sim <- kg_sim |>
    dplyr::mutate(ly_min = 5, ly_max = 5, dose = 2, weight = 1)

  vv_sim <- vv_sim |>
    dplyr::mutate(
      ly_min = dplyr::case_when(
        vaxview_type == "school"                  ~ 5,
        vaxview_type == "teen"                    ~ 14,
        vaxview_type == "child" & age == "24 months" ~ 2,
        TRUE                                      ~ 3
      ),
      ly_max = dplyr::case_when(
        vaxview_type == "school"                  ~ 5,
        vaxview_type == "teen"                    ~ 18,
        vaxview_type == "child" & age == "24 months" ~ 2,
        TRUE                                      ~ 3
      ),
      dose = dplyr::case_when(
        vaxview_type %in% c("school", "teen")     ~ 2,
        TRUE                                      ~ 1
      ),
      weight = dplyr::if_else(vaxview_type == "teen", 1 / 5, 1)
    )

  # ── Bind and normalise cohorts ───────────────────────────────────────────────
  observations <- dplyr::bind_rows(kg_sim, vv_sim |> dplyr::mutate(unit_id = 1L)) |>
    dplyr::mutate(
      # imuGAP works in cohort indices, so translate each observation's age window.
      by_max     = year - ly_min,
      by_min     = year - ly_max,
      cohort_min = by_min - min(by_min) + 1L,
      cohort_max = by_max - min(by_min) + 1L
    ) |>
    dplyr::select(-by_min, -by_max) |>
    dplyr::rename(positive = y_obs, sample_n = y_smp) |>
    dplyr::mutate(obs_id = dplyr::row_number())

  observations <- setDT(observations)

  # ── Populations ─────────────────────────────────────────────────────────────
  populations <- rbindlist(lapply(seq_len(nrow(observations)), function(i) {
    # One observation can represent several age/cohort cells; weights split it.
    data.table(
      obs_id = observations$obs_id[i],
      loc_id = observations$loc_id[i],
      cohort = observations$cohort_max[i]:observations$cohort_min[i],
      age    = observations$ly_min[i]:observations$ly_max[i],
      dose   = observations$dose[i],
      weight = observations$weight[i]
    )
  }))

  # ── Locations ────────────────────────────────────────────────────────────────
  locations <- rbindlist(list(
    data.table(loc_id = "State",        parent_id = NA_character_),
    data.table(loc_id = county_names,   parent_id = "State"),
    unique(
      observations[loc_id != "State", .(loc_id, parent_id)]
    )
  ))

  # ── Return ────────────────────────────────────────────────────────────────────
  list(
    observations = observations,
    populations  = populations,
    locations    = locations,
    params = list(
      seed                = seed,
      n_yr                = n_yr,
      n_cohort            = n_cohort,
      lambda              = lambda,
      sigma_sch           = sigma_sch,
      sigma_cnty          = sigma_cnty,
      other_vax_reduction = other_vax_reduction,
      n_schools           = n_schools,
      phi_range           = range(phi_st),
      cov_at_age5_dose2   = cov[5, 2]
    )
  )
}

# =============================================================================
#' Print a readable summary of a simulate_imuGAP_data() result
# =============================================================================
print_sim_summary <- function(sim) {
  obs  <- sim$observations
  pop  <- sim$populations
  locs <- sim$locations
  p    <- sim$params

  cat("══════════════════════════════════════════════\n")
  cat("  imuGAP synthetic data summary\n")
  cat("══════════════════════════════════════════════\n\n")

  cat("── Parameters used ───────────────────────────\n")
  cat(sprintf("  seed              : %d\n",   p$seed))
  cat(sprintf("  sigma_sch         : %.2f\n", p$sigma_sch))
  cat(sprintf("  sigma_cnty        : %.2f\n", p$sigma_cnty))
  cat(sprintf("  lambda            : %.1f, %.1f\n", p$lambda[1], p$lambda[2]))
  cat(sprintf("  other_vax_reduc.  : %.2f\n", p$other_vax_reduction))
  cat(sprintf("  phi_st range      : %.3f – %.3f\n", p$phi_range[1], p$phi_range[2]))
  cat(sprintf("  cov (age 5 dose 2): %.3f\n", p$cov_at_age5_dose2))
  cat(sprintf("  schools per county: %s  (total = %d)\n",
              paste(p$n_schools, collapse = ", "), sum(p$n_schools)))
  cat("\n")

  cat("── Hierarchy ─────────────────────────────────\n")
  cat(sprintf("  locations: %d total  (%d state, %d counties, %d schools)\n",
              nrow(locs),
              sum(is.na(locs$parent_id)),
              sum(!is.na(locs$parent_id) & locs$parent_id == "State"),
              sum(!is.na(locs$parent_id) & locs$parent_id != "State")))
  cat("\n")

  cat("── Observations ──────────────────────────────\n")
  src <- obs[, .(
    n_rows    = .N,
    yr_range  = paste(min(year), max(year), sep = "–"),
    mean_prop = round(mean(positive / sample_n, na.rm = TRUE), 3),
    min_prop  = round(min(positive  / sample_n, na.rm = TRUE), 3),
    max_prop  = round(max(positive  / sample_n, na.rm = TRUE), 3)
  ), by = vaxview_type]

  print(as.data.frame(src), row.names = FALSE)
  cat(sprintf("\n  total obs rows : %d\n", nrow(obs)))
  cat(sprintf("  total pop rows : %d\n", nrow(pop)))
  cat("\n")

  cat("── Sanity checks ─────────────────────────────\n")
  checks <- list(
    "positive <= sample_n always"    = all(obs$positive <= obs$sample_n),
    "no negative counts"             = all(obs$positive >= 0),
    "no zero sample sizes"           = all(obs$sample_n  >  0),
    "obs_id unique"                  = uniqueN(obs$obs_id) == nrow(obs),
    "pop obs_ids all valid"          = all(pop$obs_id %in% obs$obs_id),
    "cohort values positive"         = all(pop$cohort > 0),
    "no NA in observations"          = !anyNA(obs[, .(positive, sample_n, year,
                                                      loc_id, dose, weight)]),
    "school parents are counties"    = {
      school_locs <- locs[!is.na(parent_id) & parent_id != "State"]
      county_locs <- locs[!is.na(parent_id) & parent_id == "State"]
      all(school_locs$parent_id %in% county_locs$loc_id)
    }
  )

  for (nm in names(checks)) {
    flag <- if (isTRUE(checks[[nm]])) "✓" else "✗ FAILED"
    cat(sprintf("  %s  %s\n", flag, nm))
  }

  all_ok <- all(vapply(checks, isTRUE, logical(1)))
  cat("\n")
  if (all_ok) {
    cat("  ✓ All checks passed.\n")
  } else {
    cat("  ✗ Some checks failed — review parameters.\n")
  }
  cat("══════════════════════════════════════════════\n")
  invisible(sim)
}