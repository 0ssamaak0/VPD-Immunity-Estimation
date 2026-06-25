# =============================================================================
# build_data.R  —  regenerates dashboard/data.js from the repo pipeline.
#
# Run from the repo root inside the mmed conda env:
#   conda run -n mmed Rscript dashboard/build_data.R
#
# Output: dashboard/data.js  ->  window.DASHBOARD_DATA = { ... };
#
# All numbers come from the live pipeline (synthetic sim + fitted imuGAP
# posteriors), so the dashboard reflects the CORRECT values, not the stale
# report figures. Coverage profiles are State/County/School *two-dose* (MCV2)
# posteriors — the same quantity that feeds the Reed-Frost phase.
# =============================================================================

suppressMessages({
  library(imuGAP)
  library(data.table)
  library(dplyr)
  library(jsonlite)
})

base_dir <- "."
setwd(base_dir)

COUNTIES <- c("Scruggs", "Simone", "Watson")
AGES     <- 5:18

# ---------------------------------------------------------------------------
# 1. EDA input data  (default synthetic sim == baseline sim1, seed 93254)
# ---------------------------------------------------------------------------
source("raw-data/simulate_imuGAP_data.R")
sim <- simulate_imuGAP_data()
obs <- as.data.table(sim$observations)
obs[, frac := positive / sample_n]

# Per-school kindergarten records (vaxview_type is NA) -> two-dose
schools_dt <- obs[is.na(vaxview_type),
                  .(school = loc_id, county = parent_id,
                    enrollment = sample_n, frac = round(frac, 4))]

# State child stream (MCV1 / dose 1) over calendar time
child_dt <- obs[vaxview_type == "child",
                .(year = year, age = age, frac = round(frac, 4))]

scruggs_schools <- sim$locations[parent_id == "Scruggs", loc_id]

# ---------------------------------------------------------------------------
# 2. imuGAP coverage posteriors  (baseline = sim1, reduced = sim4), dose 2
# ---------------------------------------------------------------------------
summ <- function(path) {
  s <- as.data.table(summary(readRDS(path)))
  s[dose == 2 & age %in% AGES]
}
s1 <- summ("predict_sim1.rds")   # baseline
s4 <- summ("predict_sim4.rds")   # reduced

# State: average over rows per age
state_cov <- function(s) {
  s[loc_id == "State", .(
    q2_5  = round(mean(q2_5),  4),
    q50   = round(mean(q50),   4),
    q97_5 = round(mean(q97_5), 4),
    mean  = round(mean(mean),  4)
  ), by = age][order(age)]
}

# County: average over rows per (county, age)
county_cov <- function(s) {
  d <- s[loc_id %in% COUNTIES, .(
    q2_5  = round(mean(q2_5),  4),
    q50   = round(mean(q50),   4),
    q97_5 = round(mean(q97_5), 4)
  ), by = .(loc_id, age)][order(loc_id, age)]
  split(d[, .(age, q2_5, q50, q97_5)], d$loc_id)
}

# School boxplot in Scruggs: per-school median (q50) by age
school_cov <- function(s) {
  s[loc_id %in% scruggs_schools,
    .(school = loc_id, age = age, q50 = round(q50, 4))][order(age, school)]
}

coverage <- list(
  state  = list(baseline = state_cov(s1),  reduced = state_cov(s4)),
  county = list(baseline = county_cov(s1), reduced = county_cov(s4)),
  school = list(baseline = school_cov(s1), reduced = school_cov(s4))
)

# ---------------------------------------------------------------------------
# 3. Reed-Frost set-up: immunity coverage timeline + per-window immune seed
#    ave = 30 * c(baseline_by_age, reduced_by_age)   (length 28)
#    window width 14  ->  15 rolling windows           (matches histogram.R)
# ---------------------------------------------------------------------------
average1 <- readRDS("average_1.rds")
average4 <- readRDS("average_4.rds")
ave <- c(30 * average1$mean_coverage, 30 * average4$mean_coverage)

N            <- 14 * 30   # 420
WINDOW_WIDTH <- 14
N_WINDOWS    <- length(ave) - WINDOW_WIDTH + 1   # 15

windows <- lapply(seq_len(N_WINDOWS), function(t) {
  idx <- t:(t + WINDOW_WIDTH - 1)
  R0  <- round(sum(ave[idx]))
  R0  <- min(R0, N - 2)
  list(window = t, R0_seed = R0, S0 = N - R0 - 1)
})

# ---------------------------------------------------------------------------
# 4. Validation reference: per-window median final outbreak (R engine, p=0.04)
#    (printed only — used to check the JS port; NOT baked into data.js)
# ---------------------------------------------------------------------------
source("Phase_2_Dynamical_part/Reed_Frost_model.R")
set.seed(123)
ref <- sapply(windows, function(w) {
  finals <- replicate(1000, {
    rf <- simulate_reed_frost_seir(20, c(S = w$S0, E = 1, I = 0, R = w$R0_seed), 0.04)
    tail(rf$E + rf$I + rf$R - rf$R[1], 1)
  })
  c(R0 = w$R0_seed, median = median(finals), p_major = mean(finals >= 5))
})
cat("\n--- R reference (p=0.04, 1000 sims/window) ---\n")
print(round(t(ref), 2))

# ---------------------------------------------------------------------------
# 5. Assemble + write data.js
# ---------------------------------------------------------------------------
DATA <- list(
  meta = list(
    N = N, p_default = 0.04, generations = 20, E0 = 1,
    window_width = WINDOW_WIDTH, n_windows = N_WINDOWS, herd = 0.95,
    counties = COUNTIES, ages = AGES,
    generated = as.character(Sys.time())
  ),
  ave      = round(ave, 3),
  windows  = windows,
  eda      = list(child = child_dt, schools = schools_dt),
  coverage = coverage
)

json <- toJSON(DATA, dataframe = "rows", auto_unbox = TRUE, digits = 6, na = "null")
writeLines(paste0("window.DASHBOARD_DATA = ", json, ";\n"),
           "dashboard/data.js")

cat("\nWrote dashboard/data.js\n")
cat(sprintf("  eda.schools rows : %d\n", nrow(schools_dt)))
cat(sprintf("  eda.child rows   : %d\n", nrow(child_dt)))
cat(sprintf("  windows          : %d (R0 seed %d -> %d)\n",
            N_WINDOWS, windows[[1]]$R0_seed, windows[[N_WINDOWS]]$R0_seed))
cat(sprintf("  state baseline age18 q50 : %.3f | reduced : %.3f\n",
            coverage$state$baseline[age == 18, q50],
            coverage$state$reduced[age == 18, q50]))
