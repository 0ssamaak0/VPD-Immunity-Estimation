library(ggplot2)
library(tidyr)
library(dplyr)

base_dir <- file.path("~/Downloads/MMEDGit/VPD-Immunity-Estimation")
setwd(base_dir)
source("./Phase_2_Dynamical_part/Reed_Frost_model.R")

# Fix 1: extract column first, then multiply
average1 <- readRDS("average_1.rds")
average4 <- readRDS("average_4.rds")

ave1 <- 30 * average1$mean_coverage   # length 14
ave4 <- 30 * average4$mean_coverage   # length 14
ave  <- c(ave1, ave4)                 # length 28

generations_count <- 20
total_pop         <- 14 * 30          # 420

all_rf <- vector("list", 18)

for (t in 1:17) {
  # Fix 2: [ not () for indexing; guard against going past length(ave)
  idx <- t:(t + 17)
  idx <- idx[idx <= length(ave)]
  R0  <- round(sum(ave[idx]))

  # Guard: R0 cannot exceed total_pop - 2 (need at least S=1 and E=1)
  R0 <- min(R0, total_pop - 2)

  y0_rf <- c(S = total_pop - R0 - 1, E = 1, I = 0, R = R0)

  rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p = 0.04)

  # Fix 3: collect into list — Generation column already exists from the function
  all_rf[[t]] <- rf_data |>
    pivot_longer(
      cols      = c(S, E, I, R),
      names_to  = "Compartment",
      values_to = "count"
    ) |>
    mutate(
      Compartment = factor(Compartment, levels = c("S", "E", "I", "R")),
      run         = factor(t)
    )
}

all_rf_df <- bind_rows(all_rf)

# Fix 4: one ggplot over the combined data — group keeps runs separate
plot_simulation <- ggplot(
  all_rf_df,
  aes(x     = Generation,
      y     = count,
      color = Compartment,
      group = interaction(Compartment, run))
) +
  geom_line(linewidth = 1, alpha = 0.4) +
  theme_minimal(base_size = 14) +
  scale_color_manual(
    values = c("S" = "#377eb8", "E" = "#ff7f00",
               "I" = "#e41a1c", "R" = "#4daf4a")
  ) +
  labs(
    title = "Reed-Frost SEIR Chain-Binomial Simulation",
    y     = "Count",
    x     = "Generation"
  )

plot_simulation