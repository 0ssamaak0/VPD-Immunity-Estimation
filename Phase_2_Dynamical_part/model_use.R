library(ggplot2)
library(tidyr)
library(dplyr)

base_dir <- file.path("~/MMED/VPD-Immunity-Estimation")
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
p_contact         <- 0.04

all_rf <- vector("list", 18)

for (t in 1:17) {
  # Fix 2: [ not () for indexing; guard against going past length(ave)
  idx <- t:(t + 17)
  idx <- idx[idx <= length(ave)]
  R0  <- round(sum(ave[idx]))

  # Guard: R0 cannot exceed total_pop - 2 (need at least S=1 and E=1)
  R0 <- min(R0, total_pop - 2)

  y0_rf <- c(S = total_pop - R0 - 1, E = 1, I = 0, R = R0)

  rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p = p_contact)

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

all_rf_df <- bind_rows(all_rf) |>
  mutate(
    Compartment = recode(
      Compartment,
      S = "Susceptible",
      E = "Exposed",
      I = "Infectious",
      R = "Recovered"
    ),
    Compartment = factor(
      Compartment,
      levels = c("Susceptible", "Exposed", "Infectious", "Recovered")
    )
  )

median_rf_df <- all_rf_df |>
  group_by(Generation, Compartment) |>
  summarise(median_count = median(count), .groups = "drop")

compartment_colors <- c(
  "Susceptible" = "#377eb8",
  "Exposed"     = "#ff7f00",
  "Infectious"  = "#e41a1c",
  "Recovered"   = "#4daf4a"
)

plot_simulation <- ggplot(all_rf_df) +
  geom_line(
    aes(x = Generation, y = count, group = run),
    color = "gray75",
    linewidth = 0.45,
    alpha = 0.7
  ) +
  geom_line(
    data = median_rf_df,
    aes(x = Generation, y = median_count, color = Compartment),
    linewidth = 1.4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    name   = "Compartment",
    values = compartment_colors
  ) +
  scale_x_continuous(breaks = seq(0, generations_count, by = 5)) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Compartment, scales = "free_y", ncol = 2) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2, alpha = 1)
    )
  ) +
  labs(
    title    = "Reed-Frost SEIR simulations across rolling immunity windows",
    subtitle = paste0(
      "17 simulations (rolling 18-month coverage windows) · ",
      "N = ", format(total_pop, big.mark = ","),
      " · p = ", p_contact,
      " · ", generations_count, " generations"
    ),
    x        = "Generation (time step)",
    y        = "Individuals",
    caption  = "Gray lines: individual runs · Colored lines: median across runs"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "gray40", size = 11),
    plot.caption       = element_text(color = "gray50", size = 9, hjust = 0),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(face = "bold")
  )

plot_simulation














library(ggplot2)
library(tidyr)
library(dplyr)

base_dir <- file.path("~/MMED/VPD-Immunity-Estimation")
setwd(base_dir)
source("./Phase_2_Dynamical_part/Reed_Frost_model.R")

average1 <- readRDS("average_1.rds")
average4 <- readRDS("average_4.rds")

ave1 <- 30 * average1$mean_coverage   # length 14
ave4 <- 30 * average4$mean_coverage   # length 14
ave  <- c(ave1, ave4)                 # length 28

generations_count <- 20
total_pop         <- 14 * 30          # 420
p_contact         <- 0.04

all_rf <- vector("list", 18)

for (t in 1:17) {
  idx <- t:(t + 17)
  idx <- idx[idx <= length(ave)]
  R0  <- round(sum(ave[idx]))
  
  R0 <- min(R0, total_pop - 2)
  
  y0_rf <- c(S = total_pop - R0 - 1, E = 1, I = 0, R = R0)
  
  rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p = p_contact)
  
  # --- CRITICAL ADDITION: Calculate R(t) - R(0) ---
  rf_data <- rf_data |> 
    mutate(Cumulative_Infections = R - first(R))
  
  all_rf[[t]] <- rf_data |>
    pivot_longer(
      cols      = c(S, E, I, R, Cumulative_Infections), # Added here
      names_to  = "Compartment",
      values_to = "count"
    ) |>
    mutate(
      Compartment = factor(Compartment, levels = c("S", "E", "I", "R", "Cumulative_Infections")),
      run         = factor(t)
    )
}

all_rf_df <- bind_rows(all_rf) |>
  mutate(
    Compartment = recode(
      Compartment,
      S = "Susceptible",
      E = "Exposed",
      I = "Infectious",
      R = "Recovered",
      Cumulative_Infections = "Cumulative Infections" # Added label
    ),
    Compartment = factor(
      Compartment,
      levels = c("Susceptible", "Exposed", "Infectious", "Recovered", "Cumulative Infections")
    )
  )

median_rf_df <- all_rf_df |>
  group_by(Generation, Compartment) |>
  summarise(median_count = median(count), .groups = "drop")

# Added a color (purple) for your new trend
compartment_colors <- c(
  "Susceptible"           = "#377eb8",
  "Exposed"               = "#ff7f00",
  "Infectious"            = "#e41a1c",
  "Recovered"             = "#4daf4a",
  "Cumulative Infections" = "#984ea3" 
)

plot_simulation <- ggplot(all_rf_df) +
  geom_line(
    aes(x = Generation, y = count, group = run),
    color = "gray75",
    linewidth = 0.45,
    alpha = 0.7
  ) +
  geom_line(
    data = median_rf_df,
    aes(x = Generation, y = median_count, color = Compartment),
    linewidth = 1.4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    name   = "Compartment",
    values = compartment_colors
  ) +
  scale_x_continuous(breaks = seq(0, generations_count, by = 5)) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Compartment, scales = "free_y", ncol = 2) + # Automatically updates to show 5 panels
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2, alpha = 1)
    )
  ) +
  labs(
    title    = "Reed-Frost SEIR simulations across rolling immunity windows",
    subtitle = paste0(
      "17 simulations (rolling 18-month coverage windows) · ",
      "N = ", format(total_pop, big.mark = ","),
      " · p = ", p_contact,
      " · ", generations_count, " generations"
    ),
    x        = "Generation (time step)",
    y        = "Individuals",
    caption  = "Gray lines: individual runs · Colored lines: median across runs"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "gray40", size = 11),
    plot.caption       = element_text(color = "gray50", size = 9, hjust = 0),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(face = "bold")
  )

plot_simulation

#We want to know I= R0- R(t) for each R(t) generated
