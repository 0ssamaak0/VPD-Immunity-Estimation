library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)

base_dir <- "."
setwd(base_dir)
source("./Phase_2_Dynamical_part/Reed_Frost_model.R")

# Load baseline coverage data
# These RDS files are the age-level coverage summaries from the imuGAP scenarios.
average1 <- readRDS("average_1.rds")
average4 <- readRDS("average_4.rds")

# Scale by cohort size (30 students per grade) over 14 grades
ave1 <- 30 * average1$mean_coverage   # length 14
ave4 <- 30 * average4$mean_coverage   # length 14
ave  <- c(ave1, ave4)                 # length 28

generations_count <- 20
total_pop         <- 14 * 30          # 420
p_contact         <- 0.04

all_rf <- vector("list", 17) # Explicitly sized for the 17 rolling windows

for (t in 1:17) {
  # Each window turns a slice of coverage history into an initial immune count.
  idx <- t:(t + 17)
  idx <- idx[idx <= length(ave)]
  R0  <- round(sum(ave[idx]))
  
  R0 <- min(R0, total_pop - 2)
  
  # Initialize compartments based on background immunity (R0)
  y0_rf <- c(S = total_pop - R0 - 1, E = 1, I = 0, R = R0)
  
  rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p = p_contact)
  
  # Calculate cumulative infections: current R minus initial immune cluster
  rf_data <- rf_data |> 
    mutate(Cumulative_Infections = R - first(R))
  
  all_rf[[t]] <- rf_data |>
    pivot_longer(
      cols      = c(S, E, I, R, Cumulative_Infections),
      names_to  = "Compartment",
      values_to = "count"
    ) |>
    mutate(
      Compartment = factor(Compartment, levels = c("S", "E", "I", "R", "Cumulative_Infections")),
      run         = factor(t)
    )
}

# Combine all 17 runs into one master data frame
all_rf_df <- bind_rows(all_rf) |>
  mutate(
    Compartment = recode(
      Compartment,
      S = "Susceptible",
      E = "Exposed",
      I = "Infectious",
      R = "Recovered",
      Cumulative_Infections = "Cumulative Infections"
    ),
    Compartment = factor(
      Compartment,
      levels = c("Susceptible", "Exposed", "Infectious", "Recovered", "Cumulative Infections")
    )
  )

# Calculate the median behavior across the 17 runs
# The median summarizes the typical trajectory while keeping individual runs visible.
median_rf_df <- all_rf_df |>
  group_by(Generation, Compartment) |>
  summarise(median_count = median(count), .groups = "drop")


# GENERATE THE LATEX SUMMARY (I need it for the rapport)
simulation_summary <- data.frame(
  Data_Object = c("ave", "rf_data (last run)", "all_rf_df", "median_rf_df"),
  Description = c(
    "Combined historical immunity coverage vector (Sim 1 + Sim 4)",
    "Single Reed-Frost simulation run output with cumulative metrics",
    "Complete long-format dataset combining all 17 rolling runs",
    "Aggregated dataset tracking the median behavior across timelines"
  ),
  Key_Variables = c(
    "Numeric vector (coverage values)",
    "Generation, S, E, I, R, Cumulative_Infections",
    "Generation, Compartment, count, run",
    "Generation, Compartment, median_count"
  ),
  Dimensions = c(
    paste(length(ave), "values"),
    paste(nrow(rf_data), "rows"),
    paste(format(nrow(all_rf_df), big.mark=","), "rows"),
    paste(nrow(median_rf_df), "rows")
  )
)

# Print the table formatted perfectly for LaTeX booktabs
cat("\n--- Copy this LaTeX Table Code for your Report ---\n\n")
kable(simulation_summary, 
      format = "latex", 
      booktabs = TRUE, 
      col.names = c("Data Object", "Description", "Key Variables", "Dimensions / Rows"))

# =============================================================================
# PLOT GENERATION (Shows all 17 simulation timelines)
# =============================================================================

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
  facet_wrap(~ Compartment, scales = "free_y", ncol = 2) + 
  guides(
    color = guide_legend(override.aes = list(linewidth = 2, alpha = 1))
  ) +
  labs(
    title    = "Reed-Frost SEIR simulations across rolling immunity windows",
    subtitle = paste0(
      "17 simulations (rolling 18 successive years coverage windows) · ",
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

# Display plot
print(plot_simulation)
#Histogram 
all_rf_df |> subset(Compartment == "Cumulative Infections") |> group_by(run) |> summarize(max_value = max(count)) |> ggplot() + aes(x = max_value) + geom_histogram() + theme_minimal()
