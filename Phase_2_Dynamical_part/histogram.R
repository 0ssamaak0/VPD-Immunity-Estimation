library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)


base_dir <- "."
setwd(base_dir)
source("./Phase_2_Dynamical_part/Reed_Frost_model.R")

# Load baseline coverage data
# Coverage summaries become the starting immune population for each window.
average1 <- readRDS("average_1.rds")
average4 <- readRDS("average_4.rds")

# Scale by cohort size (30 students per grade) over 14 grades
ave1 <- 30 * average1$mean_coverage   # length 14
ave4 <- 30 * average4$mean_coverage   # length 14

ave  <- c(ave1, ave4)                 # length 28

generations_count <- 20
total_pop         <- 14 * 30          # 420
p_contact         <- 0.04
window_width      <- 14

number_of_simulation <- 1000
number_of_windows    <- length(ave) - window_width + 1

add_cumulative_infections <- function(rf_data) {
  rf_data |>
    mutate(
      Cumulative_Infections = .data$E + .data$I + .data$R - first(.data$R)
    )
}

all_rf <- vector("list", number_of_windows)

for (t in seq_len(number_of_windows)) {
  idx <- t:(t + window_width - 1)
  R0  <- round(sum(ave[idx]))
  
  R0 <- min(R0, total_pop - 2)
  
  # Initialize compartments based on background immunity (R0)
  y0_rf <- c(S = total_pop - R0 - 1, E = 1, I = 0, R = R0)

  step_t_rf <- vector("list", number_of_simulation)

  # Repeat the stochastic Reed-Frost run so each window has an outbreak distribution.
  for (i in seq_len(number_of_simulation)) {
    rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p = p_contact)
    
    # Count everyone infected so far, including active exposed/infectious cases.
    rf_data <- add_cumulative_infections(rf_data)
    step_t_rf[[i]] <- rf_data["Cumulative_Infections"]
  }
  all_rf[[t]] <- step_t_rf
}


# function that gets the list of final values
list_for <- function(t){
  final_values <- rep(1, length.out = number_of_simulation)
  for (j in seq_len(number_of_simulation)){
    final_values[j] <- tail(all_rf[[t]][[j]]$Cumulative_Infections, n = 1)
  }
  return(unlist(final_values))
}

# function to plot
plot_histogram_for <- function(t) {
  df <- data.frame(value = list_for(t))
  ggplot(df, aes(x = .data$value)) + 
    geom_histogram(binwidth = 1, fill = "#377eb8", color = "white", linewidth = 0.25) +
    labs(
      title = paste("Cumulative infections for rolling window", t),
      subtitle = paste(number_of_simulation, "stochastic Reed-Frost simulations"),
      x = "Cumulative infections",
      y = "Simulation count"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title       = element_text(face = "bold", size = 20),
      plot.subtitle    = element_text(color = "gray35", size = 14),
      axis.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

#change the argument to plot the histogram (argument between 1 and 15)

# How are we going to plot representative results for the 15 rolling windows
plot_histogram_for(1)
plot_histogram_for(4)
plot_histogram_for(7)
plot_histogram_for(10)
plot_histogram_for(13)

DT <- data.frame()

for(i in seq_len(number_of_windows)){
    dt <- data.frame(year = i,
    cumulative_infection = list_for(i))
    DT <- bind_rows(DT, dt) 
}

selected_years <- tail(seq_len(number_of_windows), 5)

# Show the final five rolling windows instead of plotting all panels.
histogram_plot <- ggplot(
  DT |> filter(year %in% selected_years),
  aes(x = cumulative_infection)
) +
  geom_histogram(
    binwidth = 2,
    fill     = "#377eb8",
    color    = "white",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~ year,
    nrow     = 2,
    ncol     = 3,
    scales   = "fixed",
    labeller = labeller(year = function(x) paste("Rolling window", x))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Distribution of cumulative infections across simulation runs",
    subtitle = paste0(
      "Reed-Frost SEIR · N = ", total_pop,
      " · p = ", p_contact,
      " · ", number_of_simulation, " simulations per window",
      " · final 5 rolling windows"
    ),
    x = "Cumulative infections",
    y = "Count"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title         = element_text(face = "bold", size = 24, hjust = 0.5),
    plot.subtitle      = element_text(size = 14, color = "gray35", hjust = 0.5,
                                      margin = margin(b = 8)),
    strip.text         = element_text(face = "bold", size = 15),
    strip.background   = element_rect(fill = "gray92", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    panel.spacing.x    = unit(0.8, "lines"),
    panel.spacing.y    = unit(0.7, "lines"),
    axis.title         = element_text(face = "bold", size = 14),
    axis.text          = element_text(size = 12),
    plot.margin        = margin(8, 12, 8, 10)
  )

print(histogram_plot)

ggsave(
  "./Phase_2_Dynamical_part/cumulative_infection_histograms.png",
  histogram_plot,
  width  = 14,
  height = 8,
  dpi    = 300,
  bg     = "white"
)
