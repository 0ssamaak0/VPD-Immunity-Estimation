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
average1$mean_coverage

ave  <- c(ave1, ave4)                 # length 28

generations_count <- 20
total_pop         <- 14 * 30          # 420
p_contact         <- 0.02

all_rf <- vector("list", 14) # Explicitly sized for the 17 rolling windows

for (t in 1:14) {
  idx <- t:(t + 14)
  idx <- idx[idx <= length(ave)]
  R0  <- round(sum(ave[idx]))
  
  R0 <- min(R0, total_pop - 2)
  
  # Initialize compartments based on background immunity (R0)
  y0_rf <- c(S = total_pop - R0 - 1, E = 1, I = 0, R = R0)
  
  number_of_simulation <-100

  step_t_rf <- vector("list",number_of_simulation)

   # Repeat the stochastic Reed-Frost run so each window has an outbreak distribution.
   for (i in 1:number_of_simulation){
  rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p = p_contact)
  
  # Calculate cumulative infections: current R minus initial immune cluster
  rf_data <- rf_data |> 
    mutate(Cumulative_Infections = R - first(R))
    step_t_rf[[i]] <- rf_data["Cumulative_Infections"]
  all_rf[[t]] <- step_t_rf
    }
}


#function that get the list of final value
list_for <- function(t){
  final_values <-rep(1,length.out = 100)
   for (j in 1:100){
    final_values[j]<-tail(all_rf[[t]][[j]],n= 1)[1]}
  return(unlist(final_values))
}

list_for(1)
#function lpot
plot_histogram_for <-function(t)
{   df <- data.frame(value = list_for(t))
    ggplot(df, aes(x = .data$value)) + 
  geom_histogram(binwidth = 1, fill = "#377eb8", color = "white", linewidth = 0.25) +
  labs(
    title = paste("Cumulative infections for rolling window", t),
    subtitle = "100 stochastic Reed-Frost simulations",
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

#change the argument to plot the hoistogram (argument between 1 and 13)

#How are we going to plot the result for all the 13 years
plot_histogram_for(1)
plot_histogram_for(4)
plot_histogram_for(7)
plot_histogram_for(10)
plot_histogram_for(13)

DT <- data.frame()

for(i in 1:14){
    dt<-data.frame(year=i,
    cumulative_infection = list_for(i))
    DT <- bind_rows(DT,dt) 
}
head(DT)

selected_years <- as.integer(round(seq(1, 13, length.out = 6)))

# Show representative rolling windows instead of plotting all 13 panels.
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
    scales   = "free_y",
    labeller = labeller(year = function(x) paste("Rolling window", x))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Distribution of cumulative infections across simulation runs",
    subtitle = paste0(
      "Reed-Frost SEIR · N = ", total_pop,
      " · p = ", p_contact,
      " · 100 simulations per window"
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