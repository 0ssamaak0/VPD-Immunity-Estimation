library(ggplot2)
base_dir <- file.path("~/Downloads/MMEDGit/VPD-Immunity-Estimation")
setwd(base_dir)
source("./dynamical_model/Reed_frost.R")

pop_age1 <- rep(30, length.out = 14)

average1<- readRDS("average_1.rds")

ave1 <- (pop_age1 * average1)$mean_coverage

pop_age4 <- rep(30, length.out = 14)

average4<- readRDS("average_4.rds")

ave4 <- (pop_age4 * average4)$mean_coverage

ave <- c(ave1, ave4)

ave


plot_simulation <-ggplot()
for (t in 1:18) {
   R0<-sum(ave(t:(t+17)))
   y0_rf <- c(S = 14*30-R0-1, E = 1, I = 0, R = R0)
   p_contact <- 1/2
   rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p_contact)
   rf_long <- rf_data |> 
   pivot_longer(cols = c(S, E, I, R), names_to = "Compartment", values_to = "count") |> 
   mutate(Compartment = factor(Compartment, levels = c('S', 'E', 'I', 'R')))

   plot_simulation <- plot_simulation + ggplot(rf_long, aes(x = Generation, y = count, color = Compartment)) +
   geom_line(linewidth = 1.2) +
   geom_point(size = 2) +
   theme_minimal(base_size = 14) +
   scale_color_manual(values = c("S" = "#377eb8", "E" = "#ff7f00", "I" = "#e41a1c", "R" = "#4daf4a"))
}

plot_simulation <- plot_simulation + labs(
    title = "Reed-Frost SEIR Chain-Binomial Simulation", 
    y = "Count", 
    x = "Generation Step"
  )