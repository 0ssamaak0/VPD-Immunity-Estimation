library(tidyverse)

#####################################################################
##  Reed-Frost SEIR Model 
#####################################################################
## S[t]   -> Susceptible population
## E[t]   -> Newly exposed (infected but not yet infectious)
## I[t]   -> Active infectious individuals (equal to E[t-1])
## R[t]   -> Removed/Recovered individuals (accumulated from I[t-1])


simulate_reed_frost_seir <- function(generations, y0, p) {
  # Reed-Frost is stochastic: the same inputs can produce different outbreak sizes.
  # Initialize vectors to store states over time
  S <- numeric(generations + 1)
  E <- numeric(generations + 1)
  I <- numeric(generations + 1)
  R <- numeric(generations + 1)
  
  # Set initial states at generation 1 (t = 0)
  S[1] <- y0["S"]
  E[1] <- y0["E"]
  I[1] <- y0["I"]
  R[1] <- y0["R"]
  
  for (t in 1:generations) {
    # If there are no infectious individuals, the epidemic stops expanding
    if (I[t] == 0) {
      S[t+1] <- S[t]
      E[t+1] <- 0
      I[t+1] <- E[t]       # Latent individuals progress to infectious
      R[t+1] <- R[t] + I[t] # Infectious individuals recover
      next
    }
    
    # Probability of a susceptible getting infected by the current pool of I
    # p is the per-pair effective contact probability for one generation.
    prob_infection <- 1 - (1 - p)^I[t]
    
    # Sample the number of new exposures using the binomial distribution
    new_exposures <- rbinom(1, size = S[t], prob = prob_infection)
    
    # State transitions for the next time step
    S[t+1] <- S[t] - new_exposures
    E[t+1] <- new_exposures
    I[t+1] <- E[t]       # Previous generation's E becomes active I
    R[t+1] <- R[t] + I[t] # Previous generation's I recovers
  }
  
  # Compile into a clean data frame
  results <- data.frame(
    Generation = 0:generations,
    S = S,
    E = E,
    I = I,
    R = R
  )
  
  return(results)
}

## Setup
pop <- 100               # Total population size
generations_count <- 25  # Number of steps/generations to track
p_contact <- 0.04        # Probability of effective contact per pair

# Initial state: 1 individual starts Exposed (E = 1), 0 active Infectious
y0_rf <- c(S = pop - 1, E = 1, I = 0, R = 0)

# discrete simulation
rf_data <- simulate_reed_frost_seir(generations_count, y0_rf, p_contact)

## Plot
rf_long <- rf_data |> 
  pivot_longer(cols = c(S, E, I, R), names_to = "Compartment", values_to = "count") |> 
  mutate(Compartment = factor(Compartment, levels = c('S', 'E', 'I', 'R')))

ggplot(rf_long, aes(x = Generation, y = count, color = Compartment)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Reed-Frost SEIR Chain-Binomial Simulation", 
    y = "Count", 
    x = "Generation Step"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("S" = "#377eb8", "E" = "#ff7f00", "I" = "#e41a1c", "R" = "#4daf4a"))

