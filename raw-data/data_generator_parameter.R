source("simulate_imuGAP_data.R")

# ── Default (matches the package fixture exactly) ──────────────
sim <- simulate_imuGAP_data()
print_sim_summary(sim)

# ── Try higher school-level noise ──────────────────────────────
sim2 <- simulate_imuGAP_data(sigma_sch = 1.5, sigma_cnty = 0.6)
print_sim_summary(sim2)

# ── Fewer schools, different seed ──────────────────────────────
sim3 <- simulate_imuGAP_data(seed = 42, n_schools = c(4, 4, 4))
print_sim_summary(sim3)

# ── Lower overall coverage (shift phi_st down) ─────────────────
low_phi <- sim$params  # just to see the default, then modify
sim4 <- simulate_imuGAP_data(phi_st = rep(0.75, 33))
print_sim_summary(sim4)

# ── Access the tables directly ─────────────────────────────────
sim$observations
sim$populations
sim$locations
