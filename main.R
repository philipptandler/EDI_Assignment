library(here)
library(dplyr)
# Set working directory only if not already set
setwd(here::here())

# load configuration
source("config.R")

# PRE: dataframe df with corridor | phase | habitat and levels according to config.R
# POST: returns dataframe n_df with columns n00 | n01 | n10 | n11,
# corresponding to detections of wolf and caribou
simulate_daily_observations <- function(df_obs){
  n00 <- rep(0,n_treats)
  n01 <- rep(0,n_treats)
  n10 <- rep(0,n_treats)
  n11 <- rep(0,n_treats)
  # fill dataframe for the predictors
  for (row in 1:n_treats){
    effect_corridor <- as.numeric(observation[row,]$corridor)-1
    effect_phase <- as.numeric(observation[row,]$phase)-1
    effect_habitat <- as.numeric(observation[row,]$habitat)-1
    # first for caribou (not dependent on wolf)
    linear_predictor_caribou <- 
      gamma_0 +
      gamma_corr * effect_corridor +
      gamma_phas * effect_phase +
      gamma_habi * effect_habitat +
      gamma_corr_phas * (effect_corridor*effect_phase) + 
      rnorm(1, mean = 0, sd = gamma_sigma)
    caribou <- rbinom(n = n_observation_bins, size = 1,
                      prob = plogis(linear_predictor_caribou))
    # now calculate wolf conditional on caribou detection
    linear_predictor_wolf <- 
      beta_0 +
      beta_corr * effect_corridor +
      beta_phas * effect_phase +
      beta_cari * caribou +
      beta_habi * effect_habitat +
      beta_corr_phas * (effect_corridor*effect_phase) +
      beta_corr_cari * (effect_corridor*caribou) + 
      beta_cari_phas * (caribou*effect_phase) +
      beta_cari_corr_phas * (caribou*effect_corridor*effect_phase) +
      rnorm(1, mean = 0, sd = beta_sigma)
    wolf <- rbinom(n = length(linear_predictor_wolf), size = 1,
      prob = plogis(linear_predictor_wolf))
    
    # calculate relative frequency:
    tab <- table(
      factor(caribou, levels = 0:1),
      factor(wolf, levels = 0:1)
    )
    n00[row] <- tab["0", "0"]
    n01[row] <- tab["0", "1"]
    n10[row] <- tab["1", "0"]
    n11[row] <- tab["1", "1"]
  }
  df_n <- data.frame(
    n00 = n00,
    n01 = n01,
    n10 = n10,
    n11 = n11
  )
}

# PRE: n_cameras (per treatment) >= 1
# POST: returns observation dataframe for design in configuration file,
# cameraID | corridor | sectionID | habitat | phase | n00 | n01 | n10 | n11
simulate_observations <- function(n_cameras){
  n_treats <- n_cameras*length(corridor)*n_sections*length(habitat)*length(phase)
  full_observations <- data.frame(
    cameraID = as.factor(
      rep(paste0("C", 1:(n_treats/length(phase))),
          each = length(phase))
      ),
    corridor = relevel(as.factor(
      rep(corridor, 
          each = n_treats/length(corridor))),
      ref = corridor[2]
    ),
    sectionID = as.factor(
      rep(paste0("S", 1:(n_sections*length(corridor))), 
          each = n_treats/(length(corridor)*n_sections))
    ),
    habitat = relevel(as.factor(
      rep(habitat, 
          each = n_cameras*length(phase),
          times = n_treats/(n_cameras*length(phase)*length(habitat)))),
      ref = habitat[2]
    ),
    phase = relevel(as.factor(
      rep(phase,
          times = n_treats/length(phase))),
      ref = phase[2]
    )
  )
  df_n <- simulate_daily_observations(full_observations)
  full_observations <- cbind(
    full_observations, df_n
  )
  full_observations
}

# PRE: n_cameras (per treatment) >= 1, hypothesis = c('A', 'B', 'C')
# POST: returns data frame with 
simulate_data <- function(n_cameras, hypothesis){
  df_full <- simulate_observations(n_cameras)
  # aggregate depending on hypothesis:
  
}

set.seed(1234)



