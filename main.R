library(here)
library(dplyr)
library(lme4)
# Set working directory only if not already set
setwd(here::here())

# load configuration
source("config.R")

# PRE: dataframe df_obsyear with corridor | phase | habitat and levels according to config.R,
# n_obsyears corresponding to df_obsyear
# POST: returns dataframe n_df with columns n00 | n01 | n10 | n11,
# corresponding to detections of wolf and caribou
simulate_semidaily_observations <- function(df_obsyear, n_obsyears){
  n00 <- rep(0,n_obsyears)
  n01 <- rep(0,n_obsyears)
  n10 <- rep(0,n_obsyears)
  n11 <- rep(0,n_obsyears)
  # calculate mixed effects for wolves
  beta_effect_per_sectionID <- data.frame(
    sectionIDs = unique(df_obsyear$sectionID),
    mixed_per_section = rnorm(length(unique(df_obsyear$sectionID)),
                              mean = 0,
                              sd = beta_sigma_sectionID)
    )
  beta_effect_per_cameraID <- data.frame(
    cameraIDs = unique(df_obsyear$cameraID),
    mixed_per_camera = rnorm(length(unique(df_obsyear$cameraID)),
                              mean = 0,
                              sd = beta_sigma_cameraID)
  )
  # calculate mixed effects for caribou
  gamma_effect_per_sectionID <- data.frame(
    sectionIDs = unique(df_obsyear$sectionID),
    mixed_per_section = rnorm(length(unique(df_obsyear$sectionID)),
                              mean = 0,
                              sd = gamma_sigma_sectionID)
  )
  gamma_effect_per_cameraID <- data.frame(
    cameraIDs = unique(df_obsyear$cameraID),
    mixed_per_camera = rnorm(length(unique(df_obsyear$cameraID)),
                             mean = 0,
                             sd = gamma_sigma_cameraID)
  )
  # fill dataframe for the predictors
  for (row in 1:n_obsyears){
    effect_corridor <- as.numeric(df_obsyear[row,]$corridor)-1
    effect_phase <- as.numeric(df_obsyear[row,]$phase)-1
    effect_habitat <- as.numeric(df_obsyear[row,]$habitat)-1
    
    # filter mixed effect caribou
    gamma_mixed_effect_section <- gamma_effect_per_sectionID$mixed_per_section[
      gamma_effect_per_sectionID$sectionIDs == df_obsyear[row,]$sectionID
      ]
    gamma_mixed_effect_camera <- gamma_effect_per_cameraID$mixed_per_camera[
      gamma_effect_per_cameraID$cameraIDs == df_obsyear[row,]$cameraID
    ]
    # first for caribou (not dependent on wolf)
    linear_predictor_caribou <- 
      gamma_0 +
      gamma_corr * effect_corridor +
      gamma_phas * effect_phase +
      gamma_habi * effect_habitat +
      gamma_corr_phas * (effect_corridor*effect_phase) + 
      gamma_mixed_effect_section + # mixed effect section
      gamma_mixed_effect_camera + # mixed effect camera
      rnorm(1, mean = 0, sd = gamma_sigma) # residual
    caribou <- rbinom(n = n_semidaily_observations, size = 1,
                      prob = plogis(linear_predictor_caribou))
    # filter mixed effect wolf
    beta_mixed_effect_section <- beta_effect_per_sectionID$mixed_per_section[
      beta_effect_per_sectionID$sectionIDs == df_obsyear[row,]$sectionID
    ]
    beta_mixed_effect_camera <- beta_effect_per_cameraID$mixed_per_camera[
      beta_effect_per_cameraID$cameraIDs == df_obsyear[row,]$cameraID
    ]
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
      beta_mixed_effect_section + # mixed effect section
      beta_mixed_effect_camera + # mixed effect camera
      rnorm(1, mean = 0, sd = beta_sigma) # residual
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

# PRE: n_replicates (per treatment) >= 1
# POST: returns observation dataframe for design in configuration file,
# cameraID | corridor | sectionID | habitat | phase | n00 | n01 | n10 | n11
simulate_observationyears <- function(n_replicates){
  n_obsyears <- n_replicates*length(corridor)*n_sections*length(habitat)*length(phase)
  df_observationyears <- data.frame(
    cameraID = as.factor(
      rep(paste0("C", 1:(n_obsyears/length(phase))),
          each = length(phase))
      ),
    corridor = relevel(as.factor(
      rep(corridor, 
          each = n_obsyears/length(corridor))),
      ref = corridor[2]
    ),
    sectionID = as.factor(
      rep(paste0("S", 1:(n_sections*length(corridor))), 
          each = n_obsyears/(length(corridor)*n_sections))
    ),
    habitat = relevel(as.factor(
      rep(habitat, 
          each = n_replicates*length(phase),
          times = n_obsyears/(n_replicates*length(phase)*length(habitat)))),
      ref = habitat[2]
    ),
    phase = relevel(as.factor(
      rep(phase,
          times = n_obsyears/length(phase))),
      ref = phase[2]
    )
  )
  df_n <- simulate_semidaily_observations(df_observationyears, n_obsyears)
  df_observationyears <- cbind(
    df_observationyears, df_n
  )
  df_observationyears
}

# PRE: n_replicates (per treatment) >= 1, hypothesis = c('A', 'B', 'C', 'BC')
# POST: returns data frame suitable to test hypothesis A, B, or C 
simulate_data <- function(hypothesis, n_replicates){
  df_full <- simulate_observationyears(n_replicates)
  # aggregate depending on hypothesis:
  if(hypothesis == "A"){
    df_wolf <- df_full %>%
      mutate(
        wolf = n01 + n11,
        no_wolf = n00 + n10
      ) %>%
      select(cameraID, corridor, sectionID, habitat, phase, wolf, no_wolf)
    return(df_wolf)
  }
  if(hypothesis %in% c("B", "C", "BC")){
    df_absent <- df_full %>%
      transmute(
        cameraID, corridor, sectionID, habitat, phase,
        caribou = "absent",
        wolf = n01,
        no_wolf = n00
      )
    
    df_present <- df_full %>%
      transmute(
        cameraID, corridor, sectionID, habitat, phase,
        caribou = "present",
        wolf = n11,
        no_wolf = n10
      )
    
    df_wolfcaribou <- bind_rows(df_absent, df_present) %>%
      mutate(
        caribou = factor(caribou, levels = c("absent", "present"))
      )
    return(df_wolfcaribou)
  } 
}

# PRE: model "A" or "B" or "C" or "BC", n := n replicates
# POST: returns glmer model for appropriate hypothesis with data simulated with 
# n_cameras at each treatment
glmm_model <- function(model, n){
  if(model == "A"){
    df_wolf <- simulate_data("A", n)
    model_A <- glmer(
      cbind(wolf, no_wolf) ~
        corridor * phase + habitat +
        (1|sectionID) + (1|cameraID),
      family = binomial,
      data = df_wolf
    )
    return(model_A)
  }
  if(model %in% c("B", "C", "BC")){
    df_wolfcaribou <- simulate_data("BC", n)
    model_B <- glmer(
      cbind(wolf, no_wolf) ~
        corridor * phase * caribou + habitat +
        (1|sectionID) + (1|cameraID),
      family = binomial,
      data = df_wolfcaribou
    )
    return(model_B)
  } 
}

# PRE: target effect 'target' %in% c("A", "B", "C"), 
# n := replicates, integer, n >= 1
# significance level alpha in [0,1], 
# n_sim := numer of simulations, integer, n_sim >= 1
# POST: power = (1 - errorTypeII), as probability of rejecting the Nullhypothesis of 
# target effect at significane level alpha given the alternative hypothesis
# is true, which is the assumption in the simulated data
estimate_power <- function(target, n, alpha, n_sim, print_output = FALSE){
  alternative_hypothesis_detected <- rep(0, n_sim)
  for(i in 1:n_sim){
    model <- glmm_model(target, n)
    if(target == "A"){pval = summary(model)$coefficients[5,4]}
    if(target == "B"){pval = summary(model)$coefficients[6,4]}
    if(target == "C"){pval = summary(model)$coefficients[9,4]}
    alternative_hypothesis_detected[i] <- as.numeric(pval <= alpha)
    
    if(print_output){
      cat("Simulation ", i, "\n")
      print(summary(model))
      cat("\n\n\n\n\n")
    }
  }
  mean(alternative_hypothesis_detected)
}



set.seed(1234)

#########################
#### Power Analysis #####
#########################

# creating a table with observed powers for hypothesis A, B, and C
range_of_replicates <- 1:10
hypotheses <- c("A", "B", "C")
power_table <- data.frame(
  n = range_of_replicates,
  setNames(
    replicate(length(hypotheses), numeric(length(range_of_replicates)), simplify = FALSE),
    hypotheses
  )
)
for (i in range_of_replicates) {
  cat("calculating n = ", i, ": ")
  for (j in hypotheses) {
    cat(j, "...")
    power_table[i, j] <- estimate_power(target = j, n = i, alpha = 0.05, n_sim = 100)
    cat(" done. ")
  }
  cat("\n")
}

# plot the table:
matplot(
  power_table$n,
  power_table[, c("A", "B", "C")],
  type = "p",
  pch = 16,
  col = c("burlywood", "coral", "cadetblue"),
  xlab = "Number of replicates (n)",
  ylab = "Power",
  ylim = c(0.5, 1)
)
abline(h=0.95)
abline(h=0.8, lty = 2)

legend(
  "bottomright",
  legend = c("Hypothesis A", "Hypothesis B", "Hypothesis C"),
  col = c("burlywood", "coral", "cadetblue"),
  pch = 16
)

##################################
#### Estimating Actual Effect ####
##################################

n_replicates <- 3
model_names <- c("A", "B")
n_sim <- 100

# to visualize effect, I run the models for n=3 for n_sim repetitions
modelA_ref <- glmm_model("A", n_replicates)
modelB_ref <- glmm_model("B", n_replicates)
coef_A_ref <- summary(modelA_ref)$coefficients[, 1]
coef_B_ref <- summary(modelB_ref)$coefficients[, 1]

coef_A <- matrix(
  0,
  nrow = n_sim,
  ncol = length(coef_A_ref),
  dimnames = list(NULL, names(coef_A_ref))
)

coef_B <- matrix(
  0,
  nrow = n_sim,
  ncol = length(coef_B_ref),
  dimnames = list(NULL, names(coef_B_ref))
)

# calculate effects
for (i in 1:n_sim){
  modelA_i <- glmm_model("A", n_replicates)
  modelB_i <- glmm_model("B", n_replicates)
  
  coef_A[i, ] <- summary(modelA_i)$coefficients[, 1]
  coef_B[i, ] <- summary(modelB_i)$coefficients[, 1]
  
}

# get all coefficients in one dataaframe
coef_data <- cbind(coef_A, coef_B)
positions <- c(1:5, 7:15)
box_cols <- rep("grey", ncol(coef_data))
box_cols[c(5, 11, 14)] <- c("burlywood", "coral", "cadetblue")

par(mar = c(12, 4, 2, 2))
boxplot(
  coef_data,
  at = positions,
  xaxt = "n",
  col = box_cols,
  ylab = "Coefficient estimate"
)

axis(1, at = positions, labels = FALSE)

# Add diagonal labels
usr <- par("usr")
par(xpd = NA)

text(
  x = positions,
  y = usr[3] - 0.04 * diff(usr[3:4]),
  labels = colnames(coef_data),
  srt = 45,
  adj = 1,
  cex = 0.8
)

par(xpd = FALSE)

abline(h=0)

legend(
  "bottomright",
  legend = c("Hypothesis A", "Hypothesis B", "Hypothesis C"),
  col = c("burlywood", "coral", "cadetblue"),
  pch = 15
)

