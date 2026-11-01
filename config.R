#### study design ####
#' note: second entry is reference level
corridor <- c("impact", "control") # 'control' will be reference level, must be 2 levels
phase <- c("after", "before") # 'before' will be reference level, must be 2 levels
n_sections <- 3 # integer
habitat <- c("mixed", "conifer") # 'conifer' will be reference level, must be 2 levels
n_semidaily_observations <- 730 # 365 days x 2 per day

#### expected effects on log-odds of wolf detection ####
#' note: effects are quantified logistically with logit(p) = log(OR) = linear effect 
#' where OR is the odds-ratio. For small changes in the odds-ratio
#' d (log(OR))/ d OR ~= 1
#' example: assume the OR is 10% higher because of an effect, OR=1.1/1 ~= exp(0.1)
beta_0 <- -4 # base probability ~= 0.047
beta_corr <- 0 # assume no pre-construction difference in the two corridors
beta_phas <- 0 # assume no change from pre to post construction in the control site
beta_cari <- 0.2 # increases p of detection wolf if caribou is present
beta_habi <- 0.3 # assume  higher change of detection in mixed rather than conifer forest
beta_corr_phas <- 0.9 # assume higher probability of wolf detection on road
beta_corr_cari <- 0 # assume no wolf in impact corridor pre construction respond same to caribou as in control corridor
beta_cari_phas <- 0 # assume no change in wolf detection given caribou in ref level
beta_cari_corr_phas <- 1.3 # assume wolfs can respond better to caribou in impact corridor post construction
beta_sigma <- 0.67 # residual
beta_sigma_sectionID <- 1 # mixed effect per sectionID
beta_sigma_cameraID <- 1 # mixed effect per cameraID


#### expected effects on log-odds of caribou detection ####
gamma_0 <- -3 # basic rate, roughly 0.119
gamma_corr <- 0 # assume caribou detection independent of corridor
gamma_phas <- 0 # assume caribou detection independent of pre or post construction in control
gamma_habi <- 0 # assume no habitat preference of caribou (for simplicity)
gamma_corr_phas <- -0.4 #decreased likelihood of detecting caribou post-construction
gamma_sigma <- 0.71 # residual
gamma_sigma_sectionID <- 1 # mixed effect per sectionID
gamma_sigma_cameraID <- 1 # mixed effect per cameraID