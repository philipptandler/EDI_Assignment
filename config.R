## study design ##
corridor <- c("impact", "control") # 'control' will be reference level
phase <- c("after", "before") # 'before' will be reference level
n_sections <- 3
habitat <- c("mixed", "conifer") # 'conifer' will be reference level 

## expected effects on wolf detection ##
#' note: effects are quantified logistically with logit(p) = log(OR) = linear effect 
#' where OR is the odds-ratio. For small changes in the odds-ratio
#' d (log(OR))/ d OR ~= 1
#' example: assume the OR is 10% higher because of an effect, OR=1.1/1 ~= exp(0.1)
beta_0 <- -3 # base probability ~= 0.047
beta_corr <- 0 # assume no pre-construction difference in the two corridors
beta_phas <- 0 # assume no change from pre to post construction in the control site
beta_cari <- 0.2 # increases p of detection wolf if caribou is present
beta_habi <- 0.3 # assume  higher change of detection in mixed rather than conifer forest
beta_corr_phas <- 0.9 # assume higher probability of wolf detection on road
beta_corr_cari <- 0 # assume no wolf in impact corridor pre construction respond same to caribou as in control corridor
beta_cari_phas <- 0 # assume no change in wolf detection given caribou in ref level
beta_cari_corr_phas <- 0.1 # assume wolfs can respond better to caribou in impact corridor post constuction


## expected effects on caribou detection ##
alpha_0 <- -2 # basic rate, roughly 0.119
alpha_corr <- 0 # assume caribou detection independent of corridor
alpha_phas <- 0 # assume caribou detection independent of pre or post construction in control
alpha_habi <- 0 # assume no habitat preference of caribou (for simplicity)
alpha_corr_phas <- -0.4 # assume 
