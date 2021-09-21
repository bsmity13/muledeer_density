###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#-----------------------------Step 1. Prepare data-----------------------------X
###############################################################################X

# Load packages----
library(tidyverse)
library(lubridate)

# Load data----
dat <- read.csv("data/muledeer_density.csv")

# Scale and center----
# Covariates to standardize
covs <- c("elev", "elev_doy", "elev_doy2", "doy",
          "tree", "grass_forb", "ndvi", "irg",
          "ndvi_grass_forb", "irg_grass_forb")
# Store means and SDs
cov_means <- apply(dat[, covs], 2, mean, na.rm = TRUE)
cov_sds <- apply(dat[, covs], 2, sd, na.rm = TRUE)

# Scale data
dat_sc <- dat
for (i in 1:length(covs)){
  dat_sc[, covs[i]] <- (dat[, covs[i]] - cov_means[covs[i]])/cov_sds[covs[i]]
}

# Prep lists for NIMBLE----
# ... Constants and data----
nimConst <- list()
# Observed deer
y <- tapply(dat$photos, list(dat$site, dat$t), identity)
## Set NAs to 0
y[is.na(y)] <- 0
## Store in list
nimConst$y <- y
# Observed occupancy
nimConst$z <- ifelse(nimConst$y > 0, 1, NA)
# Day-of-year
nimConst$doy <- unname(tapply(dat_sc$doy, list(dat_sc$t), unique))
# Habitat covariates
nimConst$elev <- unname(tapply(dat_sc$elev, list(dat_sc$site), unique))
nimConst$elev_doy <- unname(tapply(dat_sc$elev_doy, 
                                      list(dat_sc$site, dat_sc$t), unique))
nimConst$elev_doy2 <- unname(tapply(dat_sc$elev_doy2, 
                                      list(dat_sc$site, dat_sc$t), unique))
nimConst$tree <- unname(tapply(dat_sc$tree, list(dat_sc$site), unique))
nimConst$grass_forb <- unname(tapply(dat_sc$grass_forb, 
                                     list(dat_sc$site), unique))
nimConst$ndvi <- tapply(dat_sc$ndvi, list(dat_sc$site, dat_sc$t), identity)
nimConst$irg <- tapply(dat_sc$irg, list(dat_sc$site, dat_sc$t), identity)
nimConst$ndvi_grass_forb <- unname(tapply(dat_sc$ndvi_grass_forb, 
                                          list(dat_sc$site, dat_sc$t), unique))
nimConst$irg_grass_forb <- unname(tapply(dat_sc$irg_grass_forb, 
                                          list(dat_sc$site, dat_sc$t), unique))
nimConst$nosnow <- tapply(as.numeric(!dat_sc$snow), 
                          list(dat_sc$site, dat_sc$t), identity)
nimConst$snow <- tapply(as.numeric(dat_sc$snow), 
                          list(dat_sc$site, dat_sc$t), identity)

# Livestock
nimConst$lvstk <- tapply(dat_sc$live_pres, 
                         list(dat_sc$site, dat_sc$t), identity)
# Dimensions
nimConst$S <- nrow(nimConst$ndvi) #number of sites
nimConst$max_T <- ncol(nimConst$y) #maximum number of timesteps
# Camera okay?
nimConst$ok <- tapply(dat_sc$ok, list(dat_sc$site, dat_sc$t), identity)
# Camera brand (0 = cuddeback, 1 = reconyx)
nimConst$reconyx <- unname(tapply(dat_sc$camera, list(dat_sc$site), 
                           function(x){unique(as.integer(x == "reconyx"))}))
# Store site number for post hoc comparisons
nimConst$site <- unname(tapply(dat_sc$site, list(dat_sc$site), unique))

# ... Initial values----
nimInit <- list()
# z (occupancy status)
nimInit$z <- (nimConst$y > 0)
# Habitat covariates
nimInit$beta <- c(1.5, 0, 0, 0, 0, 0, 0, 0, 0)
# Occupancy covariates
nimInit$gamma <- c(-2, 0, 0, 0, 0, 0, 0)
# Random effect SD
nimInit$sigma_eta_psi <- 1
nimInit$sigma_eta_w <- 1
# Random effect values
nimInit$eta_psi <- rep(0, nimConst$S)
nimInit$eta_w <- rep(0, nimConst$S)

#Save----
# Lists for NIMBLE
saveRDS(nimConst, file = "data/nimConst.rds")
saveRDS(nimInit, file = "data/nimInit.rds")
# Scaled data
saveRDS(dat_sc, "data/data_scaled.rds")
# Scaling and centering parameters
saveRDS(covs, "data/covariate_names.rds")
saveRDS(cov_means, "data/covariate_means.rds")
saveRDS(cov_sds, "data/covariate_sds.rds")
