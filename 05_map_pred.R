###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#----------------------------Step 5. Map predictions---------------------------X
###############################################################################X

# Load packages----
library(raster)
library(sf)
library(tidyverse)
library(lubridate)

# Source helper functions ----
source("99_fun.R")

# Load rasters ----
# Prediction rasters
ndvi_stack <- stack("geo/NDVI_stack.tif")
irg_stack <- stack("geo/IRG_stack.tif")
elev <- raster("geo/elev.tif")
grass_forb <- raster("geo/grass_forb.tif")
tree <- raster("geo/tree.tif")

# Load MCMC samples----
# Read in the 3 chains as a list
samples <- list(
  ch1 = readRDS("models/photo_model_2021-05-13_samples_ch1.rds"),
  ch2 = readRDS("models/photo_model_2021-05-13_samples_ch2.rds"),
  ch3 = readRDS("models/photo_model_2021-05-13_samples_ch3.rds")
)

# Load scaling parameters ----
covs <- readRDS("data/covariate_names.rds")
cov_means <- readRDS("data/covariate_means.rds")
cov_sds <- readRDS("data/covariate_sds.rds")

# Prediction dates ----
pd <- as.Date(c("2019-04-01", "2019-05-01", "2019-06-01"))
# Day-of-year
jd <- yday(pd)

# Photos-to-density constant----
#(x deer pic/1day) * (2 seconds/1 pic) * (1 day/86400 seconds) * 
#   (1/220 ft2) * (107639 ft2/1 ha)
photos_to_dens <- (2/1) * (1/86400) * (1/220) * (107639/1)

# Get parameter estimates ----
# Thin to 100 samples per chain
gamma_samples <- match_samples("gamma", samples, thin_total = 100)
beta_samples <- match_samples("beta", samples, thin_total = 100)

# Convert rasters to data.frame
elev_df <- elev %>% 
  as.data.frame(xy = TRUE) %>% 
  # Need elev in ft
  mutate(elev = elev * 3.28084)

gf_df <- grass_forb %>% 
  as.data.frame(xy = TRUE)

tree_df <- as.data.frame(tree, xy = TRUE) %>% 
  mutate(tree = tree/100)

# Base data
base_dat <- elev_df %>% 
  left_join(gf_df, by = c("x", "y")) %>% 
  left_join(tree_df, by = c("x", "y")) %>% 
  mutate(lvstk = 0,
         snow = 0)

# Repeat base dat for each date
date_dat <- data.frame()
for (i in 1:length(pd)) {
  date_dat <- rbind(
    date_dat,
    mutate(base_dat, date = pd[i])
  )
}
date_dat$doy <- yday(date_dat$date)

ndvi_df <- as.data.frame(ndvi_stack, xy = TRUE) %>% 
  pivot_longer(cols = c(-x, -y),
               names_to = "name",
               values_to = "ndvi") %>% 
  mutate(doy = case_when(
    name == "NDVI_stack.1" ~ jd[1],
    name == "NDVI_stack.2" ~ jd[2],
    name == "NDVI_stack.3" ~ jd[3]
  )) %>% 
  select(-name)

# IRG data
irg_df <- as.data.frame(irg_stack, xy = TRUE) %>% 
  pivot_longer(cols = c(-x, -y),
               names_to = "name",
               values_to = "irg") %>% 
  mutate(doy = case_when(
    name == "IRG_stack.1" ~ jd[1],
    name == "IRG_stack.2" ~ jd[2],
    name == "IRG_stack.3" ~ jd[3]
  )) %>% 
  select(-name)

# All data
all_dat <- date_dat %>% 
  left_join(ndvi_df, by = c("x", "y", "doy")) %>% 
  left_join(irg_df, by = c("x", "y", "doy"))

# Composite variables
all_dat <- all_dat %>% 
  mutate(elev_doy = elev * doy,
         elev_doy2 = elev * (doy)^2,
         ndvi_grass_forb = ndvi * grass_forb,
         irg_grass_forb = irg * grass_forb)

# Scale and center
dat_sc <- all_dat
for (i in 1:length(covs)){
  dat_sc[, covs[i]] <- (all_dat[, covs[i]] - cov_means[covs[i]])/cov_sds[covs[i]]
}

# One cell has NA for NDVI
dat_sc <- dat_sc %>% 
  filter(!is.na(ndvi))

# ... predict ----
dat_sc$row <- 1:nrow(dat_sc)
## **This takes a few minutes**
pred <- lapply(1:nrow(dat_sc), function(i) {
  cat("\r", i, "of", nrow(dat_sc), "     ")
  # Get row of data
  x <- dat_sc[i, ]
  
  # Predict across data
  # psi
  psi <- plogis(gamma_samples[, 1] + 
                  gamma_samples[, 2] * x$lvstk +
                  gamma_samples[, 3] * x$snow +
                  gamma_samples[, 4] * x$elev_doy +
                  gamma_samples[, 5] * x$elev_doy2 +
                  gamma_samples[, 6] * x$doy +
                  gamma_samples[, 7] * (x$doy)^2)
  
  # w
  w <- exp(beta_samples[, 1] + 
             beta_samples[, 2] * x$doy +
             beta_samples[, 3] * (x$doy)^2 +
             beta_samples[, 4] * x$grass_forb +
             beta_samples[, 5] * x$ndvi +
             beta_samples[, 6] * x$irg + 
             beta_samples[, 7] * x$ndvi_grass_forb +
             beta_samples[, 8] * x$irg_grass_forb +
             beta_samples[, 9] * x$tree)
  
  # lambda
  lambda <- psi * w
  
  # Data.frame to return
  res <- data.frame(
    row = x$row,
    x = x$x,
    y = x$y,
    date = x$date,
    psi_lwr = quantile(psi, 0.05),
    psi_mean = mean(psi),
    psi_upr = quantile(psi, 0.90),
    w_lwr = quantile(w, 0.05),
    w_mean = mean(w),
    w_upr = quantile(w, 0.90),
    lambda_lwr = quantile(lambda, 0.05),
    lambda_mean = mean(lambda),
    lambda_upr = quantile(lambda, 0.90)
  )
  row.names(res) <- NULL
  # Return
  return(res)
}) %>% 
  bind_rows()

# Save ----
saveRDS(pred, "data/map_prediction.rds")
saveRDS(dat_sc, "data/map_data.rds")