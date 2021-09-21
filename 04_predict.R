###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#---------------------------Step 4. Model predictions--------------------------X
###############################################################################X

# Load packages----
library(tidyverse)
library(lubridate)

# Load data----
nimConst <- readRDS("data/nimConst.rds")
nimInit <- readRDS("data/nimInit.rds")
dat <- read.csv("data/muledeer_density.csv")
dat_sc <- readRDS("data/data_scaled.rds")

# Load samples----
# Read in the 3 chains as a list
samples <- list(
  ch1 = readRDS("models/photo_model_2021-05-13_samples_ch1.rds"),
  ch2 = readRDS("models/photo_model_2021-05-13_samples_ch2.rds"),
  ch3 = readRDS("models/photo_model_2021-05-13_samples_ch3.rds")
)

samples2 <- list(
  ch1 = readRDS("models/photo_model_2021-05-13_samples2_ch1.rds"),
  ch2 = readRDS("models/photo_model_2021-05-13_samples2_ch2.rds"),
  ch3 = readRDS("models/photo_model_2021-05-13_samples2_ch3.rds")
)

# Thin posterior----
# We want to use the posterior to make model predictions so that we capture the
# uncertainty in the prediction, but the entire posterior is too large, so we
# will thin it here to just 100 samples per chain.

# The chains are elements of the list "samples", so we can use lapply() to thin
# each element of the list.
post_list <- lapply(samples, function(x){
  # Thin to 300 samples
  keep <- seq(from = 1, to = nrow(x), length.out = 100)
  # Subset rows
  res <- x[keep, ]
  # Return result
  return(res)
})

post_list2 <- lapply(samples2, function(x){
  # Thin to 1000 samples
  keep <- seq(from = 1, to = nrow(x), length.out = 100)
  # eta_w columns
  w <- grep("eta_w[", colnames(x), fixed = TRUE)
  # eta_psi columns
  psi <- grep("eta_psi[", colnames(x), fixed = TRUE)
  # Subset rows
  res <- x[keep, c(w, psi)]
  # Return result
  return(res)
})

# Now rbind the lists together and convert to a data.frame
post <- as.data.frame(do.call(rbind, post_list))
post2 <- as.data.frame(do.call(rbind, post_list2))

# Predict state variables from model ----

# Give dat & dat_sc a unique row ID so we can summarize later
dat$row_id <- 1:nrow(dat)
dat_sc$row_id <- 1:nrow(dat_sc)

# Convert site to index "i"
sites <- data.frame(site = sort(unique(dat_sc$site)),
                    i = 1:length(unique(dat_sc$site)))
dat_sc <- left_join(dat_sc, sites)

# ... function to predict for 1 row ----

# x is the row of the data.frame with all the covariates we need to predict

# post is the data.frame with posterior samples of fixed effects
# post2 is the data.frame with posterior samples of random effects
#   Note: currently ignoring site-level RE for general prediction

predict_row <- function(x, post, post2) {
  # Predict Pr(occupancy)
  pred_occ <- plogis(post$`gamma[1]` + 
                       post$`gamma[2]` * x$live_pres +
                       post$`gamma[3]` * x$snow +
                       post$`gamma[4]` * x$elev_doy +
                       post$`gamma[5]` * x$elev_doy2 +
                       post$`gamma[6]` * x$doy +
                       post$`gamma[7]` * (x$doy)^2 +
                       0 # post2[[paste0("eta_psi[", x$i, "]")]]
                     )
  
  # Predict E[abundance | occupied]
  pred_abund_occ <-  exp(post$`beta[1]` + 
                           post$`beta[2]` * x$doy +
                           post$`beta[3]` * (x$doy)^2 +
                           post$`beta[4]` * x$grass_forb +
                           post$`beta[5]` * x$ndvi +
                           post$`beta[6]` * x$irg + 
                           post$`beta[7]` * x$ndvi_grass_forb +
                           post$`beta[8]` * x$irg_grass_forb +
                           post$`beta[9]` * x$tree +
                           0 # post2[[paste0("eta_w[", x$i, "]")]]
                         )
  
  # Predict E[abundance] = Pr(occupancy) * E[abundance | occupancy]
  pred_abund <- pred_occ * pred_abund_occ
  
  # Build data.frame to return predicted values with row_id
  # Start with just row_id
  res <- data.frame(row_id = rep(x$row_id, nrow(post)))
  # Now add predicted state variables to the end
  res$occ <- pred_occ
  res$abund_occ <- pred_abund_occ
  res$abund <- pred_abund
  
  # Return
  return(res)
}

# Example code:
example_x <- dat_sc[1,]

predict_row(x = example_x, post = post, post2 = post2)

# ... posterior prediction for all data----
# Split the data into a list by row id (one list element per row)
dat_sc_list <- split(dat_sc, dat_sc$row_id)

# Now use lapply to predict across all
pred <- lapply(dat_sc_list, predict_row, post = post, post2 = post2)

# Combine list into a single data.frame
pred <- bind_rows(pred)

# ... join original covariates----
# Now let's add the original covariates to the predictions so that we can decide
# how we want to summarize/plot

pred_dat <- left_join(pred, dat, by = "row_id")

# Save----
# Save pred_dat so we can access it later for plotting
saveRDS(pred_dat, "data/predicted_data.rds")
