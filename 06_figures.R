###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#--------------------------------Step 6. Figures-------------------------------X
###############################################################################X

#Load packages----
library(raster)
library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(ggspatial)
library(ragg) # for quality plotting device on Windows

# Source helper functions ----
source("99_fun.R")

#Load data----
nimConst <- readRDS("data/nimConst.rds")
nimInit <- readRDS("data/nimInit.rds")
dat <- read.csv("data/muledeer_density.csv")
dat_sc <- readRDS("data/data_scaled.rds")

# ... photos-to-density constant----
#(x deer pic/1day) * (2 seconds/1 pic) * (1 day/86400 seconds) * 
#   (1/220 ft2) * (107639 ft2/1 ha)
photos_to_dens <- (2/1) * (1/86400) * (1/220) * (107639/1)

# ... samples----
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

# ... covariate mean and SD----
covs <- readRDS("data/covariate_names.rds")
cov_means <- readRDS("data/covariate_means.rds")
cov_sds <- readRDS("data/covariate_sds.rds")

# ... prepped predictions----

# See script 04_predict.R for details
pred_dat <- readRDS("data/predicted_data.rds")
# See script 05_map_pred.R for details
map_pred <- readRDS("data/map_prediction_final.rds")

# ... rasters ----
ndvi_stack <- stack("geo/NDVI_stack.tif")
irg_stack <- stack("geo/IRG_stack.tif")
elev <- raster("geo/elev.tif")
grass_forb <- raster("geo/grass_forb.tif")
tree <- raster("geo/tree.tif")

# ... shapefiles ----
sites <- st_read("geo/sites.shp")

# Figure 1 ----
# Study area map
# A - Elevation
# B - % grass/forb

# ** Important note: for final publication figure, these rasters were
#     loaded with 30m resolution. These data are not uploaded to GitHub
#     to save space. **

elev_df <- elev %>% 
  as.data.frame(xy = TRUE)
gf_df <- grass_forb %>% 
  as.data.frame(xy = TRUE)

fig1A <- elev_df %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = elev)) +
  geom_sf(data = sites, size = 0.75) +
  annotation_north_arrow(aes(location = "tl"),
                         height = unit(1.5 * 0.75, "cm"),
                         width = unit(1 * 0.75, "cm")) +
  coord_sf(crs = st_crs(32612),
           expand = FALSE) +
  scale_fill_gradientn(name = "Elevation (m)", 
                       colors = terrain.colors(7),
                       na.value = NA) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

fig1B <- gf_df %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = grass_forb)) +
  geom_sf(data = sites, size = 0.75) +
  coord_sf(crs = st_crs(32612),
           expand = FALSE) +
  scale_fill_gradient(name = "Grass + Forb (%)", 
                      low = "white",
                      high = "forestgreen",
                      na.value = NA) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

fig1 <- fig1A + fig1B +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("fig/Figure1.tiff", plot = fig1, device = agg_tiff,
       width = 6, height = 4, units = "in",
       dpi = 600, compression = "lzw", scale = 1.4)

# Figure 2 ----
# Illustration of camera setup and example photos
# Not made in R

# Figure 3 ----
# Occupancy Figure

# ... Fig. 3A ----
# Nice labels for parameters
gamma_parm_nms <- c("Intercept",  
                    "Livestock",
                    "Snow",
                    "Elev:DoY",
                    expression("Elev:(DoY)"^2),
                    "DoY",
                    expression("(DoY)"^2))

# Get parameters from samples
gamma_samples <- match_samples("gamma", samples)
# Means
g <- apply(gamma_samples, 2, mean)

occ_params <- gamma_samples %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, 0.05),
            upr = quantile(value, 0.95)) %>% 
  mutate(order = c(1:ncol(gamma_samples)),
         param = factor(order))
levels(occ_params$param) <- gamma_parm_nms

fig3A <- ggplot(occ_params, aes(x = mean, y = param)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50",
             size = 1) +
  geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0.1) +
  geom_point() +
  scale_y_discrete(breaks = levels(occ_params$param),
                   labels = gamma_parm_nms) +
  ylab(NULL) +
  xlab(expression("Parameter Estimate" ~ (gamma))) +
  theme_bw()

# ... Fig. 3B ----
# y-limits for Fig. 3B&C
ylim_3bc <- c(0, 0.25)

# Bin doy day for summary/plot
doy_bins <- pretty(pred_dat$doy, 20)

occ_dat <- pred_dat %>% 
  # Bin doy day
  mutate(doy_bin = cut_width(doy, 
                             width = 5,
                             center = doy_bins[1],
                             dig.lab = 4),
         doy_lab = doy_bin)

#Re-do levels of doy_lab for nice plotting
levels(occ_dat$doy_lab) <- as.character(
  doy_bins[2:(length(doy_bins - 1))]
)

# Convert doy day labels to dates
occ_dat$doy_date <- ymd("2019-01-01") #dummy date
# Replace dummy date with doy day
yday(occ_dat$doy_date) <- as.numeric(as.character(occ_dat$doy_lab))

# Bin elevation
elev_bins <- pretty(occ_dat$elev, 20)
# Convert data to bins
elev_dat <- occ_dat %>% 
  # Bin elevation
  mutate(elev_bin = cut_width(elev, 
                              width = 200,
                              boundary = elev_bins[1],
                              dig.lab = 4),
         elev_lab = elev_bin)
# Give bins nice labels
levels(elev_dat$elev_lab) <- na.omit(
  as.character((elev_bins + lead(elev_bins))/2)
)

# Summarize
elev_summ <- elev_dat %>% 
  group_by(doy_date, elev_lab) %>% 
  # Summarize mean and 90% CI
  summarize(q05 = quantile(occ, 0.05),
            mean = mean(occ),
            q95 = quantile(occ, 0.95))

fig3B <- elev_summ %>% 
  filter(elev_lab %in% c("5500", "7500", "8300")) %>% 
  ggplot(aes(x = doy_date, y = mean, color = elev_lab,
             ymin = q05, ymax = q95, fill = elev_lab)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  xlab(NULL) +
  ylab("Predicted Pr(Occupancy)") +
  scale_color_brewer(name = "Elevation\nBin (ft)",
                     type = "qual",
                     palette = 1) +
  scale_fill_brewer(name = "Elevation\nBin (ft)",
                    type = "qual",
                    palette = 1)  +
  coord_cartesian(ylim = ylim_3bc) +
  theme_bw()


# ... Fig. 3C ----
snow_live_dat <- occ_dat %>% 
  # Create composite variable for snow or livestock
  mutate(snow_live = case_when(
    snow ~ "Snow",
    live_pres == 1 ~ "Livestock",
    TRUE ~ "Neither"
  ))

# Summarize
snow_live_summ <- snow_live_dat %>% 
  group_by(doy_date, snow_live) %>% 
  # Summarize several quantiles
  summarize(lwr = quantile(occ, 0.05),
            mean = mean(occ),
            upr = quantile(occ, 0.95))

fig3C <- ggplot(snow_live_summ, aes(x = doy_date, y = mean, 
                                    ymin = lwr, ymax = upr, 
                                    color = factor(snow_live),
                                    fill = factor(snow_live))) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line(size = 1) +
  xlab("Date") +
  ylab("Predicted Pr(Occupancy)") +
  scale_color_brewer(name = "Disturbance",
                     type = "qual",
                     palette = 6) +
  scale_fill_brewer(name = "Disturbance",
                    type = "qual",
                    palette = 6) +
  coord_cartesian(ylim = ylim_3bc) +
  theme_bw()

# ... Combine ----

fig3 <- fig3A + (fig3B / fig3C) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave("fig/Figure3.tiff", plot = fig3, device = agg_tiff,
       width = 6, height = 5, units = "in",
       dpi = 600, compression = "lzw", scale = 1.1)


# Figure 4 ----

# ... data prep ----
ao_dat <- occ_dat %>% 
  mutate(dens_occ = abund_occ * photos_to_dens)

# Need ranges and quantiles for grass_forb, ndvi, and irg
gf_range <- c(0, 0.6)
gf_hi <- quantile(pred_dat$grass_forb, 0.9)
gf_lo <- quantile(pred_dat$grass_forb, 0.1)

ndvi_range <- c(0, 0.9)
ndvi_hi <- quantile(pred_dat$ndvi, 0.9)
ndvi_lo <- quantile(pred_dat$ndvi, 0.1)

irg_range <- c(0, 1)
irg_hi <- quantile(pred_dat$irg, 0.9)
irg_lo <- quantile(pred_dat$irg, 0.1)

# Bin grass_forb
gf_bins <- seq(gf_range[1], gf_range[2], by = 0.1)
# Bin NDVI
ndvi_bins <- seq(ndvi_range[1], ndvi_range[2], by = 0.1)
# Bin IRG
irg_bins <- seq(irg_range[1], irg_range[2], by = 0.1)

# Which bins contain each quantile
gf_bin_lo <- max(which(gf_bins < gf_lo))
gf_bin_hi <- min(which(gf_bins > gf_hi)) - 1

ndvi_bin_lo <- max(which(ndvi_bins < ndvi_lo))
ndvi_bin_hi <- min(which(ndvi_bins > ndvi_hi)) - 1

irg_bin_lo <- max(which(irg_bins < irg_lo))
irg_bin_hi <- min(which(irg_bins > irg_hi)) - 1

# Convert data to bins
ao_dat <- ao_dat %>% 
  # Bin elevation
  mutate(gf_bin = cut_width(grass_forb, 
                            width = 0.1,
                            boundary = gf_bins[1],
                            dig.lab = 2),
         gf_lab = gf_bin,
         ndvi_bin = cut_width(ndvi, 
                              width = 0.1,
                              boundary = ndvi_bins[1],
                              dig.lab = 2),
         ndvi_lab = ndvi_bin,
         irg_bin = cut_width(irg, 
                             width = 0.1,
                             boundary = irg_bins[1],
                             dig.lab = 2),
         irg_lab = irg_bin)

# Give bins nice labels
levels(ao_dat$gf_lab) <- na.omit(
  as.character((gf_bins + lead(gf_bins))/2)
)

levels(ao_dat$ndvi_lab) <- na.omit(
  as.character((ndvi_bins + lead(ndvi_bins))/2)
)

levels(ao_dat$irg_lab) <- na.omit(
  as.character((irg_bins + lead(irg_bins))/2)
)

# Now create groups for Figs B and C
ao_dat <- ao_dat %>% 
  mutate(gf_lo = as.numeric(gf_bin) == gf_bin_lo,
         gf_hi = as.numeric(gf_bin) == gf_bin_hi,
         ndvi_lo = as.numeric(ndvi_bin) == ndvi_bin_lo,
         ndvi_hi = as.numeric(ndvi_bin) == ndvi_bin_hi,
         irg_lo = as.numeric(irg_bin) == irg_bin_lo,
         irg_hi = as.numeric(irg_bin) == irg_bin_hi,
  ) %>% 
  mutate(group_B = case_when(
    gf_lo & ndvi_lo ~ "Low - Low",
    gf_lo & ndvi_hi ~ "Low - High",
    gf_hi & ndvi_lo ~ "High - Low",
    gf_hi & ndvi_hi ~ "High - High",
    TRUE ~ NA_character_),
    group_C = case_when(
      gf_lo & irg_lo ~ "Low - Low",
      gf_lo & irg_hi ~ "Low - High",
      gf_hi & irg_lo ~ "High - Low",
      gf_hi & irg_hi ~ "High - High",
      TRUE ~ NA_character_))

# ... Fig. 4A ----
# Nice labels for parameters
beta_parm_nms <- c("Intercept",
                   "DoY",
                   expression("(DoY)"^2),
                   "Grass+Forb",
                   "NDVI",
                   "IRG",
                   "NDVI:Grass+Forb",
                   "IRG:Grass+Forb",
                   "Trees")

# Get parameters from samples
beta_samples <- match_samples("beta", samples)

# Means
b <- apply(beta_samples, 2, mean)

hs_params <- beta_samples %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, 0.05),
            upr = quantile(value, 0.95)) %>% 
  mutate(order = c(1:ncol(beta_samples)),
         param = factor(order))
levels(hs_params$param) <- beta_parm_nms

fig4A <- ggplot(hs_params, aes(y = param, x = mean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50",
             size = 1) +
  geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0.1) +
  geom_point() +
  ylab(NULL) +
  xlab(expression("Parameter Estimate" ~ (beta))) +
  scale_y_discrete(breaks = levels(hs_params$param),
                   labels = beta_parm_nms) +
  theme_bw()

# ... Fig. 4B ----
# y-limits for Fig. 4B&C
ylim_4bc <- c(0, 2)

summ_4B <- ao_dat %>% 
  filter(!is.na(group_B)) %>% 
  group_by(doy_date, group_B) %>% 
  # Summarize mean and 90% CI
  summarize(q05 = quantile(dens_occ, 0.05),
            mean = mean(dens_occ),
            q95 = quantile(dens_occ, 0.95))

fig4B <- summ_4B %>% 
  ggplot(aes(x = doy_date, y = mean, color = group_B,
             ymin = q05, ymax = q95, fill = group_B)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  xlab(NULL) +
  ylab("E[Density | Occupancy] (deer/ha)") +
  scale_color_brewer(name = "G+F - NDVI",
                     type = "qual",
                     palette = 1) +
  scale_fill_brewer(name = "G+F - NDVI",
                    type = "qual",
                    palette = 1)  +
  coord_cartesian(ylim = ylim_4bc) +
  theme_bw()

# ... Fig. 4C ----
summ_4C <- ao_dat %>% 
  filter(!is.na(group_C)) %>% 
  group_by(doy_date, group_C) %>% 
  # Summarize mean and 90% CI
  summarize(q05 = quantile(dens_occ, 0.05),
            mean = mean(dens_occ),
            q95 = quantile(dens_occ, 0.95))

fig4C <- summ_4C %>% 
  ggplot(aes(x = doy_date, y = mean, color = group_C,
             ymin = q05, ymax = q95, fill = group_C)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  xlab("Date") +
  ylab("E[Density | Occupancy] (deer/ha)") +
  scale_color_brewer(name = "G+F - IRG",
                     type = "qual",
                     palette = 2) +
  scale_fill_brewer(name = "G+F - IRG",
                    type = "qual",
                    palette = 2)  +
  coord_cartesian(ylim = ylim_4bc) +
  theme_bw()

# ... combine ----

fig4 <- fig4A + (fig4B / fig4C) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave("fig/Figure4.tiff", plot = fig4, device = agg_tiff,
       width = 6, height = 5, units = "in",
       dpi = 600, compression = "lzw", scale = 1.2)

# Figure 5 ----

# ... prep data ----
map_dat <- map_pred %>% 
  mutate(dens = lambda_mean * photos_to_dens,
         log_dens = log(dens)) %>% 
  split(.$date)

# Prediction dates
pd <- as.Date(c("2019-04-01", "2019-05-01", "2019-06-01"))
# Day-of-year
jd <- yday(pd)

# ... Fig. 5A ----
# Function to make ggplot
fun_5A <- function(data) {
  ggplot(data, aes(x = x, y = y, fill = ndvi)) +
    geom_raster() +
    coord_sf(crs = st_crs(32612),
             expand = FALSE) +
    scale_fill_gradient(name = "NDVI", 
                        low = "white",
                        high = "forestgreen",
                        na.value = NA,
                        limits = c(0, 1)) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# ... ... panel 1 ----
dat5A_1 <- as.data.frame(ndvi_stack[[1]], xy = TRUE)
names(dat5A_1)[3] <- "ndvi"

fig5A_1 <- fun_5A(dat5A_1) +
  annotation_north_arrow(aes(location = "tl"),
                         height = unit(1.5 * 0.75, "cm"),
                         width = unit(1 * 0.75, "cm")) +
  theme(axis.text.x = element_blank()) +
  ggtitle("April 1, 2019")

# ... ... panel 2 ----
dat5A_2 <- as.data.frame(ndvi_stack[[2]], xy = TRUE)
names(dat5A_2)[3] <- "ndvi"

fig5A_2 <- fun_5A(dat5A_2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("May 1, 2019")

# ... ... panel 3 ----
dat5A_3 <- as.data.frame(ndvi_stack[[3]], xy = TRUE)
names(dat5A_3)[3] <- "ndvi"

fig5A_3 <- fun_5A(dat5A_3) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("June 1, 2019")

# ... ... combine ----
fig5A <- fig5A_1 + fig5A_2 + fig5A_3 +
  plot_layout(nrow = 1, guides = "auto")

# ... Fig. 5B ----
# Function to make ggplot
fun_5B <- function(data) {
  ggplot(data, aes(x = x, y = y, fill = irg)) +
    geom_raster() +
    coord_sf(crs = st_crs(32612),
             expand = FALSE) +
    scale_fill_gradientn(name = "IRG", 
                         colors = rev(terrain.colors(7)),
                         na.value = NA,
                         limits = c(0, 1)) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# ... ... panel 1 ----
dat5B_1 <- as.data.frame(irg_stack[[1]], xy = TRUE)
names(dat5B_1)[3] <- "irg"

fig5B_1 <- fun_5B(dat5B_1) +
  theme(axis.text.x = element_blank())

# ... ... panel 2 ----
dat5B_2 <- as.data.frame(irg_stack[[2]], xy = TRUE)
names(dat5B_2)[3] <- "irg"

fig5B_2 <- fun_5B(dat5B_2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# ... ... panel 3 ----
dat5B_3 <- as.data.frame(irg_stack[[3]], xy = TRUE)
names(dat5B_3)[3] <- "irg"

fig5B_3 <- fun_5B(dat5B_3) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# ... ... combine ----
fig5B <- fig5B_1 + fig5B_2 + fig5B_3 +
  plot_layout(nrow = 1, guides = "auto")

# ... Fig. 5C ----
# Function to make ggplot
fun_5C <- function(data) {
  ggplot(data, aes(x = x, y = y, fill = dens)) +
    geom_raster() +
    coord_sf(crs = st_crs(32612),
             expand = FALSE) +
    scale_fill_viridis_c(name = " Density \n (deer/ha) ", 
                         na.value = NA,
                         limits = c(0, 0.085)) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
}

# ... ... panel 1 ----
fig5C_1 <- fun_5C(map_dat[[1]])

# ... ... panel 2 ----
fig5C_2 <- fun_5C(map_dat[[2]]) +
  theme(axis.text.y = element_blank())

# ... ... panel 3 ----
fig5C_3 <- fun_5C(map_dat[[3]]) +
  theme(axis.text.y = element_blank())

# ... ... combine ----
fig5C <- fig5C_1 + fig5C_2 + fig5C_3 +
  plot_layout(nrow = 1, guides = "auto")

# ... combine ----
# Adds a dark border for making sure we filled the page area.
# Not intended for use in final figure.
# theme_border <- theme_gray() + 
#   theme(plot.background = element_rect(fill = NA, colour = 'black', size = 3))

fig5 <- (fig5A / fig5B / fig5C)  +
  # plot_annotation(theme = theme_border) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

ggsave("fig/Figure5.tiff", plot = fig5, device = agg_tiff,
       width = 6, height = 8, units = "in",
       dpi = 600, compression = "lzw", scale = 1.5)

# Supplemental Figures ----


# ... nuisance parameters----
# Nice labels for parameters
e_parm_nms <- c("SD_w",
                "SD_psi")

# Get parameters from samples
eta_samples <- cbind(
  match_samples("sigma_eta_w", samples, exact = TRUE),
  match_samples("sigma_eta_psi", samples)
)
e <- apply(eta_samples, 2, mean)

eta_params <- eta_samples %>%
  as.data.frame() %>%
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, 0.05),
            upr = quantile(value, 0.95)) %>% 
  mutate(order = c(1:ncol(eta_samples)),
         param = factor(order))
levels(eta_params$param) <- e_parm_nms

eta_parm_fig <- ggplot(eta_params, aes(x = param, y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50",
             size = 1) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
  geom_point() +
  xlab(NULL) +
  ylab("Parameter Estimate") +
  scale_x_discrete(breaks = e_parm_nms, 
                   labels = c(expression(sigma[w]), 
                              expression(sigma[psi]))) +
  theme_bw()

ggsave("fig/FigureS1-1.tiff", plot = eta_parm_fig, device = agg_tiff,
       width = 6, height = 4, units = "in",
       dpi = 300, compression = "lzw", scale = 0.8)
