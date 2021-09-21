###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#---------------------Step 3. Model convergence diagnostics--------------------X
###############################################################################X

# Load packages----
library(nimble)
library(coda)

# Load samples----
# Read in the 3 chains as a list
samples <- list(
  ch1 = readRDS("models/photo_model_2021-05-13_samples_ch1.rds"),
  ch2 = readRDS("models/photo_model_2021-05-13_samples_ch2.rds"),
  ch3 = readRDS("models/photo_model_2021-05-13_samples_ch3.rds")
  )

# Convert to 'coda' object----
samps_coda <- as.mcmc.list(lapply(samples, as.mcmc))
# Gelman-Rubin statistic
(rhat <- gelman.diag(samps_coda))

# Effective sample size
effectiveSize(samps_coda)

# Create figure directory
dir.create("fig", showWarnings = FALSE)

# Traceplots
{
  pdf("fig/traceplot.pdf", width = 8, height = 6,
      onefile = TRUE, title = "Deer Model Traceplots")
  traceplot(samps_coda)
  dev.off()
}

# Density plots
{
  pdf("fig/densityplot.pdf", width = 8, height = 6,
      onefile = TRUE, title = "Deer Model Posterior Densities")
  densplot(samps_coda)
  dev.off()
}

