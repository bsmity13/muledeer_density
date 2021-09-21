###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#-------------------------------Step 2. Fit model------------------------------X
###############################################################################X

# Load packages----
library(nimble)

# Load data----
nimConst <- readRDS("data/nimConst.rds")
nimInit <- readRDS("data/nimInit.rds")

# Create model directory ----
dir.create("models", showWarnings = FALSE)

# Source definition----
source("99_NIMBLE_defn.R")

# Create model definition----
photo <- nimbleModel(code = photoCode, name = "photo", 
                     constants = nimConst, inits = nimInit, 
                     check = FALSE, calculate = FALSE)

# Check for initialized values
photo$initializeInfo()

# Configure the MCMC
photoConf <- configureMCMC(photo, print = TRUE,
                           monitors = c("beta", 
                                        "gamma",
                                        "sigma_eta_w",
                                        "sigma_eta_psi"),
                           thin = 15,
                           monitors2 = c("eta_w",
                                         "eta_psi",
                                         "psi",
                                         "lambda"),
                           thin2 = 150)

# Manually edit MCMC samplers
# ... block sample beta parms----
photoConf$removeSampler("beta")
photoConf$addSampler("beta", "AF_slice", 
                     control = list(sliceWidths = rep(0.01, 9),
                                    sliceAdaptFactorInterval = 100,
                                    sliceMaxSteps = 300,
                                    maxContractions = 300,
                                    maxContractionsWarning = FALSE))

# ... block sample gamma parms----
photoConf$removeSampler("gamma")
photoConf$addSampler("gamma", "AF_slice", 
                     control = list(sliceWidths = rep(0.01, 7),
                                    sliceAdaptFactorInterval = 100,
                                    sliceMaxSteps = 300,
                                    maxContractions = 300,
                                    maxContractionsWarning = FALSE))

# ... block sample RE SD----
photoConf$removeSampler(c("sigma_eta_psi", "sigma_eta_w"))
photoConf$addSampler(c("sigma_eta_psi", "sigma_eta_w"), 
                     "AF_slice",
                     control = list(sliceWidths = c(0.01, 0.01),
                                    sliceAdaptFactorInterval = 100,
                                    sliceMaxSteps = 300,
                                    maxContractions = 300,
                                    maxContractionsWarning = FALSE))

# ... slice sample REs----
photoConf$removeSampler("eta_psi")
photoConf$addSampler("eta_psi", type = "slice", 
                     scalarComponents = TRUE,
                     control = list(adaptive = TRUE,
                                    adaptInterval = 100,
                                    sliceWidth = 0.1,
                                    sliceMaxSteps = 300,
                                    maxContractions = 300,
                                    maxContractionsWarning = FALSE))

photoConf$removeSampler("eta_w")
photoConf$addSampler("eta_w", type = "slice", 
                     scalarComponents = TRUE,
                     control = list(adaptive = TRUE,
                                    adaptInterval = 100,
                                    sliceWidth = 0.1,
                                    sliceMaxSteps = 300,
                                    maxContractions = 300,
                                    maxContractionsWarning = FALSE))

# ... Build the MCMC ----
photoMCMC <- buildMCMC(photoConf)

# ... Compile ----
# Compile the model
Cphoto <- compileNimble(photo)

# Compile with MCMC
CphotoMCMC <- compileNimble(photoMCMC)

# MCMC parameters----
nchains <- 3
niter <- 100000
nburnin <- 25000
nthin <- 15

# Run the MCMC----
# Lists to hold timing
mcmc_start <- list()
mcmc_end <- list()
mcmc_time <- list()
# Loop over chains
for (i in 1:nchains) {
  # Report status
  cat("Running chain", i, "... \n")
  # Start timer
  mcmc_start[[i]] <- Sys.time()
  # Call run()
  CphotoMCMC$run(niter = niter, nburnin = nburnin, thin = nthin, reset = TRUE)
  
  # Retrieve samples
  samples <- as.matrix(CphotoMCMC$mvSamples)
  samples2 <- as.matrix(CphotoMCMC$mvSamples2)
  
  # Save samples
  saveRDS(samples, paste0("models/photo_model_", 
                          Sys.Date(), "_samples_ch", i, ".rds"))
  saveRDS(samples2, paste0("models/photo_model_", 
                           Sys.Date(), "_samples2_ch", i, ".rds"))
  # End timer
  mcmc_end[[i]] <- Sys.time()
  # Total time:
  print(mcmc_time[[i]] <- mcmc_end[[i]] - mcmc_start[[i]])
}
