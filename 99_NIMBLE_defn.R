###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#----------------------------NIMBLE Model Definition---------------------------X
###############################################################################X

library(nimble)

# Model definition----
photoCode <- nimbleCode({
  # ... priors ----
  # Abundance | Occupancy
  for (i in 1:9){
    beta[i] ~ dnorm(mean = 0, sd = 10)
  }
  
  # Occupancy
  for (i in 1:7){
    gamma[i] ~ dnorm(mean = 0, sd = 10)
  }
  
  # Site-level REs
  sigma_eta_psi ~ T(dnorm(mean = 1, sd = 1), 0, )
  sigma_eta_w ~ T(dnorm(mean = 1, sd = 1), 0, )
  
  # ... Likelihood ----
  
  # ... ... Process model ----
  for (i in 1:S){
    # Site-level random effect for occupancy
    eta_psi[i] ~ dnorm(0, sd = sigma_eta_psi)
    # Site-level random effect for abundance|occupancy
    eta_w[i] ~ dnorm(0, sd = sigma_eta_w)
    
    for (t in 1:max_T){
      # Occupancy
      psi[i, t] <- ilogit(gamma[1] + 
                            gamma[2] * lvstk[i, t] +
                            gamma[3] * snow[i, t] +
                            gamma[4] * elev_doy[i, t] +
                            gamma[5] * elev_doy2[i, t] +
                            gamma[6] * doy[t] +
                            gamma[7] * (doy[t])^2 +
                            eta_psi[i])
      #Sample occupancy
      z[i, t] ~ dbern(psi[i, t])
      
      # E[N] = lambda
      lambda[i, t] <- 
        # Occupancy
        z[i, t] *
        # Abundance if occupied
        exp(beta[1] + 
              beta[2] * doy[t] +
              beta[3] * (doy[t])^2 +
              beta[4] * grass_forb[i] +
              beta[5] * ndvi[i, t] +
              beta[6] * irg[i, t] + 
              beta[7] * ndvi_grass_forb[i, t] +
              beta[8] * irg_grass_forb[i, t] +
              beta[9] * tree[i] +
              eta_w[i])
    }
      
      for (t in 1:max_T){
      # Likelihood of count
      y[i, t] ~ dpois(lambda = lambda[i, t])
    }
  }
})