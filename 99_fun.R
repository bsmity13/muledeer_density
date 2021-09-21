###############################################################################X
#-----------Code for "An  Eulerian Perspective on Green-Wave Surfing"----------X
#---------------------T. Del Bosco, B. J. Smith, & T. Avgar--------------------X
#-----------------------------Submitted to Ecology-----------------------------X
#----------------------------Last updated 2021-09-21---------------------------X
###############################################################################X
#-------------------------------Helper Functions-------------------------------X
###############################################################################X

# Function to get matching parameter names from all chains ----
match_samples <- function(pattern, samples, fixed = FALSE,
                          exact = FALSE, thin_total = NULL) {
  
  # Possibly thin each chain to desired length
  if (is.null(thin_total)) {
    r <- 1:nrow(samples[[1]])
  } else {
    r <- seq(1, nrow(samples[[1]]), length.out = thin_total)
  }
  
  if (exact){
    y <- do.call(rbind,
                 lapply(samples, function(x){
                   x[r, pattern, drop = FALSE]
                 }))
  } else {
    y <- do.call(rbind,
                 lapply(samples, function(x){
                   x[r, grep(pattern, colnames(x), fixed = fixed), drop = FALSE]
                 }))
  }
  return(y)
}
