library(nimble)
library(nimbleSCR)
library(parallel)
library(ggplot2)
library(coda)

#### Helper functions ####
### Data prep function (returns data and constants
source("../Source/prep_panthera_data.R")

### Inits function
source("../Source/init_funs.R")

### Function to fit chain (used for parallel model fitting)
source("../Source/run_mcmc_chain.R")

#### Model Inputs ####
### Specify site name
site = "Kafue_South"

### Main data list
load("../Data/Rogan_data_2025.Rdata")

inputs <- prep_panthera_data(secr_cptr_fls, site, max_density = 14)

data_list <- inputs$data
constants <- inputs$constants

### Specify model specification ("static", "rw", "trend")
mod_type = "trend"

# derive upper limit of population size on log scale from M
if(mod_type == "trend"){
  constants$bmax <- log(constants$M)
  constants$time <- (1:constants$nt) - mean(1:constants$nt)
  constants$future <- constants$time[constants$nt] + c(1, 2)
  constants$n_future <- 2
} 


### Source model code
source(paste0("../Source/nimble_scr_",
              mod_type,
              "_model.R"))



### Data list
# data_list <- list(
#   Y, # capture histories
#   X, # trap coordinates
#   K, # sampling effort
#   Sex, # partially observed sex
#   Z # partially observed state of being in true population
# )

### Constants

#   M, # Total number of individauls in augmented population for each session
#   nt = 5, # number of sessions
#   nj, # vector of the number of traps in each session 
#   A, # Area in sq km
#   xlim, # horizontal limits of the state space
#   ylim, # vertical limits of the state space

#### Fit chains in parallel ####
# Detect number of cores and use 3 (or fewer if less than 3 available)
n_cores <- min(4, detectCores() - 1)  # Leave one core free
cat("Using", n_cores, "cores for parallel processing\n")

# Set up cluster
cl <- makeCluster(n_cores)

# Export necessary objects to cluster
clusterExport(cl, c("mod_code", "constants", "data_list", "inits_fun", "mod_type"))

# Load required libraries on each worker
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleSCR)
})

# Run 3 chains in parallel
{
  total_run_time <- Sys.time()
  
  results <- parLapply(cl, 1:4, run_mcmc_chain, 
                       niter = 2.25e4)
  
  total_run_time[2] <- Sys.time()
  total_time <- round(difftime(total_run_time[2], 
                               total_run_time[1], 
                               units = "hours"), 
                      1)
  cat("Total parallel execution time:", total_time, "hours\n")
  
  # Stop cluster
  stopCluster(cl)
  
  save.image(paste0("../Outputs/",
                    mod_type,
                    "_model_results_workspace.rdata"))
}


