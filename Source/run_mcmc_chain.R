#### Function to run single MCMC chain ####
run_mcmc_chain <- function(chain_id,
                           niter = 5e4) {
  
  ## Create nimble model for this chain
  mod <- nimbleModel(code = mod_code, 
                     constants = constants,
                     data = data_list, 
                     inits = inits_fun(data_list, constants, mod_type),
                     calculate = T)
  
  # specify parameters to monitors
  monitor_params = c("D", "N", "psi", "mu.p0", "sex.trap", 
                     "sigma.trap", "sigma", "p.sex", 
                     "D_proj")
  
  if (mod_type == "rw"){
    monitor_params = c(monitor_params,
                       "psi_sd",
                       "logit_psi_proj",
                       "psi_proj")
  }
  
  if (mod_type == "trend"){
    monitor_params = c(monitor_params,
                       "beta0",
                       "beta1",
                       "lambda_proj")
  } 
  
  # Configure MCMC
  conf_mod <- configureMCMC(mod, 
                            monitors = monitor_params,
                            enableWAIC = TRUE)
  
  # Add slice samplers for bounded priors
  slice_pars <- c("mu.p0")
  
  if(mod_type == "trend"){
    slice_pars <- c(slice_pars, "beta0")
  }
  
  for(p in slice_pars){
    conf_mod$removeSampler(p)
    conf_mod$addSampler(target = p, type = "slice")
  }
  
  # 1. Adaptive RW samplers (adjust proposal variance automatically)
  conf_mod$addSampler(target = "sex.trap", type = "RW",
                      control = list(adaptive = TRUE), silent = T)
  
  
  # Add block samplers for S
  for(t in 1:constants$nt){
    for(i in 1:constants$M){
      nodes <- c(paste0("S[", i, ", ", t, ", ", "1]"), paste0("S[", i, ", ", t, ", ", "2]"))
      conf_mod$removeSamplers(nodes)
      conf_mod$addSampler(target = nodes, type = "RW_block", silent = T)
    }
  }
  
  # Build and compile
  mcmc_build <- buildMCMC(conf_mod)
  Cmod <- compileNimble(mod)
  Cmcmc <- compileNimble(mcmc_build, project = Cmod)
  
  
  # Run MCMC
  start_time <- Sys.time()
  cat("Starting chain", chain_id, "\n")
  
  out <- try(runMCMC(Cmcmc, 
                     niter = niter, 
                     nchains = 1,
                     inits = inits_fun(data_list, constants, mod_type),
                     thin = 1,
                     samplesAsCodaMCMC = F,
                     WAIC = T))
  
  end_time <- Sys.time()
  run_time <- round(difftime(end_time, start_time, units = "hours"), 1)
  cat("Chain", chain_id, "completed in", 
      run_time, "hours\n")
  
  return(list(samples = out, chain_id = chain_id, run_time = run_time))
}
