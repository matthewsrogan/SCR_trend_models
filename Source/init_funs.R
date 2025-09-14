inits_fun <- function(dat, cnts, model){
  
  if(!model %in% c("static", "rw", "trend")){
    stop("The model parameterisation must be specified as 'static', 'rw' for random walk, or 'trend'.")
  }
  
  if(model == "static"){
    inits_list <- list(psi = rbeta(1, 2, 3))
  }
  
  if(model == "rw"){
    inits_list <- list(logit_psi0 = logit(rbeta(1, 2, 3)),
                       psi_sd = rgamma(1, 1, 2))
    
    inits_list$logit_psi <- inits_list$logit_psi0
    
    for(t in 2:cnts$nt){
      inits_list$logit_psi[t] = rnorm(1, inits_list$logit_psi[t-1], inits_list$psi_sd)
    }
    
    inits_list$psi <- secr::invlogit(inits_list$logit_psi)
  }
  
  if(model == "trend"){
    inits_list <- list(beta0 = log(runif(1, 50, cnts$M - 50)), #allow a buffer for population growth
                       beta1 = rnorm(1, 0, 0.5))
    inits_list$lambda <- exp(inits_list$beta0 + 
                               cnts$time * 
                               inits_list$beta1)
    inits_list$psi = inits_list$lambda/cnts$M
    
    # inits for future projections
    inits_list$lambda_proj <- exp(inits_list$beta0 + 
                                    cnts$future *
                                    inits_list$beta1)
    inits_list$D_proj <- inits_list$lambda_proj/cnts$A * 100
  }
  
  # Z inits: randomly assign 0 or 1 to augmented individuals and NA to observed individuals
  z_strt <- matrix(rbinom(cnts$M * cnts$nt, 1, inits_list$psi[1]), 
                   nrow = cnts$M, 
                   ncol = cnts$nt)
  z_strt[is.na(dat$Z) == F] = as.integer(NA)
  
  # Sex inits
  inits_list$p.sex = rbeta(1, 3, 4)
  inits_list$lp.sex = logit(inits_list$p.sex)
  
  
  sex_strt <- matrix(rbinom(cnts$M * cnts$nt, 1, inits_list$p.sex), 
                     nrow = cnts$M,
                     ncol = cnts$nt)
  sex_strt[is.na(dat$Sex) == F] = as.integer(NA)
  
  inits_list$Sex <- sex_strt
  inits_list$Sex2 <- inits_list$Sex + 1
  
  # Random activity center starting locations
  s_strt <- round(array(c(matrix(runif(cnts$M * cnts$nt, cnts$xlim[1], cnts$xlim[2])),
                          matrix(runif(cnts$M * cnts$nt, cnts$ylim[1], cnts$ylim[2]))), 
                        dim = c(cnts$M, cnts$nt, 2)), 2)
  inits_list$S = s_strt
  
  # Detection parameters
  inits_list$log_sigma = rnorm(2, log(3), 0.5)
  inits_list$sigma = exp(inits_list$log_sigma)
  inits_list$alpha =  -1/(2*inits_list$sigma^2)
  
  inits_list$mu.p0 = rbeta(1, 1, 10)
  inits_list$lmu.p0 = logit(inits_list$mu.p0)
  
  inits_list$sigma.trap = rgamma(1, 0.5, 2)
  inits_list$sex.trap = rnorm(1, 0, sd = 0.5)
  
  
  inits_list$lp0.trap = matrix(rnorm(cnts$nt * max(cnts$nj), 
                                     inits_list$lmu.p0, 
                                     inits_list$sigma.trap), 
                               nrow = max(cnts$nj), ncol = cnts$nt)
  inits_list$p0 = array(c(
    inits_list$lp0.trap - 0.5*inits_list$sex.trap,
    inits_list$lp0.trap + 0.5*inits_list$sex.trap
  ), dim = c(max(cnts$nj), cnts$nt, 2))
  
  return(inits_list)
  
}
