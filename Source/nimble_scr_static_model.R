require(nimble)
require(nimbleSCR)

#### Model code ####
mod_code <- nimbleCode({
  
  ### State parameters
  
  # prior for psi
  psi ~ dbeta(2, 3) 
  
  ### Detection parameters
  
  #prior for probability of being male
  p.sex ~ dbeta(3, 4) # based on known biology
  lp.sex <- logit(p.sex) # convert to logit scale
  
    #priors for baseline detection probability of sexes with random trap effect
  mu.p0 ~ dbeta(1, 10) #beta distribution of the mean of baseline detection probability per day. Based on personal knowledge of trapping 
  lmu.p0 <- logit(mu.p0) # transform to logit scale
  sigma.trap ~ dgamma(0.5, 2) #sigma of the normal distribution for the random trap effect on logit scale 
  sex.trap ~ dnorm(0, sd = 0.5) #half the difference in baseline detection probability between the sexes on the logit
  
  #specify prior for sex-specific sigma
  for(sx in 1:2){
    log_sigma[sx] ~ dnorm(log(3), sd = 0.5) #log-transformed spatial decay parameter in km
    sigma[sx] <- exp(log_sigma[sx]) #convert to real scale
    alpha[sx] <- -1/(2*sigma[sx]^2) #half-normal rate of decay
  }
  
  for(t in 1:nt){ #nt is the number of sessions
    
    # specify trap effects
    for(j in 1:nj[t]){ #nj[t] is the number of traps in session t
      lp0.trap[j, t] ~ dnorm(lmu.p0, sd = sigma.trap) # trap-specific baseline detection probability on logit scale
      logit(p0[j, t, 1]) <- lp0.trap[j, t] - (0.5*sex.trap) # female trap-specific baseline detection rate
      logit(p0[j, t, 2]) <- lp0.trap[j, t] + (0.5*sex.trap) # male trap-specific baseline detection rate
    }
    
    for(i in 1:M){
      Sex[i, t] ~ dbern(p.sex)
      Sex2[i, t] <- Sex[i, t] + 1 #convert to categorical variable
      
      S[i, t, 1] ~ dunif(xlim[1], xlim[2])
      S[i, t, 2] ~ dunif(ylim[1], ylim[2])
      
      Z[i, t] ~ dbern(psi) #realization of probability individual i is in the true population
      
      #Calculate squared distances between individual i's activity center and each trap
      d2[i, 1:nj[t], t] <- (S[i, t, 1] - X[1:nj[t], t, 1])^2 + 
        (S[i, t, 2] - X[1:nj[t], t, 2])^2 
      
      # Expected detection rate
      p[i, 1:nj[t], t] <- p0[1:nj[t], t, Sex2[i, t]] * 
        exp(alpha[Sex2[i, t]] * d2[i, 1:nj[t], t])

      Y[i, 1:nj[t], t] ~ dbinom_vector(size = K[1:nj[t], t], 
                                    prob = p[i, 1:nj[t], t]*Z[i, t]) #observation process
      
    }
    
    N[t] <- sum(Z[1:M, t]) #realized number of individuals in true population at time t
  }
  
  D <- (psi*M)/A * 100 # convert psi to density
  
})
