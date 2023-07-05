# log prior

lnprior = function(tau, theta, time, mu_tau, mu_theta) {
  
  lnp = -(tau / mu_tau) - (theta / mu_theta) # Exponential prior on tau with mean mu_tau
                                             # Exponential prior on theta with mean mu_theta
  
  eq = log(2/theta) - 2 / theta * time
  
  eq = sum(eq)
  
  return(lnp + eq)
}

# log likelihood

lnlikelihood = function(tau, time) {
  
  p = 3/4 - 3/4 * exp(-8/3 * (tau + time))
  
  lnp = data$xi * log(p) + (data$ni - data$xi) * log(1 - p)
  
  lnp = sum(lnp)
  
  return(lnp)
  
}

# Acceptance ratio

ln_accept_ratio_time_j = function(j, tau, theta, time_j, new.time_j){
  
  p = 3/4 - 3/4 * exp(-8/3 * (tau + time_j))
  
  p.update = 3/4 - 3/4 * exp(-8/3 * (tau + new.time_j))
  
  accept.ratio = -2 / theta*(new.time_j - time_j) + data$xi[j] * log(p.update / p) + 
    (data$ni[j] - data$xi[j]) * log((1 - p.update) / (1 - p))
  
  return(accept.ratio)
}


# Bayesian Markov Chain Monte Carlo implementation of the multispecies coalescent model for 2 species 

mcmc_MSC = function(N, tau, theta, mu_tau, mu_theta, time_0, w_tau, w_theta, w_t, nloci = 1000){
  
  run.time = Sys.time() 
  
  time = rep(time_0, nloci)    # Create time vector
  
  sample_tau = sample_theta = numeric(N+1)        # Initialise
  accept.tau = accept.theta = accept.time = 0  
  sample_tau[1] = tau
  sample_theta[1] = theta  
  lnp = lnprior(tau, theta, time, mu_tau, mu_theta)  
  L = lnlikelihood(tau, time) 
  
  for(i in 1:N){
    
    tau.new = tau + w_tau*(runif(1)-0.5)
    if (tau.new<0) tau.new = -tau.new; 
    
    lnp.new = lnprior(tau.new, theta, time, mu_tau, mu_theta)        # Update prior with new tau
    L.new = lnlikelihood(tau.new, time)                              # Update likelihood with new tau 
    lnaccept = lnp.new + L.new - lnp - L 
    
    if(lnaccept >= 0 || runif(1) < exp(lnaccept)) {  # If proposal meets requirements:
      
      tau = tau.new                                  # - Update tau
      lnp = lnp.new                                  # - Update log prior 
      L = L.new                                      # - Update Likelihood 
      accept.tau = accept.tau + 1                    # - Add 1 to acceptance vector 
      
    }
    
    theta.new = theta + w_theta*(runif(1)-0.5)
    if (theta.new<0) theta.new = -theta.new
    
    lnp.new = lnprior(tau, theta.new, time, mu_tau, mu_theta)  # Only update the log prior because theta is not a present 
                                                               # in the log likelihood calculation. 
    lnaccept = lnp.new - lnp     
    
    if(lnaccept >= 0 || runif(1) < exp(lnaccept)) {
      
      theta = theta.new
      lnp = lnp.new 
      accept.theta = accept.theta + 1
      
    }
    
    for(j in 1:nloci){
      
      time.j = time[j]    # Store time[j] in time.j. Makes it easier to keep track of everything
      
      time.j.new = abs(time.j + w_t * (0.5 - runif(1))) # Update time.j with new proposal 
      lnaccept = ln_accept_ratio_time_j(j, tau, theta, time.j, time.j.new)
      
      if(lnaccept >= 0 || runif(1) < exp(lnaccept)){
        
        time[j] = time.j.new
        accept.time = accept.time + 1
        
      } 
    }
    
    sample_tau[i+1] = tau             # Add accepted tau estimate
    sample_theta[i+1] = theta         # Add accepted theta estimate
    
    lnp = lnprior(tau, theta, time, mu_tau, mu_theta)   # Update log prior 
    L = lnlikelihood(tau, time)                         # Update log likelihood
    
    
    if(i %% 100 == 0) { # progress to ensure the program is running and not dead
      
      print(paste(i, "iterations done. acceptance:  ", 
                  "tau", accept.tau/i, "  ",
                  "theta", accept.theta/i, "  ",
                  "time", accept.time/(i*nloci)))
      
    }
    
    
  }
  
  accept.tau = accept.tau/N
  accept.theta = accept.theta/N
  accept.time = accept.time/(N * nloci)  # Proportion of accepted time proposals.
  run.time = Sys.time() - run.time         # Total runtime

  return(list(sample_tau, sample_theta, accept.tau, accept.theta, accept.time, run.time))
  
}
