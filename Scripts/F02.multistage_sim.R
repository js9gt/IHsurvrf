
source("F01.one_stage_sim.R")
#~/survrf/Scripts/F01.one_stage_sim.R




# Function to simulate data for 10 patients with specified conditions
## test
# n.sample = 10
# max_stages = 5
# tau = 100
# p = 1




simulate_patients <- function(n.sample, max_stages, tau,
                              ## initially all patients are at risk
                              at.risk = rep(1, n.sample),
                              ## Life so far lived at the beginning of the stage
                              cumulative.time = rep(0, n.sample),
                              
                              ## inital length of previous visits 
                              prior.visit.length = rep(0, n.sample),
                              ## dimensions of the state vector generated at each stage
                              p = 1) {
  
  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(at.risk) == 1) at.risk = rep(at.risk, n.sample)
  
  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(cumulative.time) == 1) cumulative.time = rep(cumulative.time, n.sample)
  
  ## if prior visit length is 0, then replicate the value until it matches the length of n.sample
  if (length(prior.visit.length) == 1) prior.visit.length = rep(prior.visit.length, n.sample)
  
  require(dplyr)
  
  ######################################### TO DELETE AFTER CLARIFICATION
  ## generating initial baseline state
  ## each patient has a randomly generated baseline covariate of State (S_0)
  #baseline_covariates <- simulate_baseline_covariates(n.sample, p)
  
  #colnames(baseline_covariates) = c("S_0")
  
  
  ## generates predicted probability of treatment for each patient as a vector
  #propensity_scores <- propensity_function(state = baseline_covariates, p)
  
  ########################################## 

  
  # skeleton
  ## uses the dynamics.vec() function to create an initial structure to set up an initial state for 
  ## later dynamics of the multi-stage
  tmp <- 
    one_stage.vec(nstages = 0, cumulative_length = 0, prior.visit.length = 0, at.risk = 1,
                  ## using any action for now just to get column names for the structure
                  terminal.stage = F,
                  a1 = -2, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
                  a2 = -1, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1, tau = tau, p = p) %>% t
  
  
  
  ### test
  #n.sample = 10
  
  # n.stages for HC, max_stages for JS
  #n.stages = 5
  
  #max_stages = 5
  # number of covariates
  #p = 6
  
  ####### COMPARING #####
  # HC generate covariate matrix
  #covariate = mvrnorm(n.sample, rep(0, p), diag(p) + 0.2 - diag(p) * 0.2)
  #colnames(covariate) = paste0("Z", 1:p)
  ###
  ######################
  
  ## initializes an array called "output" filled with NA values with dimensions:
  ##        number of samples being simulated (input)
  ##        dim(tmp)[2]: number of columns in the tmp opject (AKA the output variables from dynamics.vec)
  ##        n.stages: number of stages (input) in simulation
  output <-  array(NA, dim = c(n.sample, 2 + dim(tmp)[2] + 2, max_stages), 
                   ## row and column names for the array
                   dimnames = list(1:n.sample, c("subj.id", "rep.id", colnames(tmp), "at.risk", "cumulative.time"),
                                                 #colnames(baseline_covariates)), 
                                   1:max_stages))
  
  
  ## initializes the subj.id column we just created in the output from 1 to n.sample
  output[, "subj.id", ] <- 1:n.sample
  
  ## initialize rep.id column with 0s
  output[, "rep.id", ] <- 0           # rep.id is reserved for later use (repeated random draws)
  
  
  ## iterate through each stage
  for (stage in 1:max_stages) {
    
    ## print the current stage number
    cat("stage ", stage, "\n")
    
    
    ## update the at.risk column in the output array for the current stage using the value from at.risk
    output[, "at.risk", stage]           <- at.risk
    
    
    ## update the cumulative time column in the output array for the current stage using the cumulative vector
    output[, "cumulative.time", stage]     <- cumulative.time
    
    ## updates corresponding matching column names of the baseline covariates with the values in there
    #output[, colnames(baseline_covariates), stage] <- as.matrix(baseline_covariates)
    
    
    
    stage.output <- 
      
      ## generate values for current stage
      ##### why is this outputting 10 values of each of the rates instead of one?
      one_stage.vec(nstages = stage - 1, cumulative_length = cumulative.time, at.risk = at.risk, 
                    prior.visit.length = prior.visit.length,
                    #action = as.numeric(output[, "action", stage]),  ### remove this after clarification
                    terminal.stage = (stage == max_stages),
                    a1 = -2, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
                    a2 = -1, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1, tau = tau) %>% t
    
    
      
      ## generate values for current stage
      #one_stage.vec(nstages = 2, 
      #              # each patient's cumulative time is 0
      #              cumulative_length = cumulative.time, 
      #              # each patient is at.risk at first
      #              at.risk = at.risk,
      #              ## for action vector use the one generated in the "action" column earlier
      #              ## these are different actions FOR EACH PATIENT, so we need to pull the values across each stage
      #              action = rep(1, n.sample),  
      #              terminal.stage = F,
      #              a1 = -7, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, r1 = -0.8,
      #              a2 = -0.1, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, r2 = -1, tau = tau) %>% t
    

    
    ## combining data from the stage.output with "action" and "at.risk" into a matrix format
    ##      combined data includes info on dynamics of simulation, actions, status of individuals
    output[, c(colnames(tmp),"at.risk"), stage] <- cbind(stage.output, at.risk)
    
    ## whether a patient is at risk is updated based on values of gamma and delta. gamma = 0 means at risk
    # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
    ## once delta is 0, the patient is no longer at risk also, only those who are non-censored (delta == 1) are at are risk
    at.risk <- (stage.output[, "gamma"] == 0 & stage.output[, "delta"] == 1)  
    
    ## replaces NA values in at.risk with 0 (not at risk)
    at.risk[is.na(at.risk)] = 0
    
    ##
    ### replace not at risk with NA in action 
    ##
    ##
    
    ## update the cumulative to prepare for next state calculations by adding the current stage's event time
    cumulative.time <- cumulative.time + stage.output[, "event.time"] # updating for the next stage
    
    ## update the prior visit length column in the output array for the previous stage
    prior.visit.length    <- stage.output[, "event.time"]
  }
  
  return(output)
}


set.seed(123)
simulate_patients(n.sample=10, max_stages = 3 , tau = 1000,
                  ## initially all patients are at risk
                  at.risk = 1,
                  ## Life so far lived at the beginning of the stage
                  cumulative.time = 0,
                  prior.visit.length = 0, 
                  p = 1) 

