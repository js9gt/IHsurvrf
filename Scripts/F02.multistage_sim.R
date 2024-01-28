
source("F01.one_stage_sim.R")
#~/survrf/Scripts/F01.one_stage_sim.R


## fit propensity model on the baseline covariates
##  propensity_model <- glm(Ak ~ Sk + nstages + cumulative_length, family = binomial(link = "logit"), data = patient_data)


# Function to simulate baseline covariates for stage 0
simulate_baseline_covariates <- function(n.sample) {
  state <- matrix(rnorm(n.sample, mean = 2, sd = 1), ncol = 1)
  action <- matrix( replicate(n.sample, sample(c(0, 1), size = 1, replace = TRUE, prob = c(0.5, 0.5))), ncol = 1)
  return( data.frame(state, action))
}

# Function to simulate data for 10 patients with specified conditions
## test
# n.sample = 10
# max_stages = 10



simulate_patients <- function(n.sample, max_stages, tau,
                              ## initially all patients are at risk
                              at.risk = rep(1, n.sample),
                              ## Life so far lived at the beginning of the stage
                              cumulative.time = rep(0, n.sample)) {
  
  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(at.risk) == 1) at.risk = rep(at.risk, n.sample)
  
  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(cumulative.time) == 1) cumulative.time = rep(cumulative.time, n.sample)
  
  require(dplyr)
  
  
  ## generating initial baseline state
  ##### TO DO: we can make this p-dimentional later
  ## each patient has a randomly generated baseline covariate of State (S_0) and action (A_0), this is the same across all stages)
  baseline_covariates <- simulate_baseline_covariates(n.sample)
  
  ## fitting initial propensity score model from the baseline data
  propensity_model <- glm(action ~ state, family = binomial(link = "logit"), data = baseline_covariates)
  
  
  colnames(baseline_covariates) = c("S_0", "A_0")
  
  # skeleton
  ## uses the dynamics.vec() function to create an initial structure to set up an initial state for 
  ## later dynamics of the multi-stage
  tmp <- 
    one_stage.vec(nstages = 0, cumulative_length = 0, at.risk = 1,
                  propensity_model = propensity_model, terminal.stage = F,
                  a1 = -5, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, r1 = -0.8,
                  a2 = -3, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, r2 = -1, tau = tau) %>% t
  
  
  
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
  output <-  array(NA, dim = c(n.sample, 2 + dim(tmp)[2] + 2 + 2, max_stages), 
                   ## row and column names for the array
                   dimnames = list(1:n.sample, c("subj.id", "rep.id", colnames(tmp),"at.risk", "cumulative.time",
                                                 colnames(baseline_covariates)), 
                                   1:max_stages))
  
  
  ####### COMPARING #####
  #hc_temp <- dynamics.vec(time.max = 1, tick = 0.01, rho.plus = rep(0.5, 2), omega.plus = rep(0.5, 2), 
  #             pred.hazard = 0, pred.censor = 0, admin.censor = Inf,
  #             at.risk = 1, corr = -0.5, full_data = ifelse(FALSE, 2, 0)) %>% t
  
  
  ## initializes an array called "output" filled with NA values with dimensions:
  ##        number of samples being simulated (input)
  ##        dim(tmp)[2]: number of columns in the tmp opject (AKA the output variables from dynamics.vec)
  ##        n.stages: number of stages (input) in simulation
  #hc_output <- 
  #  array(NA, dim = c(n.sample, 2 + dim(hc_temp)[2] + 6 + p, n.stages), 
  #        ## row and column names for the array
  #        dimnames = list(1:n.sample, c("subj.id", "rep.id", colnames(hc_temp), "action", "at.risk", 
  #                                      "surv.previous", "lB", "rho.0", "omega.0", colnames(covariate)), 
  #                        1:n.stages))
  
  ######################
  
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
    output[, colnames(baseline_covariates), stage] <- as.matrix(baseline_covariates)
    
    
    stage.output <- 
      
      ## generate values for current stage
      one_stage.vec(nstages = stage - 1, cumulative_length = cumulative.time, at.risk = at.risk,
                    propensity_model = propensity_model, 
                    terminal.stage = (stage == max_stages),
                    a1 = -7, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, r1 = -0.8,
                    a2 = -0.1, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, r2 = -1, tau = tau) %>% t
    
    ## combining data from the stage.output with "action" and "at.risk" into a matrix format
    ##      combined data includes info on dynamics of simulation, actions, status of individuals
    output[, c(colnames(tmp), "at.risk"), stage] <- cbind(stage.output, at.risk)
    
    ## whether a patient is at risk is updated based on values of gamma and delta. gamma = 0 means at risk
    # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
    ## once delta is 0, the patient is no longer at risk also, only those who are non-censored (delta == 1) are at are risk
    at.risk <- (stage.output[, "gamma"] == 0 & stage.output[, "delta"] == 1)  
    
    ## replaces NA values in at.risk with 0 (not at risk)
    at.risk[is.na(at.risk)] = 0
    
    ## once delta is 0, the patient is no longer at risk also, only those who are non-censored (delta == 1) are at are risk
    #at.risk <- (stage.output[, "delta"] == 1)
    
    ## replaces NA values in at.risk with 0 (not at risk)
    #at.risk[is.na(at.risk)] = 0
    
    ## update the cumulative to prepare for next state calculations by adding the current stage's event time
    cumulative.time <- cumulative.time + stage.output[, "event.time"] # updating for the next stage
  }
  
  return(output)
}


set.seed(123)
simulate_patients(n.sample=4, max_stages = 10, tau = 100,
                  ## initially all patients are at risk
                  at.risk = 1,
                  ## Life so far lived at the beginning of the stage
                  cumulative.time = 0) 


#########################################

## starting from stage 1
for (stage in 1:max_stages) {
  patient_data <- data.frame(
    Stage = integer(),
    Ak = integer(),
    failure.time = numeric(),
    treatment.time = numeric(),
    event.time = numeric(),
    Gamma = integer(),
    nstages = integer(),
    cumulative_length = numeric()
  )
  
  ## looping through all patients
  for (patient_index in 1:n.sample) {
    
    ## initialize cumulative_length variable and number of stages variable
    cumulative_length <- rep(0, n.sample)
    nstages <- rep(0, n.sample)
    
    # Simulate data for each stage using dynamics.vec
    dynamics_output <- dynamics(
      nstages = nstages,
      cumulative_length = cumulative_length,
      propensity_model = propensity_model,
      at.risk = 1,
      terminal.stage = (stage == max_stages),
      a1 = -0.3, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, r1 = -0.8,
      a2 = 1.2, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, r2 = -1
    )
    
    # Extract relevant information from dynamics output
    stage_data <- data.frame(
      Stage = stage,
      Ak = ifelse(nstages == 0, dynamics_output$statistics[3], NA),
      failure.time = ifelse(nstages == 0, dynamics_output$times[1], NA),
      treatment.time = ifelse(nstages == 0, dynamics_output$times[2], NA),
      event.time = ifelse(nstages == 0, dynamics_output$statistics[1], NA),
      Gamma = dynamics_output$statistics[2],
      nstages = nstages,
      cumulative_length = cumulative_length
    )
    
    patient_data <- rbind(patient_data, stage_data)
    
    # Check conditions to break the loop
    ## if the cumulative event time is greater than tau, break
    ## if gamma = 1 (subject fails during this visit)
    ## if we have reached the maximum number of stages, break
    if (dynamics_output$statistics[2] > tau || dynamics_output$statistics[2] == 1 || nstages >= max_stages) {
      break
    }
    
    ## update the value for cumulative length
    cumulative_length <- cumulative_length + dynamics_output$statistics[1]
    
    ## update the number of stages
    nstages <- nstages + 1
  }
  
  all_patient_data[[stage]] <- patient_data  # Store data for each stage
}

return(all_patient_data)
}

# Example: Simulate data for 10 patients with a maximum of 5 stages and tau = 1000
simulated_data <- simulate_patients(n = 10, max_stages = 5, tau = 1000)