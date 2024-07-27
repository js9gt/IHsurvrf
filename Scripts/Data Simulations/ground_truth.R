

### we want to simulate the truth. Let's start with a patient with 3 stages.
## these are their possible treatment combinations:

## 1. stage 1: 0; stage 2: 0, stage 3: 0
## 2. stage 1: 1; stage 2: 0, stage 3: 0
## 3. stage 1: 0, stage 2: 1, stage 3: 0
## 4. stage 1: 0, stage 2: 0, stage 3: 1
## 5. stage 1: 1, stage 2: 1, stage 3: 1
## 6: stage 1: 0, stage 2: 1, stage 3: 1
## 7: stage 1: 1, stage 2: 0, stage 3: 1
## 8: stage 1: 1, stage 1: 1, stage 3: 0

### we can simulate a list of all these treatment combinations in the following manner:
treatment_combinations <- expand.grid(stage1 = c(0, 1), stage2 = c(0, 1), stage3 = c(0, 1))



simulate_truth <- function(..., n.sample, max_stages, tau,
                           ## initially all patients are at risk
                           at.risk = rep(1, n.sample),
                           ## Life so far lived at the beginning of the stage
                           cumulative.time = rep(0, n.sample),
                           
                           ## inital length of previous visits
                           prior.visit.length = rep(0, n.sample),
                           ## dimensions of the state vector generated at each stage-- 20MAR2024 this is a global variable now
                           #p = 1,
                           a1 = -2, b1 = 0.1, c1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
                           a2 = -1, b2 = -0.05, c2, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
                           a3 = -4.5, b3 = -1, c3, z3 = -2.5, p3 = -0.1, g3 = 0.2, h3 = -0.6, r3 = -1,
                           ## "covariate" value for later state generation
                           rho = 0.5,
                           # action-specific effects
                           D0 = 0,
                           D1 = 1,
                           g = 1,
                           
                           ## a logical for if we want to include censoring (besides administrative censoring)
                           censoringyesno = TRUE,
                           
                           summary = TRUE,
                           treatment_combo = c(0, 0, 0)) {
  
  
  ## if the stage.start is input,we generate that number of stages instead of max_stages (AKA change value of max_stages)
  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(at.risk) == 1) at.risk = rep(at.risk, n.sample)
  
  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(cumulative.time) == 1) cumulative.time = rep(cumulative.time, n.sample)
  
  ## if prior visit length is 0, then replicate the value until it matches the length of n.sample
  if (length(prior.visit.length) == 1) prior.visit.length = rep(prior.visit.length, n.sample)
  
  require(dplyr)
  
  
  
  # skeleton
  ## uses the dynamics.vec() function to create an initial structure to set up an initial state for
  ## later dynamics of the multi-stage
  tmp <- one_stage.vec(nstages = 0, cumulative.length = 0, prior.visit.length = 0, at.risk = 1, time.max = 1000,
                       ## using any action for now just to get column names for the structure
                       terminal.stage = F, input.state.value = 0, input.state.value2 = 0,
                       a1 = -4.5, b1 = -1, c1 = -0.5, z1 = -0.07, p1 = -0.05, g1 = 0.7, h1 = -0.2, r1 = -0.05,
                       ## 2 for time to next visit
                       a2 = -3.9, b2 = -1, c2 = -0.5, z2 = -0.008, p2 = -0.01, g2 = -0.7, h2 = -0.2, r2 = -0.005,
                       a3 = -4.5, b3 = -2, c3 = -0.5, z3 = -0.1, p3 = -0.1, g3 = 0.2, h3 = -0.4, r3 = -0.05,
                       tau = tau, censoringyesno = TRUE, input.policy.action = NA, input_opt = NA) %>% t
  
  
  
  
  ## initializes an array called "output" filled with NA values with dimensions:
  ##        number of samples being simulated (input)
  ##        dim(tmp)[2]: number of columns in the tmp opject (AKA the output variables from dynamics.vec)
  ##        n.stages: number of stages (input) in simulation
  output <-  array(NA, dim = c(n.sample, 2 + 2 + dim(tmp)[2] + 2 + 2 + 2, max_stages),
                   ## row and column names for the array
                   dimnames = list(1:n.sample, c("subj.id", "rep.id", "baseline1", "baseline2", colnames(tmp), "at.risk",
                                                 ## new column to track number of each treatment
                                                 "action.1.count",
                                                 "action.0.count",
                                                 "cumulative.time", "IHsurvrf.A",
                                                 "trt.diff"),
                                   #colnames(baseline_covariates)),
                                   1:max_stages))
  
  
  ## initializes the subj.id column we just created in the output from 1 to n.sample
  output[, "subj.id", ] <- 1:n.sample
  
  ## initialize rep.id column with 0s
  output[, "rep.id", ] <- 0           # rep.id is reserved for later use (repeated random draws)
  
  
  ## initialize an input.state vector of 0
  ## previous was runif()
  input.state1 <- rnorm(n.sample, 0, 0.5)
  input.state2 <- rnorm(n.sample, 0, 0.5)
  
  ## initialize action counts as 0's
  action.1.count <- rep(0, n.sample)
  action.0.count <- rep(0, n.sample)
  
  ### note: we need to generate all the covariates used in the model before inputting these into the policy
  
  
  ## initialize a value of baseline that's drawn from a U(0, 1) distribution
  output[, "baseline1", ] <- runif(n.sample,0, 1)
  
  ## initialize a second with random draw from a binomial distribution
  output[, "baseline2", ] <- as.factor(rbinom(n.sample, size = 1, prob = 0.5))
  
  ## initialize time.max as a vector with the value being tau at the first stage
  time.max <- rep(tau,n.sample)
  
  ## initialize an input action vector with NA values, with the same length as number of samples
  input.policy.action = rep(NA, n.sample)
  
  
  ## iterate through each stage
  for (stage in 1:max_stages) {
    
    ## print the current stage number
    cat("stage ", stage, "\n")
    
    
    ## update the at.risk column in the output array for the current stage using the value from at.risk
    output[, "at.risk", stage]           <- at.risk
    
    
    ## update the cumulative time column in the output array for the current stage using the cumulative vector
    ## to avoid it being overly large, have this be a proportion of tau
    output[, "cumulative.time", stage]     <- cumulative.time
    
    ## calculates maximum time allowed in simulation based on diff between tau (limit) and existing time lived
    ## this should be updated at each stage: recall that cumulative.time is a proportion of tau
    ## initialize this as the input
    output[, "time.max", stage] = time.max
    
    ## initialize the counts of actions with 0's
    output[,"action.1.count", stage] <- action.1.count
    output[, "action.0.count", stage] <- action.0.count
    
    ## initialize state, prior.visit.length, nstages as 0
    output[, "state1", stage] <- input.state1
    output[, "state2", stage] <- input.state2
    output[, "prior.visit.length", stage] <- prior.visit.length
    output[, "nstages", stage] <- rep(0,n.sample)
    
    
    input.policy.action[at.risk != 0] <- treatment_combo[, stage]
    
    ### generate observed data following the input treatment sequence
    ## for each stage, we index another one
    stage.output <-
      
      ## generate values for current stage
      ## for nstages, to avoid the value being overly large, have it also be a proportion of the total number of stages
      one_stage.vec(nstages = (stage - 1)/max_stages, cumulative.length = cumulative.time, at.risk = at.risk,time.max = time.max,
                    prior.visit.length = prior.visit.length,
                    #action = as.numeric(output[, "action", stage]),  ### remove this after clarification
                    terminal.stage = (stage == max_stages), input.state.value = input.state1, input.state.value2 = input.state2,
                    a1 = a1, b1 = b1, c1 = c1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
                    a2 = a2, b2 = b2, c2 = c2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
                    a3 = a3, b3 = b3, c3 = c3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau, censoringyesno = censoringyesno,
                    
                    ## inputting action == 1
                    input.policy.action = input.policy.action) %>% t
    
    ## for people who aren't at risk anymore, this is set to NA
    input.policy.action[at.risk == 0] <- NA
    
    ## combining data from the stage.output with "action" and "at.risk" into a matrix format
    ##      combined data includes info on dynamics of simulation, actions, status of individuals
    output[, c(colnames(tmp),"at.risk"), stage] <- cbind(stage.output, at.risk)
    
    ## whether a patient is at risk is updated based on values of gamma and delta. gamma = 0 means at risk
    # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
    ## once delta is 0, the patient is no longer at risk also, only those who are non-censored (delta == 1) are at are risk
    at.risk <- (stage.output[, "gamma"] == 0 & stage.output[, "delta"] == 1)
    
    
    ## replaces NA values in at.risk with 0 (not at risk)
    at.risk[is.na(at.risk)] = 0
    
    
    ## update the cumulative to prepare for next state calculations by adding the current stage's event time
    ## to avoid it being overly large, have this be a proportion of tau
    cumulative.time <- (cumulative.time*tau + stage.output[, "event.time"]) / tau # updating for the next stage
    
    
    time.max <- tau - (cumulative.time*tau)
    
    ## update the prior visit length column in the output array for the previous stage
    prior.visit.length    <- stage.output[, "event.time"]
    
    # Update action count vectors for action 1
    action.1.count <- action.1.count + (stage.output[, "action"] == 1)
    
    # Update action count vectors for action 0
    action.0.count <- action.0.count + (stage.output[, "action"] == 0)
    
    # if it's not the initial stage, input the value from the updated state (using the formulas)
    # this will use values from the previous state
    # NOTE: this lags behind by one iteration, that's why we can update the prior.visit.length like we do above
    
    ## we have first D1 term with the I() be opposite effects to find optimal rule immediately
    input.state1 = (rho*stage.output[, "state1"] + ifelse(stage.output[, "state1"] > 0, 1, 0)*(D0 - stage.output[, "action"]*D1) +
                      ifelse(stage.output[, "state1"] < 0, 1, 0)*(D0^2 + stage.output[, "action"]*(D1)^2) +
                      0.5*g*(1-(rho)^2)^0.5 * replicate(n.sample, rnorm(n = 1, 0, 0.1))-
                      
                      ## subtracting a constant value to get an average
                      0)
    
    
    input.state2 = (rho*stage.output[, "state2"] + ifelse(stage.output[, "state2"] > 0, 1, 0)*(D0 - stage.output[, "action"]*D1) +
                      ifelse(stage.output[, "state2"] < 0.5, 1, 0)*(D0^2 + stage.output[, "action"]*(D1)^2) +
                      0.5*g*(1-(rho)^2)^0.5 * replicate(n.sample, rnorm(n = 1, 0, 0.1))-
                      
                      ## subtracting a constant value to get an average
                      0)
    
    
  }
  
  
  if (summary) {
    
    ## selecting only stages 1:ss
    ## NOTE if ss input is null, it's set to the maximum number of stages
    cum.event.time = apply(output[, "event.time",1:max_stages], 1, cumsum) %>% t
    output.summary <-
      tibble(subj.id = output[, "subj.id", 1],
             rep.id = output[, "rep.id", 1],
             terminal.stage = apply(output[, "at.risk", 1:max_stages], 1, sum),
             cumulative.event.time = cum.event.time[cbind(1:n.sample, terminal.stage)],
             censor.status = output[cbind(1:n.sample, "delta", terminal.stage)],
             actions       = apply(output[, "action", ], 1,
                                   function(s) gsub("NA", "*", paste0(s, collapse = ""))))
    return(list(output = output, summary = output.summary))
  } else {
    return(output)
  }
  
}


### now, we want to create a for loop, where we loop through each row of the possible treatment combinations
### we simulate patients according to those.
### then, we pick the best regimen (for each patient), as the one that gives the longest cumulative time
### we pick the worst regimen (for each patient), as the one that gets the shortest cumulative time


# Initialize an empty list to store the results
results <- list()

# Loop through each row of treatment_combinations and call simulate_truth, storing the results
for (i in 1:nrow(treatment_combinations)) {
  result <- simulate_truth(n.sample = 5, max_stages = 3, tau = 1000,
                           at.risk = rep(1, 5),
                           cumulative.time = rep(0, 5),
                           prior.visit.length = rep(0, 5),
                           a1 = -2, b1 = 0.1, c1 = -0.5, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
                           a2 = -1, b2 = -0.05, c2 = -0.5, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
                           a3 = -4.5, b3 = -1, c3 = -0.5, z3 = -2.5, p3 = -0.1, g3 = 0.2, h3 = -0.6, r3 = -1,
                           rho = 0.5,
                           D0 = 0,
                           D1 = 1,
                           g = 1,
                           censoringyesno = TRUE,
                           summary = TRUE,
                           treatment_combo = treatment_combinations[i, ])
  results[[i]] <- result$summary
}



# Initialize lists to store the indices of the maximum and minimum cumulative.event.time for each patient
max_index <- rep(NA, length(results[[1]]$cumulative.event.time))
min_index <- rep(NA, length(results[[1]]$cumulative.event.time))

# Initialize lists to store the maximum and minimum cumulative.event.time for each patient
max_cumulative_event_time <- rep(-Inf, length(results[[1]]$cumulative.event.time))
min_cumulative_event_time <- rep(Inf, length(results[[1]]$cumulative.event.time))

# Loop through the results and update the indices of the maximum and minimum values for each patient
for (i in 1:length(results)) {
  result <- results[[i]]
  for (j in 1:length(result$cumulative.event.time)) {
    if (result$cumulative.event.time[j] > max_cumulative_event_time[j]) {
      max_cumulative_event_time[j] <- result$cumulative.event.time[j]
      max_index[j] <- i
      
      max_cumulative_event_time[j] <- max(max_cumulative_event_time[j], result$cumulative.event.time[j])
    }
    if (result$cumulative.event.time[j] < min_cumulative_event_time[j]) {
      min_cumulative_event_time[j] <- result$cumulative.event.time[j]
      min_cumulative_event_time[j] <- min(min_cumulative_event_time[j], result$cumulative.event.time[j])
      min_index[j] <- i
    }
  }
}

# Combine the results into a data frame for easier viewing
comparison_df <- data.frame(
  patient = 1:length(max_index),
  max_index = max_index,
  max_time = max_cumulative_event_time,
  min_index = min_index,
  min_time = min_cumulative_event_time
)

comparison_df$best_seq <- treatment_combinations[max_index, ]
comparison_df$worst_seq <- treatment_combinations[min_index, ]

