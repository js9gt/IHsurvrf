
source("F01.one_stage_sim.R")
source("F03.Figs.R")
#~/survrf/Scripts/F01.one_stage_sim.R
#~/survrf/Scripts/F03.Figs.R



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
                              p = 1,
                              a1 = -2, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
                              a2 = -1, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
                              a3 = -4.5, b3 = -1, z3 = -2.5, p3 = -0.1, g3 = 0.2, h3 = -0.6, r3 = -1,
                              ## "covariate" value for later state generation
                              rho = 0.5,
                              # action-specific effects 
                              D0 = 0,
                              D1 = 1,
                              g = 1) {
  
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
                  terminal.stage = F, initial.stage = T, input.state.value = 0, 
                  a1 = -4.5, b1 = -1, z1 = -0.07, p1 = -0.05, g1 = 0.7, h1 = -0.2, r1 = -0.05,
                  ## 2 for time to next visit
                  a2 = -3.9, b2 = -1, z2 = -0.008, p2 = -0.01, g2 = -0.7, h2 = -0.2, r2 = -0.005, 
                  a3 = -4.5, b3 = -2, z3 = -0.1, p3 = -0.1, g3 = 0.2, h3 = -0.4, r3 = -0.05,
                  tau = tau, p = p) %>% t
  
  
  
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
  
  
  ## initialize an input.state vector of 0
  input.state <- rep(0, n.sample)
  
  ## iterate through each stage
  for (stage in 1:max_stages) {
    
    ## print the current stage number
    cat("stage ", stage, "\n")
    
    
    ## update the at.risk column in the output array for the current stage using the value from at.risk
    output[, "at.risk", stage]           <- at.risk
    
    
    ## update the cumulative time column in the output array for the current stage using the cumulative vector
    ## to avoid it being overly large, have this be a proportion of tau
    output[, "cumulative.time", stage]     <- cumulative.time
    
    ## updates corresponding matching column names of the baseline covariates with the values in there
    #output[, colnames(baseline_covariates), stage] <- as.matrix(baseline_covariates)
    
    
    
    stage.output <- 
      
      ## generate values for current stage
      ## for nstages, to avoid the value being overly large, have it also be a proportion of the total number of stages
      one_stage.vec(nstages = (stage - 1)/max_stages, cumulative_length = cumulative.time, at.risk = at.risk, 
                    prior.visit.length = prior.visit.length,
                    #action = as.numeric(output[, "action", stage]),  ### remove this after clarification
                    terminal.stage = (stage == max_stages), input.state.value = input.state, 
                    ## if initial stage is true (at stage 1), then use uniform(0, 1), otherwise, use the input
                    initial.stage = (stage == 1), 
                    a1 = a1, b1 = b1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
                    a2 = a2, b2 = b2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2, 
                    a3 = a3, b3 = b3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau) %>% t
    
    
      
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
    ## to avoid it being overly large, have this be a proportion of tau
    cumulative.time <- (cumulative.time*tau + stage.output[, "event.time"]) / tau # updating for the next stage
    
    ## update the prior visit length column in the output array for the previous stage
    prior.visit.length    <- stage.output[, "event.time"]
    

      # if it's not the initial stage, input the value from the updated state (using the formulas)
      # this will use values from the previous state
      # NOTE: this lags behind by one iteration, that's why we can update the prior.visit.length like we do above 
    
        ## we have first D1 term with the I() be opposite effects to find optimal rule immediately 
      input.state = rho*stage.output[, "state"] + ifelse(stage.output[, "state"] > 0.5, 1, 0)*(D0 - stage.output[, "action"]*D1) +
        ifelse(stage.output[, "state"] < 0.5, 1, 0)*(D0^2 + stage.output[, "action"]*(D1)^2) +
        0.5*g*(1-(rho)^2)^0.5 * replicate(n.sample, runif(n = 1, min = 0, max = 1) - 
                                            
                                            ## subtracting a constant value to get an average 
                                            0.5)
      
      
    ## if the cumulative.time >1, we reset it to 0 and use that to plug in to the next stage
      ## we also need to subtract the difference from the previous stage's event.time
      
      #if (any(!is.na(cumulative.time) & cumulative.time > 1)) {
        
        
      #  cumulative.time[cumulative.time > 1] <- 1
        
        
      #}
    
      

  }
  
  
  return(output)
}


set.seed(123)
pts <- simulate_patients(n.sample = 30, max_stages = 20, tau = 1000,
                  ## initially all patients are at risk
                  at.risk = 1,
                  ## Life so far lived at the beginning of the stage
                  cumulative.time = 0,
                  prior.visit.length = 0, 
                  p = 1, ## "covariate" value for later state generation,
                  ## 1 for failure rate: smaller values mean larger visit times so we want this for "longer survival"
                  
                  #### Failure rate parameters get really small fairly quickly
                  ## larger magnitude values of p1 will cause the rate to shift too dramatically (holding all others constant)
                  ## if we want "normal" visit lengths which can be large values, the magnitude of h (prior visit length) must be small to lessen perturbation
                  ## changing the intercept of failure.time to be larger makes the most difference in increasing overall number of ppl making it to later stages
                  
                  ## -4.5, b1 = -1, z1 = -0.07, p1 = -0.05, g1 = 0.7, h1 = -0.2, r1 = -0.05
                  
                  a1 = -6, b1 = -0.3, z1 = -0.025, p1 = 0.02, g1 = -0.2, h1 = 0.008, r1 = 0.01,
                  ## 2 for time to next visit
                  
                  ## -1.5, b2 = -1, z2 = -0.008, p2 = -0.01, g2 = -0.7, h2 = -0.2, r2 = -0.005
                  ## it looks like this rate has the most variation, but that's because the rate is the LARGEST and most visible on the plot
                  ## I wonder if due to this, there should be less variation?
                  
                  a2 = -3, b2 = -0.3, z2 = -0.015, p2 = 0.025, g2 = -0.2, h2 = 0.008, r2 = 0.01, 
                  a3 = -8, b3 = -0.3, z3 = -0.04, p3 = 0.02, g3 = -0.2, h3 = 0.008, r3 = 0.01,
                  rho = 0.75,
                  # action-specific effects 
                  D0 = 0,
                  D1 = 1,
                  g = 0.25) 

## NOTE: works for -4.5, a2 = -3.9

pts[,, 1:5]


######### generating visualizations

generate_plots_and_summary(pts = pts, num_patients = "30", num_stages = "20", other_text = "tau1000")


##### recovering the betas from propensity scores

a <- matrix(pts[, "action", ], ncol = 1)
s <- matrix(pts[, "state", ], ncol = 1)

as_df <- as.data.frame(cbind(a, s))
names(as_df) <- c("a", "s")

glm(a ~ s, data = as_df, family = "binomial")




