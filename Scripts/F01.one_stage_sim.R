
# use this to generate propensity score
source("C21.sim_params.R")

## Function starts with w/ these @return which shows what the function is supposed to output
## This serves as documentation fo the function

#' @return event.time, "visit length"; X = min(Tk, Uk), earlier of time to death or time to next visit
#' @return failure.time: Tk: failure time
#' @return gamma 1(T <= U); T = failure time, U = treatment time, which is a failure indicator. 1 means subject fails before next stage.
#' @return state; based on generated value from normal distribution
#' @return delta: delta is 1 if a patient is not censored, 0 if a patient is censored AKA their cumulative time is greater than tau


#' @return rate_failure: to track rate (mean) which generates the exponential distribution failure time
#' @return rate_time_next_visit: to track rate (mean) which generates the exponential distribution


## no longer:
# #' @return treatment.time: Uk: time to next visit (treatment time)


#################
##   Inputs    ##
#################

# nstage = initialize at 0
# cumulative_length = initialize at 0
# terminal.stage
# an at.risk marker (not at risk = NA values)
# prior_vis_length = initialize at 0

# Parameters for the failure time distribution (Tk)
#a1 <- -0.3
#b1 <- 0.1
#z1<- -0.3
#p1 <- -1
#g1 <- -0.2
#r1 <- -0.8

# Parameters for the time to next visit distribution
#a2 <- 1.2
#b2 <- -0.05
#z2 <- -2.5
#p2 <- 0.1
#g2 <- -2
#r2 <- -1

one_stage <- function(
    ## nstages is number of previous stages
  nstages = 0,
  cumulative_length = 0,
  ## length of of the prior visit (not cumulative)
  prior.visit.length = 0,
  at.risk = 1,
  #action = 0, #only used when we first generate prop scores for all stages at once, but this doesn't make sense
  terminal.stage = FALSE, 
  a1 = -0.3, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
  a2 = 1.2, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
  tau = 50,
  #dimensions of covariates for state vector generated
  p = 1) {
  
  while (TRUE) {
    if (!at.risk) {
      # If not at risk, assign NA values to output variables
      output = c(event.time = NA, gamma = NA, delta = NA,
                 failure.time = NA, treatment.time = NA, action = NA, state = NA, prior.visit.length = prior.visit.length, rate.failure = NA,
                 rate.next.visit = NA)
      return(output)
    }
    
    valid_failure_time = FALSE
    while (!valid_failure_time) {
      
      # Generate state data from a uniform(0, 1) distribution to control values
      state = runif(1, 0, 1)
      
      ## generates predicted probability of treatment for each patient as a vector, inputting generated state
      propensity_scores <- propensity_function(state = state, p)
      
      action <- rbinom(n = 1, size = 1, prob = propensity_scores)
      
      
      # Rate of failure time (Tk)
      rate_failure = exp(
        a1 + b1 * state + z1 * nstages + p1 * cumulative_length + g1 * action + h1*prior.visit.length + 
          r1 * action * state * nstages * cumulative_length*prior.visit.length
      )
      
      # Ensure rate_failure is positive
      #if (rate_failure <= 0) next
      
      # Rate of time to next visit (Uk)
      rate_time_next_visit = exp(
        a2 + b2 * state + z2 * nstages + p2 * cumulative_length + g2 * action + h2*prior.visit.length +
          r2 * action * state * nstages * cumulative_length*prior.visit.length
      )
      
      # Ensure rate_time_next_visit is positive
      #if (rate_time_next_visit <= 0) next
      
      # Generating failure time (Tk) using the rate from an exponential distribution
      failure.time = rexp(n = 1, rate = rate_failure)
      
      # Generating time to next visit (treatment.time) (Uk)using the rate from an exponential distribution
      treatment.time = rexp(n = 1, rate = rate_time_next_visit)
      
      ## delta is 1 if a patient is not censored, 0 if a patient is censored AKA their cumulative time is greater than tau
      delta = (cumulative_length <= tau)
      
      
      # Check if failure.time and treatment.time are NA or non-positive
      #if (is.na(failure.time) || failure.time <= 0 || is.na(treatment.time) || treatment.time <= 0) next
      
      ### this while loop with the conditions of bounds is what causes the process to take a long time
      if (!is.na(failure.time) && !is.na(treatment.time) &&  failure.time < 100 && treatment.time < 100) {
        valid_failure_time = TRUE
      }
    }
    
    # Terminal stage (input) == T; simulation has reached the final phase, no further treatment
    if (terminal.stage) {
      trt.time = Inf
      
    } else {
      # Terminal stage = F, find the index of the first occurrence where rho <= 1
      trt.time = treatment.time
    }
    
    
    # Summary statistics for event times
    X = min(failure.time, trt.time)
    gamma = failure.time <= trt.time
    
    # if the patient is censored (delta == 0), then their gamma is UNKNOWN, so it should be NA
    if (!delta) {
      gamma = NA
      
    } else {
      # otherwise, make no change to gamma
      gamma = gamma
    }
    
    # Construct output vector of relevant variables and their values
    output = c(event.time = X, gamma = gamma, delta = delta,
               failure.time = failure.time, treatment.time = trt.time, action = action, state = state, prior.visit.length = prior.visit.length, rate.failure = rate_failure,
               rate.next.visit = rate_time_next_visit)
    
    # Return a list containing statistics
    return(output)
  }
}

## vectorizing the following arguments
one_stage.vec <- Vectorize(one_stage, vectorize.args = c("nstages", "cumulative_length", "at.risk", "prior.visit.length"))

### example code after fitting propensity model in F02.multistage_sim.R
#dynamics.vec(nstages = rep(0, 5),
#         cumulative_length = rep(0, 5),
#         at.risk = rep(1, 5),
#         propensity_model = propensity_model,
#         terminal.stage = FALSE, 
#         a1 = -0.3, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, r1 = -0.8,
#         a2 = 1.2, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, r2 = -1)
