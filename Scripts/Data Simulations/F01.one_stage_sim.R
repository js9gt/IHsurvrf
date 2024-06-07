
# use this to generate propensity score
#source("C21.sim_params.R")
### UPDATE: 20Mar2024, propensity score is now generated as a global variable in the simulation_run.R setting

## Function starts with w/ these @return which shows what the function is supposed to output
## This serves as documentation fo the function

#' @return event.time, "visit length"; X = min(Tk, Uk), earlier of time to death or time to next visit
#' @return gamma 1(T <= U); T = failure time, U = treatment time, which is a failure indicator. 1 means subject fails before next stage.
#' @return delta: delta is 1 if a patient is not censored, 0 if a patient is censored AKA their cumulative time is <= tau
#                        this can also happen if a patient's generated censoring time happens first
#' @return failure.time: Tk: failure time
#' @return treatment.time: Uk: treatment time (if this is min pt moves to next stage)
#' @return censor.time: Uk: censoring time (if this is min pt is censored at the current stage)
#' @return time.max: remaining life AKA that the patient has-- this is needed to bound their time up to tau
                ## if a patient's generated time exceeds time.max, it is truncated at time.max
#' @return action: binary-- if a patient fails at a stage, they should not receive a treatment?
#' @return state; based on generated value from a uniform distribution and then updated based on dependency
#' @return prior.visit.length: length of the previous stage (in terms of days)
#' @return nstages: number of previous stages (as a proportion of the maximum number of stages)



#' @return rate.failure: to track rate (mean) which generates the exponential distribution failure time
#' @return rate.next.visit: to track rate (mean) which generates the exponential distribution-- time to next visit
#' @return rate.censoring: to track rate (mean) which generates the exponential distribution-- time to censoring

## output from the multistage_sim section:

#' @return at.risk: 1 if the patient is still at risk, 0 if not (then put all their variables to NA) bc they're ineligible
#' @return cumulative.time: a proportion of tau (so that our rates don't get too small)= sum of previous visits / tau


## no longer:
# #' @return treatment.time: Uk: time to next visit (treatment time)


#################
##   Inputs    ##
#################

# nstage = initialize at 0
# cumulative.length = initialize at 0
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
  cumulative.length = 0,
  ## length of of the prior visit (not cumulative)
  prior.visit.length = 0,
  at.risk = 1,
  time.max = 10000,
  terminal.stage = FALSE,
  input.state.value = 0,
  input.state.value2 = 0,
  a1 = -0.3, b1 = 0.1, c1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
  a2 = 1.2, b2 = -0.05, c2 = 0.1, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
  a3 = 1.2, b3 = -0.05, c3 = 0.1, z3 = -2.5, p3 = 0.1, g3 = -2, h3 = 0.6, r3 = -1,
  tau = 50,
  #dimensions of covariates for state vector generated-- 20MAR2024 not needed bc p is a global variable now
  #p = 1,
  ## a logical for if we want to include censoring (besides administrative censoring)
  ## this is used as an input into multistage sims
  censoringyesno = TRUE,


  ## should be a vector of actions to take at the current stage
  ## the actions to take are decided based on the optimal policy
  input.policy.action = NA,


  ## the optimal action to take at that stage
  input_opt = NA) {

  ## a logical for if we want to use propensity scores (FALSE), or use a policy (TRUE)
  ## if ALL of the input.policy.action values are NA, then we generate based on propensity score. Otherwise, we ue the input policy
  policyTF = ifelse(sum(!is.na(input.policy.action)) != 0, T, F)

  while (TRUE) {
    if (!at.risk) {
      # If not at risk, assign NA values to output variables
      output = c(event.time = NA, gamma = NA, delta = NA,
                 failure.time = NA, treatment.time = NA, censor.time = NA, time.max = NA,
                 action = NA, optimal.action = NA, state1 = NA, state2 = NA, prior.visit.length = prior.visit.length, nstages = NA, rate.failure = NA,
                 rate.next.visit = NA, rate.censoring = NA)
      return(output)
    }

    valid_failure_time = FALSE
    while (!valid_failure_time) {


        ## we input values for the state. For the first state, the value is generated from U(0, 1) distribution
        ## this is input here for single stage data generation, then will be updated based on depepdency in the next stage
      ## then, the updates state values from the dependency will be input here
       state1 = input.state.value
       state2 = input.state.value2


      ### if policy is NULL, we just generate propensity scores from the state
      if (policyTF == FALSE) {
      ## generates predicted probability of treatment for each patient as a vector, inputting generated state
      propensity_scores = predictedPropensity(state1 = state1, state2 = state2)

      ## allow action to be -1 or 1 instead of 0 and 1 so that there's a treatment effect
      ## suppressWarnings(rbinom(1, 1, propensity_scores) * 2 - 1)

      action = rbinom(1, 1, propensity_scores)

      optimal.action= input_opt


      ### if policy is NOT null, we use the input action
      } else{

        action = input.policy.action

        optimal.action= input_opt


      }


      # Rate of failure time (Tk)
      rate_failure = exp(
        a1 + b1 * state1 + c1*state2 + z1 * nstages + p1 * cumulative.length + g1 * action*state1 + h1*prior.visit.length +
          r1 * action * state1*state2 * nstages * cumulative.length*prior.visit.length
      )


      # Rate of time to next visit (Uk)
      rate_time_next_visit = exp(
        a2 + b2 * state1 + c2*state2 + z2 * nstages + p2 * cumulative.length + g2 * action*state1 + h2*prior.visit.length +
          r2 * action * state1*state2 * nstages * cumulative.length*prior.visit.length
      )


      if (censoringyesno) {

      # Rate of time to censoring (Ck)
      #delta is 1 if a patient is not censored, 0 if a patient is censored
      rate_censoring = exp(
        a3 + b3 * state1 + c3*state2 + z3 * nstages + p3 * cumulative.length + g3 * action*state1 + h3*prior.visit.length +
          r3 * action * state1*state2 * nstages * cumulative.length*prior.visit.length
      )

      # Generating time to censoring (censor.time) (Ck) using the rate from an exponential distribution
      censor.time = rexp(n = 1, rate = rate_censoring)
      } else {
        rate_censoring = NA
        censor.time = Inf
      }


      # Generating failure time (Tk) using the rate from an exponential distribution
      failure.time = rexp(n = 1, rate = rate_failure)

      # Generating time to next visit (treatment.time) (Uk) using the rate from an exponential distribution
      treatment.time = rexp(n = 1, rate = rate_time_next_visit)


      ## delta is 1 if a patient is not censored, 0 if a patient is censored AKA their cumulative time is greater than tau
      ## delta is also 1 if censoring time is smaller than the treatment.time and failure.time

      ## cumualative length < 1 since it's a proportion of tau
      delta = (censor.time > min(failure.time, treatment.time) & cumulative.length < 1)



      ### this while loop with the conditions of bounds is what causes the process to take a long time
      if (!is.na(failure.time) && failure.time!= 0 && !is.na(treatment.time) && treatment.time != 0 && !is.na(censor.time)
          && censor.time !=0){
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

    # Adjust X if greater than time.max
    X = ifelse(X > time.max, time.max, X)

    ## if we truncate the value of X at time.max, delta would also need to be 0
    #### this is for administrative censoring
    delta = ifelse(X == time.max, 0, delta)

    gamma = failure.time <= trt.time


    # if the patient is censored (delta == 0), then their gamma is UNKNOWN, so it should be NA
    ## we also need to change their event.time to be either be 0, or tau - x (they cannot go beyond tau)
    if (!delta) {
      gamma = NA


    } else {
      # otherwise, make no change to gamma
      gamma = gamma
    }

    # Construct output vector of relevant variables and their values
    output = c(event.time = X, gamma = gamma, delta = delta,
               failure.time = failure.time, treatment.time = trt.time, censor.time = censor.time, time.max = time.max,
               action = action, optimal.action = optimal.action, state1 = state1, state2 = state2, prior.visit.length = prior.visit.length, nstages = nstages, rate.failure = rate_failure,
               rate.next.visit = rate_time_next_visit, rate.censoring = rate_censoring)

    # Return a list containing statistics
    return(output)
  }
}

## vectorizing the following arguments
one_stage.vec <- Vectorize(one_stage, vectorize.args = c("nstages", "cumulative.length", "input.state.value", "input.state.value2", "at.risk",
                                                         "prior.visit.length", "time.max", "input.policy.action", "input_opt"))

### example code after fitting propensity model in F02.multistage_sim.R
#dynamics.vec(nstages = rep(0, 5),
#         cumulative.length = rep(0, 5),
#         at.risk = rep(1, 5),
#         propensity_model = propensity_model,
#         terminal.stage = FALSE,
#         initial.stage = FALSE,
#         input.state.value = rep(0, 5),
#         a1 = -0.3, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
#         a2 = 1.2, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
#         a3 = 1.2, b3 = -0.05, z3 = -2.5, p3 = 0.1, g3 = -2, h3 = 0.6, r3 = -1,
#         p =1, tau = 70)
