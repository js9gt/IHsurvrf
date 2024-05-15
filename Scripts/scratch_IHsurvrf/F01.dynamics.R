

library(MASS); library(ggplot2); library(dplyr); library(Rcpp)


## Function starts with w/ these @return which shows what the function is supposed to output
## This serves as documentation fo the function

#' @return event.time observed time. Earlier of the failure/treatment or censoring time
#' @return gamma 1(T <= U); T = failure time, U = treatment time, T_tilde = min(T, U)
#' @return delta 1(T_tilde <= C); C = censoring time
#' @return rho tumor size
#' @return omega wellness

#################
##   Inputs    ##
#################

## time.max: maximum time in  terms of stages
## tick: time intervals
## rho.plus: tumor size?
#    if rho is greater than the upper bound, the patient is in serious state and cannot make it to the next stage.
## omega.plus
## pred.hazard
## pred.censor
## admin.censor: remaining study length. For stage 1, = tau, but for later stages, = tau - time survived.
## at.risk: a vector of who's eligible for this stage. (patients who already died or censored will not be at risk.)
## corr: correlatino of OU process (rho, omega)
## multiplier
## full_data: full_data: present all hidden event times (failure, treatment, censoring) and trajectories
#            if (full_data == 2) present all hidden event times without trajectories.
## terminal.stage
## omega.threshold

#################
##   Function  ##
#################

## simulates dynamic process over a specified time frame related to the progression of a medical condition/disease.
## Uses stochastic processes & parameters to model events & measures by generating brownian motion process
#      with specified correlations representing stochastic fluctuations in tumor size and wellness over time
## Then calculates trajectories of rho and omega over the given time frame, hazard function based on these
#      then simulates survival probabilities and determines event times (failure & treatment time) based on these

dynamics <- function(time.max = 3, tick = 0.01, rho.plus, omega.plus, 
                     pred.hazard = 0, pred.censor = 0, admin.censor = Inf,
                     at.risk = 1, corr = -0.5, multiplier = 2, full_data = FALSE, terminal.stage = FALSE,
                     omega.threshold = 0.1) {
  # at.risk  : a vector of who's eligible for this stage. (patients who already died or censored will not be at risk.)
  # full_data: present all hidden event times (failure, treatment, censoring) and trajectories
  #            if (full_data == 2) present all hidden event times without trajectories.
  # admin.censor: remaining study length. For stage 1, = tau, but for later stages, = tau - time survived.
  # corr: correlation of OU process (rho, omega)
  # if rho is greater than the upper bound, the patient is in serious state and cannot make it to the next stage.
  
  ## outputs :
  # event.time = X = min(T, D), gamma = 1(D <= T), delta = 1(not censored)
  # rho.x / omega.x = tumor size / wellness measure at the event time.
  if (!at.risk) {
    
    # If not at risk, assign NA values to output variables (specified earlier through @return)
    
    output = c(event.time = NA, gamma = NA, delta = NA,  rho.x = NA, omega.x = NA)
    if (full_data == 1) {

      ## for full_data = 1, present all hidden event times AND trajectories (rho, omega, surv) through a list
      return(list(statistics = output, 
                  times = c(failure.time = NA, treatment.time = NA,
                            censor.time = NA),
                  rho = NA, omega = NA, surv = NA))
    } else if (full_data == 2) {
      ## for full_data = 2, present all hidden event times without trajectories through a vector
      return(c(output, failure.time = NA, treatment.time = NA,
               censor.time = NA))
    } else {
      
      # Return only output if full_data is neither 1 nor 2
      
      return(output)
    }
  }
  
  ## initializes a timeframe sequence
  timeframe = seq(0, time.max, by = tick)
  ## calculates length of sequence
  n.time = length(timeframe)
  
  ## error generation
  # brownian_motion <- mvrnorm(n.time, rep(0,2), matrix(c(1,corr, corr, 1), 2) / n.time)
  
  ## Generate correlated errors using multivariate normal distribution
  ## Brownian motion is a widely used gaussian random process
  brownian_motion <- mvrnorm(n.time, c(0.1, 0.05) * tick, matrix(c(1, corr, corr, 1), 2) * tick^2)
  
  ## Generate Brownian motion with specific mean and covariance matrix, then apply a cumulative Ornstein-Uhlenbeck process
  ## cum_ou function is separately written in RCPP on the bottom.
  #      iterates through each vector element to update each subsequent element based on the OU process
  #      based on the formula: x_k = r_k + x_{k-1} * (1 - theta * dt)
  #                 where r_k is current value and x_{k-1} is previous value in vector, theta is reverting constant
  # then multiply generated brownian motion with a constant input parameter
  
  brownian_motion <- apply(brownian_motion, 2, cum_ou, reverting_const = 2 * tick) * multiplier
  
  # ggplot(data.frame(brownian_motion, col = 1:n.time), aes(X1, X2, col=col)) + geom_path()
  
  ## trajectory of rho and omega: simulating those values across various timeframes 
  ## brownian_motion[, 1]: values in the first column of that matrix
  ## then calculates a vector based on arithmetic operations involving rho.plus scalar (input) and timeframe
  ## pmax() compares each element of the resulting vector w/ 0 & returns the maximum between 0 and each element
  rho = pmax(0, brownian_motion[, 1] + 2 * rho.plus * timeframe + rho.plus)
  
  ## omega.plus is an input scalar value
  ## element-wise multiplication between a scalar and a vector, results in a vector
  ## brownian_motion[, 1]: values in the 2nd column of that matrix
  ## then compares each element of the resulting vector with 0.1, and returns the maximum between 0.1 and vector elements
  omega = pmax(0.1, omega.plus + (1 - omega.plus) * (1 - exp(- timeframe/2)) + brownian_motion[, 2])
  
  
# tmp.bm <<- brownian_motion
  ## event times
  # surv = exp(-cumsum(timeframe * tick *4 /3 * exp(-brownian_motion[, 2]/3 + pred.hazard) * rho.plus))
  
  ## calculating hazard and survival probability
  ## rho/omega: maybe relative risk of impact of tumor size on the hazard rate, then multiplied by expon. term
  ## adding either 0 or 1 based on whether omega <= the threshold (input value)
  ## terms scaled by tick/5
  ## formula for hazard for FAILURE TIME
  
  
  hazard = (rho / omega * exp(pred.hazard) + (omega <= omega.threshold)) * tick / 5  
  surv = exp(-cumsum(hazard))
# print(surv)
# tmp.haz <<- hazard
  
  ## determine failure time based on survival probability and randomness
  ## surv: vector of survival probabilities
  ## runinf(1) generates random number from uniform(0, 1) distribution.
  ## we find the first index in the surv vector w/ survival probability <= randomly generated value
  ## then scale by the time increment 
  failure.time = which(surv <= runif(1))[1] * tick # NA if administratively censored
  
  ## if failure time missing (doesn't occur within observed timeframe or admin censoring), set to infinity
  ## AKA, the failure event did not occur within the observed period
    if (is.na(failure.time)) failure.time <- Inf
  
  ## determine treatment time based on conditions
  # critical <- which(omega <= omega.threshold)[1]
  # if (!is.na(critical)) if (critical * tick < failure.time) failure.time = critical * tick
  
  ## terminal.stage (input) == T; simulation has reached final phase, no further treatment
  #       then trt.time is set to infinity to denote no additional treatment time
  if (terminal.stage) {
    trt.time = Inf
  } else {
    ## terminal.stage = F, then find index of first occurrence where rho <= 1, scale by time increment
    #       to determine time point at which tumor size reaches a specific threshold as trigger for trtmt
    trt.time = which(rho >= 1)[1] * tick           # NA if administratively censored
    
    ## if trt.time is missing from adminitrative censoring, set this to infinity to show info regarding trtmt
    #      is unavailable or administratively censored or condition for treatment not met
    if (is.na(trt.time)) trt.time <- Inf
  }
  
  ## summary statistics for event times
  ## min between failure time and treatment time; earliest event occurrence between failure and trtmt
  X = min(failure.time, trt.time)
  
  ## if T: failure event occurs before or simultaneously w/ treatment event
  gamma = failure.time <= trt.time
  
  ## finds index in timeframe vector that the earlier event (X) occurs
  x.index = which(timeframe == X)
  
  ## censoring time: minimum between admin.censor (remaining study length) and a random # from exp distribution w/ rate exp(pred.censor)
  censor.time = min(admin.censor, rexp(1, exp(pred.censor)))
# tmp.cens <<- pred.censor
# tmp.cens.time <<- censor.time
  
  ## computes min between earlier event time and censoring time
  XX = min(X, censor.time)
  
  ## determined whether failure or treatment occurs at the same time as censoring 
  ## TRUE: event occurs before or simultaneously as censoring AKA not censored 
  delta = (X <= censor.time)
  
  ## checks if delta is 0 (censored).  
  if (!delta) { #if censored (when delta = 0)
    
    ## if event is censored, we output information related to the censored event
    ## here, delta would be 0 which indicates censoring
    output = c(event.time = XX, gamma = NA, delta = delta, rho.x = NA, omega.x = NA)
    
    ## if full_data (input) == 1, return this list of summary statistics 
    if (full_data == 1) {
      return(list(statistics = output, 
                  times = c(failure.time = failure.time, treatment.time = trt.time,
                            censor.time = censor.time),
                  rho = rho, omega = omega, surv = surv))
      
    ## otherwise if full_data == 2, return vector of summary statistics without rho, omega, survival
    } else if (full_data == 2) {
      return(c(output, failure.time = failure.time, treatment.time = trt.time,
                censor.time = censor.time))
    } else {
    ## if neither condition is met, return only output summary statistics 
      return(output)
    }
  }
  
  
  # cat(failure.time, " ", trt.time, " ",  gamma,"\n")  
  
  ## checks length of index where earliest event occurs (x.index). If length is 0 (no event),
  if (!length(x.index)) {
    ## then set these to NA
    rho.x = NA
    omega.x = NA
  } else {
    
    ## otherwise if length isn't 0, assign these their respective values from the rho and omega vectors at that time index
    rho.x = rho[x.index]
    omega.x = omega[x.index]
  }
  
  
  ## construct output vector of relevant variables and their values
  output = c(event.time = X, gamma = gamma, delta = delta, rho.x = rho.x, omega.x = omega.x)
  
  ## these determine the format of the output
  
  if (full_data == 1) {
    
    ## return a list containing statistics
    return(list(statistics = output, 
                times = c(failure.time = failure.time, treatment.time = trt.time,
                          censor.time = censor.time),
                rho = rho, omega = omega, surv = surv))
    
    ## return a vector containing statistics
  } else if (full_data == 2) {
    return(c(output, failure.time = failure.time, treatment.time = trt.time,
             censor.time = censor.time))
  } else {
    
    ## otherwise, return only the output vector
    return(output)
  }
}

## create a vectorized version of the dynamics function to vectorize a scalar function.
## this allows it to handle input vectors or arrays for the following arguments
## now, this new function, dynamics.vec can process multiple values for these arguments

dynamics.vec <- Vectorize(dynamics, vectorize.args = c("rho.plus", "omega.plus", "pred.hazard", 
                                                       "pred.censor", "at.risk", "admin.censor"))

## OU process with parameter theta -> replaced to an external Rcpp function
  # dx_t = -theta x_t dt + sigma dW_t
  # cum_ou <- function(vec, reverting_const) { #reverting_const = theta dt
  #   output = vec
  #   for (k in  2:length(vec)) {
  #     output[k] = output[k] + output[k-1] * (1 - reverting_const)
  #     # x_k = r_k + x_{k-1} (1-theta dt)
  #   }
  #   output
  # }
  # brownian_motion <- apply(brownian_motion, 2, cumsum)

  #cum_ou in Rcpp
  cppFunction(
    "NumericVector cum_ou(NumericVector vec, double reverting_const) {
    const int n = vec.size();
    NumericVector y(n);
    y[0] = vec[0];
    for (int i=1; i < n; ++i) {
    y[i] = y[i] + y[i-1] * (1 - reverting_const);
    }
    return y;
    }")
