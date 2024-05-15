

source("F01.dynamics.R")


## documentation describing some input parameters of the function

#' @param p number of covariates
#' @param Sig covariate error covariance
#' @param corr cor of two error processes (rho and lambda)

## documentation describing the outputs of the function

#' @return event.time observed time
#' @return gamma 1(death <= treatment): 1 means died first, 0 means got treatment first
#' @return delta 1(min(death, treatment) <= censor)
#' @rho.x tumor size at the event
#' @omega.x wellness at the event
#' @at.risk Availability at the beginning of the stage
#' @surv.previous Life so far lived at the beginning of the stage
#' @covariates(Z1-Zp) covariate values at the beginning of the stage


## relies on dynamics.vec function from F01.dynamics which simulates survival times in dif scenarios for param values
## executes multi-stage simulation by iterating through each stage
## for each stage:
##     determines an action based on a policy or propensity scores
##     updates current stage & progresses to next state
##     calculate covariates for subsequent state

## outputs a summary of the simulations w/ info about subjects, stages, action taken, risk status, event time,
##        censoring status, terminal stage


multiStageDynamics <- 
  function(n.sample = 100, n.stages, tau = 10, tick = 0.01,          # structural parameters
           rho = NULL, omega = NULL, surv.previous = rep(0, n.sample), # initial state vector
           covariate = NULL, at.risk = rep(1, n.sample),             
           p = 5, Sig = diag(p) + 0.2 - diag(p) * 0.2,               # covariate structure
           corr = -0.5,                                               # cor of two error processes
           predHazardFn, predPropensityFn, predCensorFn,             # list of predictor functions
           hidden_data = FALSE, summary = TRUE,                      # output control
           policy = NULL,                                            # optimal rule (if !is.null, propensity scores are ignored.) for value calculation
           printFlag = TRUE                        # policy is an object of rsf.obj list.
           ) {
    
    
    ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
    if (length(at.risk) == 1) at.risk = rep(at.risk, n.sample)
    
    ## if length of surv.previous (input) doesn't match number of samples, raises an error
    if (length(surv.previous) != n.sample)  stop ("length of surv.previous should match.")
    
    ## if length of at.risk (input) doesn't match number of samples, raises an error
    if (length(at.risk) != n.sample)  stop ("length of at.risk should match.")
    
    ## if value of omega is input (not null) and length doesn't match the number of samples, raises error
    if (!is.null(omega)) if (length(omega) != n.sample)  stop ("length of omega should match.")
    
    ## if value of covariate is input (not null) and number of rows (observations) doesn't match number of samples, raises error
    if (!is.null(covariate)) if (dim(covariate)[1] != n.sample) stop ("length of covariate should match.")
    
    ## if value of rho is input (not null) and length doesn't match the number of samples, raises error
    if (!is.null(rho)) if(length(rho) != n.sample)  stop ("length of rho should match.")
    
    ## calculates maximum time allowed in simulation based on diff between tau (limit) and existing time lived
    time.max = max(tau - surv.previous)
    
    # initial state
    
    ## if rho is null (not provided externally), generate random values between 0.5 and 1
    if (is.null(rho))           rho = runif(n.sample) * 0.5 + 0.5
    
    ## if omega is null (not provided externally), generate random values between 0.5 and 1
    # if (is.null(rho))           rho = rep(1, n.sample)
    if (is.null(omega))         omega = runif(n.sample) * 0.5 + 0.5
    
    ## if surv.previous is null (not provided externally), initializes as a vector of 0's 
    ##         duration of survival up to the start of each stage
    if (is.null(surv.previous)) surv.previous = rep(0, n.sample)
    
    ## if covariate is null (not provided externally), generate MVN random values w/ mean 0 and covariance matrix sig (input)
    if (is.null(covariate))     covariate = mvrnorm(n.sample, rep(0, p), Sig)
    
    ## if there are no column names for covariate, assigns names Z1, Z2.. ZP
    if (is.null(colnames(covariate))) colnames(covariate) = paste0("Z", 1:p)
    
    # skeleton
    ## uses the dynamics.vec() function to create an initial structure to set up an initial state for 
    ## later dynamics of the multi-stage
    tmp <-
      dynamics.vec(time.max = 1, tick = tick, rho.plus = rep(0.5, 2), omega.plus = rep(0.5, 2), 
                   pred.hazard = 0, pred.censor = 0, admin.censor = Inf,
                   at.risk = 1, corr = corr, full_data = ifelse(hidden_data, 2, 0)) %>% t

    ## initializes an array called "output" filled with NA values with dimensions:
    ##        number of samples being simulated (input)
    ##        dim(tmp)[2]: number of columns in the tmp opject (AKA the output variables from dynamics.vec)
    ##        n.stages: number of stages (input) in simulation
    output <- 
      array(NA, dim = c(n.sample, 2 + dim(tmp)[2] + 6 + p, n.stages), 
            ## row and column names for the array
            dimnames = list(1:n.sample, c("subj.id", "rep.id", colnames(tmp), "action", "at.risk", 
                                          "surv.previous", "lB", "rho.0", "omega.0", colnames(covariate)), 
                            1:n.stages))
    
    ## initializes the subj.id column we just created in the output from 1 to n.sample
    output[, "subj.id", ] <- 1:n.sample
    
    ## initialize rep.id column with 0s
    output[, "rep.id", ] <- 0           # rep.id is reserved for later use (repeated random draws)
    
    ## initialize action vector with NA values, will be same length as number of samples
    action = rep(NA, n.sample)           # initialize action vector
    
# tmp.tmp <<- list()
    
    ## iterate through each stage
    for (stage in 1:n.stages) {
# print(paste0("stage ", stage))
      
      ## if printFlag (input) == TRUE, then pring the current stage number
      if (printFlag) cat("stage ", stage, "\n")
      
      ## update the at.risk column in the output array for the current stage using the value from at.risk
      output[, "at.risk", stage]           <- at.risk
      
      ## update the surv.previous column in the output array for the current stage using the surv.previous vector
      output[, "surv.previous", stage]     <- surv.previous
      
      ## calculate the log the survival + 1 and uses this to update the lB column for the current stage with these values
      output[, "lB", stage]                <- log(surv.previous + 1)
      
      ## updates corresponding matching olumn names of the covariate() matrix with the values in there
      output[, colnames(covariate), stage] <- covariate
      
      ## updates rho.0 column in output array with values from rho vector
      output[, "rho.0", stage]             <- rho
      
      ## updates omega.0 column in output array with values from omega vector
      output[, "omega.0", stage]           <- omega
# tmp.out <<- output 
      ## action                (t=0)
      
      ## checks whether a policy has been defined (input)
      if (!is.null(policy)) {
#         if (attr(policy, "class") %in% c("Goldberg-lm")) {
#           vars <- policy[[stage]]$model %>% names
#           vars <- vars[!grepl("(weights)", vars)]  # exclude (weights)
#           args <- c(list(policy[[stage]],
#                          newdata = output[as.logical(at.risk), c("subj.id", vars), stage] %>% as.data.frame,
#                          return.optimal.Q = FALSE))
#           action[at.risk != 0] <- do.call(decision.rule.gk, args)$rule.fix
#         } else if (attr(policy, "class") %in% c("Goldberg-rf")) {
        #   args <- c(list(policy[[stage]],
        #                  newdata = output[as.logical(at.risk), c("subj.id", policy[[stage]]$xvar.names), stage] %>% as.data.frame,
        #                  return.optimal.Q = FALSE))
        #   action[at.risk != 0] <- do.call(decision.rule.gk, args)$rule.fix
        # } else 
        
        
        ## if the policy class is "dwSurv"
        if (attr(policy, "class") %in% c("dwSurv")) {
    # tmp.out <<- output      
          
          ## prepares the output data in a wide format
          dat.wide <- dwTrans(output, n.stages = n.stages, p = p)
          # vars.tf   <- names(get_all_vars(policy$tf.mod[[stage]], dat.wide))
          # vars.blip <- names(get_all_vars(policy$blip.mod[[stage]], dat.wide))
          # vars <- unique(c(vars.tf, vars.blip))
          
          ## construct arguments using the wide format data for the decision.rule.dw function
          args <- c(list(policy,
                         newdata = dat.wide[as.logical(at.risk), ],
                           # output[as.logical(at.risk), c("subj.id", vars), stage] %>% as.data.frame,
                         stage = stage,
                         return.optimal.Q = FALSE))
          action[at.risk != 0] <- do.call(decision.rule.dw, args)$rule.fix
          
          ## if the policy class is "cho"
        } else if (attr(policy, "class") %in% c("cho")) {  #### old code
          
          ## prepare arguments by extracting specific columns from the output data to input into decision rule
          args <- c(list(policy[[stage]],
                         newdata = output[as.logical(at.risk), c("subj.id", policy[[stage]][[1]]$xvar.names), stage] %>% as.data.frame,
                         return.optimal.s = FALSE, tau = tau),
                    attr(policy[[stage]], "criterion"))
          
          ## uses decision.rule.list with the specified arguments to determine action for subjects at risk
          action[at.risk != 0] <- do.call(decision.rule.list, args)$rule.fix
          
          
          
          
          
          
          ## when the policy class is DTRSurv
          ######## INPUT POLICY for evaluation
        } else if (attr(policy, "class") %in% c("DTRSurv")) {
          
          ######### JANE'S CODE
          
          ## puts output data into a specific format by setting specific columns to dummy values
          ## uses output2observable() function from F00.generic.R which rearranges data into each row is subject,
          ##     columns have time, delta, action covariate info for each stage
          ##     also computes new delta column representing censoring based on values of 0 in each stage's delta column
          
          x = output[as.logical(at.risk),,,drop = FALSE] %>% output2observable()
          
          
          
          
          ## selects columns in x  by concatenating T_ and delta_ with stage:
          ### ex) T_1, delta_1
          ### then assign these to 0
          ## NOTE: we are looping across each stage
          x[, paste0(c("T_", "delta_"), stage)] = 0  # dummy values in order for the package not to drop the NA rows.
          
          ## recall that the input "policy" is a DTRSurv object from the completed RF output
          ## we want to use the FINAL forest (at the last step)-- this is in slot pool1_final
          ## now, we retrieve the model used at each stage (the current stage from the loop): Surv(T_1, delta_1) ~ Z1_1 + Z2_1 + Z3_1 + Z4_1 + Z5_1
          ## get_all_vars() retrieves all the variables from this stage's model, then updates the dataframe to only include these
          ### we want to retrieve vars_stage variables, since we have a common model across all treatment times
          
          
          # Extract the response variable from the formula
          response_var <- as.character(formula(policy@FinalForest@model)[[2]][-1])
          
          # Append the current stage number to each term
          response_with_stage <- paste0(response_var, "_", stage)
          
          # Extract the terms of the formula excluding the response variable
          terms <- attr(terms(policy@FinalForest@model), "term.labels")
          
          # Append the current stage number to each term
          terms_with_stage <- paste0(terms, "_", stage)
          
          
          # Reconstruct the formula
          updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
          
          x = get_all_vars(updated_formula, x)
          
          ### for the new data, we need to get rid of the stage labels since our model doesn't have any stage labels:
          # Remove "_stage" suffix from all column names
          new_col_names <- gsub(paste0("_", stage, "$"), "", colnames(x))
          colnames(x) <- new_col_names
          
          ## constructs argument to use in the predict() function in the DTRSurv policy
          ## new data is the observed data generated (pts have not yet been fed through the forest)
          ## policy is the optimal estimated policy (DTRSurv Object), with a stage
          ##### NOTE: we need to change the "predict" to work using the FinalForest slot
          
          args <- list(policy, newdata = x)
          
          ## for at.risk != 0 (pt still at risk), retrieve the optimal treatment after feeding this new data into the trained forest
          ## feeds this into predict() function which is defined in class_DTRSurv.R
          ## acts on objects of class DTRSurv
          ## the new data we are feeding through the forest is the covariates at the current stage
          
          
          ## this then subsets to objects of class DTRSurvStep, and calls .Predict() in class_IH.DTRSurv.R
          ## .Predict() on objects of DTRSurvStep is defined in class_IH.DTRSurvStep.R
          ## this subsets to the "SurvRF" object to act on which is defined in class_IH.SurvRF.R
          
          ## follow the predicted optimal policy based on the input policy
          action[at.risk != 0] <- do.call(predict, args)$optimal@optimalTx
          
          
          
          
        } else if (attr(policy, "class") %in% c("GKLM")) { #### new function!
          
          ## uses same formatting using output2observable()
          x = output[as.logical(at.risk),,] %>% output2observable()
          
          ## uses function predict.opt.sep to predict and assign an optimal treatment(optimal.Tx)
          action[at.risk != 0] <- 
            predict.opt.sep(policy[[stage]],  newdata = x, Tx.label = paste0("A_", stage))$optimal.Tx %>% 
            ## transformation to modify the treatment label
            {ifelse(.==1, -1, 1)}
          
          ## if the policy class is GKRF
        } else if (attr(policy, "class") %in% c("GKRF")) { #### new function!
          
          ## uses same formatting using output2observable()
          x = output[as.logical(at.risk),,] %>% output2observable()
          x = x[, c(paste0("A_", stage), policy[[stage]][[1]]$xvar.names)]
          x[, paste0("A_", stage)] = 0 # dummy values
          
          action[at.risk != 0] <- 
            predict.opt.sep(policy[[stage]],  newdata = x, Tx.label = paste0("A_", stage))$optimal.Tx %>% 
            {ifelse(.==1, -1, 1)}
        }
        
# if (stage == 2) stop("")
        
        ## if a policy is provided, actions for subjects not at risk are set to NA based on the policy
    
        action[at.risk == 0] <- NA
      } else {
        
        ## if a policy was not assigned (null input), use predPropensityFn function for the current stage
        ##     calculates propensity values for each subject based on the parameters:
        ##     surv.previous, rho, omega, covariate
        
        propensity = predPropensityFn[[stage]](surv.previous = surv.previous, rho = rho, 
                                               omega = omega, covariate = covariate)
        
        ## for action,generate random binary outcome based on propensity values, then scale this to 0
        action = rbinom(n.sample, 1, propensity) # 1 for aggressive and 0 for gentle
        # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
        
        ## if a policy isn't assigned (input), then those not at risk are still assigned 0
        action[at.risk == 0] <- NA
      }
      
        
      ## instantaneous state   (t=0+)
      
      ## updates rho.plus. If action > 0, rho.plus is set to rho/omega/10
      ##      otherwise, it's set to rho/omega/4
      rho.plus = rho / omega / ifelse (action > 0, 10, 4)
      ## updates omega.plus. If action > 0, omega.plus set to omega - 0.5
      ##      otherwise, it's set to omega - 0.25
      omega.plus = omega - ifelse(action > 0, 0.5, 0.25)
      
      ## progress              (0 < t < T_k)
      
      ## calculates predicted hazard function for current stage for each subject
      pred.hazard = predHazardFn[[stage]](surv.previous = surv.previous, rho = rho, 
                                             omega = omega, action = action, covariate = covariate)
      
      ## calculates predicted censoring probability for current stage for each subject
      #### if there is no censoring, this just returns a matrix of -10 values of dimension surv.prev
      ## basically this is used in dynamics.vec to generate a really large value of censoring (>>> tau) 
      ## aka there's no censoring except administrative censoring
      pred.censor = predCensorFn[[stage]](surv.previous = surv.previous, rho = rho, 
                                             omega = omega, action = action, covariate = covariate)
      
      ## gathers and stores calculated variables and metrics using new pred.censor
tmp3 <<- list(surv.previous = surv.previous, rho = rho, 
              omega = omega, action = action, covariate = covariate,
              pred.censor = pred.censor,
              cens.fn = predCensorFn[[stage]])
# print(92)
# tmp3 <<- list(time.max = time.max, tick = tick, rho.plus = rho.plus, omega.plus = omega.plus, 
#               pred.hazard = pred.hazard, pred.censor = pred.censor, at.risk = at.risk,
#               corr = 0.3, admin.censor = tau - surv.previous, 
#               full_data = ifelse(hidden_data, 2, 0),
#               terminal.stage = (stage == n.stages))
      stage.output <- 
        
        ## use dynamics.vec function with updated pred.hazard, pred.censor, rho.plus,omega.plus, time.max,
        ##     admin.censor, at.risk, terminal.stage to generate values for current stage
        dynamics.vec(time.max = time.max, tick = tick, rho.plus = rho.plus, omega.plus = omega.plus, 
                     pred.hazard = pred.hazard, pred.censor = pred.censor, at.risk = at.risk,
                     corr = corr, admin.censor = tau - surv.previous, 
                     full_data = ifelse(hidden_data, 2, 0),
                     terminal.stage = (stage == n.stages)) %>% t
      
# tmp3 <<- stage.output
# print(104)      
      ## bookkeeping and reassigning
      
      ## combining data from the stag.output with "action" and "at.risk" into a matrix format
      ##      combined data includes info on dynamics of simulation, actions, status of individuals
      output[, c(colnames(tmp), "action", "at.risk"), stage] <- cbind(stage.output, action, at.risk)
      
      ## whether a patient is at risk is updated based on values of gamma. gamma = 0 means at risk
      at.risk <- (stage.output[, "gamma"] == 0)  # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
      ## replaces NA values in at.risk with 0 (not at risk)
      at.risk[is.na(at.risk)] = 0
      
      ## updates values of rho and omega based on the values generated in stage.output
      rho   <- stage.output[, "rho.x"  ]
      ## updates values of omega
      omega <- stage.output[, "omega.x"]
      
      ## covariate for the next stage. (actually not needed for the terminal stage)
      
      ## existing covariate multiplied element-wise by 0.5 and a MVN number 
      covariate = covariate * 0.5 + mvrnorm(n.sample, rep(0, p), Sig) * 0.5
      ## then update the first two columns using weighted averages by the newly generated vals
      ##      and then weighting with rho and omega values as well
      covariate[, 1] <- covariate[, 1] * 0.2 + rho * 0.8
      covariate[, 2] <- covariate[, 2] * 0.2 + omega * 0.8
      
      ## update the survival time to prepare for next state calculations by adding the current stage's event time
      surv.previous <- surv.previous + stage.output[, "event.time"] # updating for the next stage
    }
    
    ## end of for loop going through each stage
    
    
    ## tau and p attributes (additional metadata) are assigned to "output" object
    attr(output, "tau") = tau
    attr(output, "p") = p
    
    ## summary calculation where if summary = T (an input into function), calculate summary statistics
    if (summary) {
      ## cumulative event time is cum sum of the event time for each subject
      cum.event.time = apply(output[, "event.time",], 1, cumsum) %>% t
      
      ## constructs a tibble ("output.summary") containing the following columns
      output.summary <- 
        tibble(subj.id = output[, "subj.id", 1],
               rep.id = output[, "rep.id", 1],
               terminal.stage = apply(output[, "at.risk", ], 1, sum), 
               cumulative.event.time = cum.event.time[cbind(1:n.sample, terminal.stage)],
               censor.status = output[cbind(1:n.sample, "delta", terminal.stage)],
               actions       = apply(output[, "action", ], 1, 
                                     function(s) gsub("NA", "*", paste0(s, collapse = ""))))
      
      ## return a list containing both the "output" data as well as the "output.summary" requested above
      return(list(output = output, summary = output.summary))
    } else {
      
      ## if summary (input) is FALSE, only return the "output" data
      return(output)
    }
}