
#source("F01.one_stage_sim.R")
source("~/survrf/Scripts/Data Simulations/F01.one_stage_sim.R")
#source("F03.Figs.R")
source("~/survrf/Scripts/Data Simulations/F03.Figs.R")
source("~/survrf/Scripts/Data Simulations/F00.generic.R")
#~/survrf/Scripts/F01.one_stage_sim.R
#~/survrf/Scripts/F03.Figs.R



# Function to simulate data for 10 patients with specified conditions
## test
# n.sample = 100
# max_stages = 50
# tau = 100
# p = 1




simulate_patients <- function(..., n.sample, max_stages, tau,
                              ## initially all patients are at risk
                              at.risk = rep(1, n.sample),
                              ## Life so far lived at the beginning of the stage
                              cumulative.time = rep(0, n.sample),

                              ## inital length of previous visits
                              prior.visit.length = rep(0, n.sample),
                              ## dimensions of the state vector generated at each stage-- 20MAR2024 this is a global variable now
                              #p = 1,
                              a1 = -2, b1 = 0.1, z1 = -0.3, p1 = -1, g1 = -0.2, h1 = 0.2, r1 = -0.8,
                              a2 = -1, b2 = -0.05, z2 = -2.5, p2 = 0.1, g2 = -2, h2 = 0.6, r2 = -1,
                              a3 = -4.5, b3 = -1, z3 = -2.5, p3 = -0.1, g3 = 0.2, h3 = -0.6, r3 = -1,
                              ## "covariate" value for later state generation
                              rho = 0.5,
                              # action-specific effects
                              D0 = 0,
                              D1 = 1,
                              g = 1,

                              ## a logical for if we want to include censoring (besides administrative censoring)
                              censoringyesno = TRUE,

                              ## for inputting dif policies used to generate data
                              ## this is a DTRSurv object
                              policy = NULL,
                              summary = TRUE,
                              ## stage start-- if this is input, then we only calculate up to that for the summary function
                              stage.start = NULL) {

  ## if the stage.start is input,we generate that number of stages instead of max_stages (AKA change value of max_stages)

  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(at.risk) == 1) at.risk = rep(at.risk, n.sample)

  ## if length of at.risk (input) is 1, then replicate the value until it matches the length of n.sample
  if (length(cumulative.time) == 1) cumulative.time = rep(cumulative.time, n.sample)

  ## if prior visit length is 0, then replicate the value until it matches the length of n.sample
  if (length(prior.visit.length) == 1) prior.visit.length = rep(prior.visit.length, n.sample)

  #### create TF variable for usage of policy to be plugged into single stage sims
  ## this tells us whether or not to use propensity score to generate action vs one that's input
  ## input meaning the optimal policy from a method:
  ## if a policy is null, then policyTF == FALSE
  ## if a policy is input, then policyTF == TRUE
  #policyTF <- !is.null(policy)

  ## 28MAY2024: this has been modified to populate based on if the input actions are all NA (then generate w prop score) or not (then use input policy)

  ## if SS is NULL, then we use the maximum number of stages (we don't consider going back any stages)
  stage.start <- ifelse(is.null(stage.start), max_stages, stage.start)



  require(dplyr)



  # skeleton
  ## uses the dynamics.vec() function to create an initial structure to set up an initial state for
  ## later dynamics of the multi-stage
  tmp <- one_stage.vec(nstages = 0, cumulative.length = 0, prior.visit.length = 0, at.risk = 1, time.max = 1000,
                  ## using any action for now just to get column names for the structure
                  terminal.stage = F, input.state.value = 0, input.state.value2 = 0,
                  a1 = -4.5, b1 = -1, z1 = -0.07, p1 = -0.05, g1 = 0.7, h1 = -0.2, r1 = -0.05,
                  ## 2 for time to next visit
                  a2 = -3.9, b2 = -1, z2 = -0.008, p2 = -0.01, g2 = -0.7, h2 = -0.2, r2 = -0.005,
                  a3 = -4.5, b3 = -2, z3 = -0.1, p3 = -0.1, g3 = 0.2, h3 = -0.4, r3 = -0.05,
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

  true.opt <- rep(NA, n.sample)
  IHsurvrf_optimal <- rep(NA, n.sample)




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

    ### to track true optimal and predicted optimal at each stage
    #output[, "true.opt.A", stage] <- true.opt
    output[, "IHsurvrf.A", stage] <- IHsurvrf_optimal
    

    ## we run this chunk if there was a policy input
    #########################

    if (!is.null(policy)) {

      ### if policy == DTRSurv object, then we get the optimal policy for ppl at risk in this stage
      ## the extracted actions are then input into the one_stage.vec

      if (attr(policy, "class") %in% c("DTRSurvRes")) {

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
      
      
      ### NOTE: depending on the strata, we should use the forest from that strata.
      ## for the very first stage, since the cumulative time is 0, we should use strata 2 forest
      
      if (stage == 1) {
        
        forest <- policy@Forest2
        
      } else {
        
        ## Calculate the cumulative time
        ## since cumulative time is a proportion with tau as the denom, we need to multiply by tau
        cumulative_time <- output[, "cumulative.time", stage] * tau
        #print(cumulative_time)
        #print(policy@cutoff)
        
        ## Assign the appropriate forest based on the cumulative time
        if (!is.na(cumulative_time)) {
          if (cumulative_time < policy@cutoff) {
            forest <- policy@Forest2
            
            message("using strata 2 forest")
            
          } else {
            forest <- policy@Forest1
            
            message("using strata 1 forest")
          }
        }
        
      }


      # Extract the response variable from the formula
      response_var <- as.character(formula(forest@FinalForest@model)[[2]][-1])

      # Append the current stage number to each term
      response_with_stage <- paste0(response_var, "_", stage)

      # Extract the terms of the formula excluding the response variable
      terms <- attr(terms(forest@FinalForest@model), "term.labels")

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

      args <- list(forest, newdata = x)

      ## for at.risk != 0 (pt still at risk), retrieve the optimal treatment after feeding this new data into the trained forest
      ## feeds this into predict() function which is defined in class_DTRSurv.R
      ## acts on objects of class DTRSurv
      ## the new data we are feeding through the forest is the covariates at the current stage


      ## this then subsets to objects of class DTRSurvStep, and calls .Predict() in class_IH.DTRSurv.R
      ## .Predict() on objects of DTRSurvStep is defined in class_IH.DTRSurvStep.R
      ## this subsets to the "SurvRF" object to act on which is defined in class_IH.SurvRF.R

      ## follow the predicted optimal policy based on the input policy
      input.policy.action[at.risk != 0] <- do.call(predict, args)$optimal@optimalTx

      ### ADDED EXTRA
      ## for people who aren't at risk anymore, this is set to NA
      input.policy.action[at.risk == 0] <- NA

      IHsurvrf_optimal <- input.policy.action

      ##### now, we track the values for the average change in the last stage integral: this is tracked in the policy itself

      prop_last_iter_change <- tail(policy@avgKM_diff, 1)

      conv_iterations <- policy@n_it

      } else if (attr(policy, "class") %in% c("trt1")) {

        ## we give patient action == 1 at every stage

        input.policy.action[at.risk != 0] <- 1

        prop_last_iter_change <- NA

        conv_iterations <- NA


      }else if (attr(policy, "class") %in% c("trt0")){

        ## we give patient action == 0 at every stage
        input.policy.action[at.risk != 0] <- 0

        prop_last_iter_change <- NA

        conv_iterations <- NA

      }else if (attr(policy, "class") %in% c("trueopt")){

        ### if we want to know what the true optimal is and generate observational data following that,


        ### if we want to know what the true optimal is and generate observational data following that,
        stage.output1 <-

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
                        input.policy.action = 1) %>% t


        stage.output2 <-

          ## generate values for current stage
          ## for nstages, to avoid the value being overly large, have it also be a proportion of the total number of stages
          one_stage.vec(nstages = (stage - 1)/max_stages, cumulative.length = cumulative.time, at.risk = at.risk,time.max = time.max,
                        prior.visit.length = prior.visit.length,
                        #action = as.numeric(output[, "action", stage]),  ### remove this after clarification
                        terminal.stage = (stage == max_stages), input.state.value = input.state1, input.state.value2 = input.state2,
                        a1 = a1, b1 = b1, c1 = c1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
                        a2 = a2, b2 = b2, c2 = c2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
                        a3 = a3, b3 = b3, c3 = c3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau, censoringyesno = censoringyesno,

                        ## inputting action == 0
                        input.policy.action = 0) %>% t

        ## find the maximum of the event times-- for whatever the maximum is, we should use that action
        true.opt <- ifelse(stage.output1[, "event.time"] > stage.output2[, "event.time"], stage.output1[, "action"], stage.output2[, "action"])
        
        #### we want to print out what the survival times are, and take the average 

        ## use the true policy as the input policy for people at risk
        ### we omit the NA values since that's what's included for patients who are ineligible
        input.policy.action[at.risk != 0] = na.omit(true.opt)

        prop_last_iter_change <- NA

        conv_iterations <- NA

      }
      ## for people who aren't at risk anymore, this is set to NA
      input.policy.action[at.risk == 0] <- NA
    }else{

      ## if there's no policy input,
      # otherwise, it remains as NA
      input.policy.action <- input.policy.action

      prop_last_iter_change <- NA

      conv_iterations <- NA

      IHsurvrf_optimal <- NA



      }


  #  # Call the function and store the returned value
  #  true.opt2 <- compute_true_opt(policy, stage, max_stages, cumulative.time, at.risk, time.max, prior.visit.length, input.state1, input.state2,
  #                               a1, b1, c1, z1, p1, g1, h1, r1, a2, b2, c2, z2, p2, g2, h2, r2, a3, b3, c3, z3, p3, g3, h3, r3, tau, censoringyesno)
#

#    # Assume true.opt and true.opt2 have been computed previously
#    if (is.null(policy) || !attr(policy, "class") %in% c("trueopt")) {
#      selected.opt <- true.opt2
#    } else {
#      selected.opt <- true.opt
#    }

    ### if we want to know what the true optimal is and generate observational data following that,
    stage.output1 <-
      
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
                    input.policy.action = 1, input_opt = NA) %>% t
    
    
    stage.output2 <-
      
      ## generate values for current stage
      ## for nstages, to avoid the value being overly large, have it also be a proportion of the total number of stages
      one_stage.vec(nstages = (stage - 1)/max_stages, cumulative.length = cumulative.time, at.risk = at.risk,time.max = time.max,
                    prior.visit.length = prior.visit.length,
                    #action = as.numeric(output[, "action", stage]),  ### remove this after clarification
                    terminal.stage = (stage == max_stages), input.state.value = input.state1, input.state.value2 = input.state2,
                    a1 = a1, b1 = b1, c1 = c1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
                    a2 = a2, b2 = b2, c2 = c2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
                    a3 = a3, b3 = b3, c3 = c3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau, censoringyesno = censoringyesno,
                    
                    ## inputting action == 0
                    input.policy.action = 0, input_opt = NA) %>% t
    
    
#    ## avg difference in survival time between the two treatments
#    # Identify the indices where 'event.time' is NA in stage.output1
#    na_indices1 <- which(is.na(stage.output1[, "event.time"]))
#    
#    # Identify the indices where 'event.time' is NA in stage.output2
#    na_indices2 <- which(is.na(stage.output2[, "event.time"]))
#    
#    # Combine the indices and get unique ones, or create an empty vector if there are no NAs
#    all_na_indices <- if(length(na_indices1) == 0 & length(na_indices2) == 0) {
#      integer(0)
#    } else {
#      unique(c(na_indices1, na_indices2))
#    }
#    
#    # Remove the rows with these indices from both dataframes, if there are any indices to remove
#    if(length(all_na_indices) > 0) {
#      stage.output1_clean <- stage.output1[, "event.time"][-all_na_indices]
#      stage.output2_clean <- stage.output2[, "event.time"][-all_na_indices]
#    } else {
#      stage.output1_clean <- stage.output1[, "event.time"]
#      stage.output2_clean <- stage.output2[, "event.time"]
#    }
    
    
    # Now selected.opt is either true.opt or true.opt2 based on the condition
    # You can use selected.opt in the subsequent code
    stage.output <- one_stage.vec(
      nstages = (stage - 1) / max_stages, cumulative.length = cumulative.time, at.risk = at.risk, time.max = time.max,
      prior.visit.length = prior.visit.length,
      terminal.stage = (stage == max_stages), input.state.value = input.state1, input.state.value2 = input.state2,
      a1 = a1, b1 = b1, c1 = c1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
      a2 = a2, b2 = b2, c2 = c2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
      a3 = a3, b3 = b3, c3 = c3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau, censoringyesno = censoringyesno,
      input.policy.action = input.policy.action, input_opt = NA
    ) %>% t




    ## combining data from the stage.output with "action" and "at.risk" into a matrix format
    ##      combined data includes info on dynamics of simulation, actions, status of individuals
    output[, c(colnames(tmp),"at.risk"), stage] <- cbind(stage.output, at.risk)

    ## whether a patient is at risk is updated based on values of gamma and delta. gamma = 0 means at risk
    # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
    ## once delta is 0, the patient is no longer at risk also, only those who are non-censored (delta == 1) are at are risk
    at.risk <- (stage.output[, "gamma"] == 0 & stage.output[, "delta"] == 1)
    
    
    #### calculating the difference in treatment
    ####
    ####
    
    # Compute the absolute difference
    absolute_difference <- abs(stage.output1[, "event.time"] - stage.output2[, "event.time"])
    
    # Replace with NA where stage.output[, "delta"] == 0
    result <- ifelse(stage.output[, "delta"] == 0, NA, absolute_difference)
    
    ## output the treatment differences; but if delta == 0 (pt was censored, we don't calculate a treatment diff and put NA instead)
    output[, "trt.diff", stage] <- result
    
    ####
    ####
    
    ### if there are any treatment differences that a
    

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

    ## update true optimal:
    IHsurvrf_optimal <- IHsurvrf_optimal




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
    cum.event.time = apply(output[, "event.time",1:stage.start], 1, cumsum) %>% t
    output.summary <-
      tibble(subj.id = output[, "subj.id", 1],
             rep.id = output[, "rep.id", 1],
             terminal.stage = apply(output[, "at.risk", 1:stage.start], 1, sum),
             cumulative.event.time = cum.event.time[cbind(1:n.sample, terminal.stage)],
             censor.status = output[cbind(1:n.sample, "delta", terminal.stage)],
             actions       = apply(output[, "action", ], 1,
                                   function(s) gsub("NA", "*", paste0(s, collapse = ""))),
             true.opt  = apply(output[, "optimal.action", ], 1,
                               function(s) gsub("NA", "*", paste0(s, collapse = ""))),
             IHsurvrf_actions = apply(output[, "IHsurvrf.A", ], 1,
                                      function(s) gsub("NA", "*", paste0(s, collapse = ""))),
             avg.last.iter.change = prop_last_iter_change,
             convergence.iterations = conv_iterations)
    return(list(output = output, summary = output.summary))
  } else {
    return(output)
  }
}

## computing true opt is only each stage's maximum so it's not actually useful
#compute_true_opt <- function(policy, stage, max_stages, cumulative.time, at.risk, time.max, prior.visit.length, input.state1, input.state2,
#                             a1, b1, c1, z1, p1, g1, h1, r1, a2, b2, c2, z2, p2, g2, h2, r2, a3, b3, c3, z3, p3, g3, h3, r3, tau, censoringyesno) {
#  true.opt <- NULL  # Initialize true.opt
#
#  if (is.null(policy) || !attr(policy, "class") %in% c("trueopt")) {
#    # Generate the values of stage.output1 and stage.output2
#    stage.output1 <- one_stage.vec(
#      nstages = (stage - 1) / max_stages, cumulative.length = cumulative.time, at.risk = at.risk, time.max = time.max,
#      prior.visit.length = prior.visit.length,
#      terminal.stage = (stage == max_stages), input.state.value = input.state1, input.state.value2 = input.state2,
#      a1 = a1, b1 = b1, c1 = c1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
#      a2 = a2, b2 = b2, c2 = c2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
#      a3 = a3, b3 = b3, c3 = c3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau, censoringyesno = censoringyesno,
#      input.policy.action = 1
#    ) %>% t
#
#    stage.output2 <- one_stage.vec(
#      nstages = (stage - 1) / max_stages, cumulative.length = cumulative.time, at.risk = at.risk, time.max = time.max,
#      prior.visit.length = prior.visit.length,
#      terminal.stage = (stage == max_stages), input.state.value = input.state1, input.state.value2 = input.state2,
#      a1 = a1, b1 = b1, c1 = c1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
#      a2 = a2, b2 = b2, c2 = c2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
#      a3 = a3, b3 = b3, c3 = c3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3, tau = tau, censoringyesno = censoringyesno,
#      input.policy.action = 0
#    ) %>% t
#
#    # Compute true.opt based on stage.output1 and stage.output2
#    true.opt <- ifelse(stage.output1[, "event.time"] > stage.output2[, "event.time"], stage.output1[, "action"], stage.output2[, "action"])
#  }
#
#  return(true.opt)  # Return true.opt
#}


#set.seed(123)
#pts <- simulate_patients(n.sample = 50, max_stages = 25, tau = 1000,
#                  ## initially all patients are at risk
#                  at.risk = 1,
#                  ## Life so far lived at the beginning of the stage
#                  cumulative.time = 0,
#                  prior.visit.length = 0,
#                  ## 1 for failure rate: smaller values mean larger visit times so we want this for "longer survival"
#
#                  #### Failure rate parameters get really small fairly quickly
#                  ## larger magnitude values of p1 will cause the rate to shift too dramatically (holding all others constant)
#                  ## if we want "normal" visit lengths which can be large values, the magnitude of h (prior visit length) must be small to lessen perturbation
#                  ## changing the intercept of failure.time to be larger makes the most difference in increasing overall number of ppl making it to later stages
#
#                  ## -4.5, b1 = -1, c1 = -0.5, z1 = -0.07, p1 = -0.05, g1 = 0.7, h1 = -0.2, r1 = -0.05
#
#                  a1 = -6, b1 = -0.3, c1 = -0.5, z1 = -0.025, p1 = 0.02, g1 = -0.2, h1 = 0.008, r1 = 0.01,
#                  ## 2 for time to next visit
#
#                  ## -1.5, b2 = -1, z2 = -0.008, p2 = -0.01, g2 = -0.7, h2 = -0.2, r2 = -0.005
#                  ## it looks like this rate has the most variation, but that's because the rate is the LARGEST and most visible on the plot
#                  ## I wonder if due to this, there should be less variation?
#
#                  a2 = -3, b2 = -0.3, c2 = -0.05, z2 = -0.015, p2 = 0.025, g2 = -0.2, h2 = 0.008, r2 = 0.01,
#                  a3 = -8, b3 = -0.3, c3 = -0.6,z3 = -0.04, p3 = 0.02, g3 = -0.2, h3 = 0.008, r3 = 0.01,
#                  rho = 0.75,
#                  # action-specific effects
#                  D0 = 0,
#                  D1 = 1,
#                  g = 0.25,
#                  ## FALSE means we only have administrative censoring
#                  censoringyesno = FALSE,
#                  policy = NULL,
#                  summary = TRUE)
#



### showMethods(".Predict")


## NOTE: works for -4.5, a2 = -3.9



######### generating visualizations

#generate_plots_and_summary(pts = IHsurvrf.data.rep$output, num_patients = "300", num_stages = "25", other_text = "tau1000")


##### recovering the betas from propensity scores

#a <- matrix(pts[, "action", ], ncol = 1)
#s <- matrix(pts[, "state", ], ncol = 1)
#
#as_df <- as.data.frame(cbind(a, s))
#names(as_df) <- c("a", "s")
#
#glm(a ~ s, data = as_df, family = "binomial")
#



