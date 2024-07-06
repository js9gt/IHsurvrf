

IHdtrConv_otherstrata <- function(data,
                      
                      ## a numeric to indicate the strata we are assessing for convergence
                      strata,
                      
                      ## the input forest results from the previous iteration
                      prev.iteration = NULL, nDP, params, nTimes, models, mTry, ..., long_data,
                      prev_probs){
  
  
  ######
  ###### NOTE: also input long_data & output long data
  
  ###########
  ########### have a stage-wise eligibility matrix so that each patient has same # stages
  ########### for patients who are not eligible, their matrices should just be 0
  
  
  ## using the long data that was input
  data <- data
  
  prev.iteration <- prev.iteration
  
  nDP <- nDP
  
  params <- params
  
  nTimes <- nTimes
  
  models <- models
  
  mTry <- mTry
  
  long_data <- long_data
  
  last_strata_prob <- prev_probs
  
  
  
  ########## for the last stage (stage 10) we don't need to append so we do this outside the loop
  ## start at stage nDP:
  message("Convergence Strata ", strata, " starting prediction from stage ", nDP)
  
  ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
  ### if the current iteration is nDP - 2, we want to predict for nDP - 1
  ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24
  
  ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
  response_var <- as.character(formula(prev.iteration@FinalForest@model)[[2]][-1])
  response_with_stage <- paste0(response_var, "_", (nDP))
  terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")
  terms_with_stage <- paste0(terms, "_", (nDP))
  updated_formula <- paste("Surv(", paste(response_with_stage, collapse = " , "), ") ~ ", paste(terms_with_stage, collapse = " + "))
  x = get_all_vars(updated_formula, data %>% filter(!!sym(paste0("strata", strata,"_",nDP)) == 1))
  ## remove stage suffixto use in prediction
  new_col_names <- gsub(paste0("_", (nDP), "$"), "", colnames(x))
  colnames(x) <- new_col_names
  last.stage.pred <- .Predict(object = prev.iteration@FinalForest,
                              newdata = x,
                              params = params,
                              findOptimal = T)
  
  
  ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
  # Iterate over each stage of interest
  
  # Extract optimal treatments for the current stage
  optimal_treatments <- last.stage.pred$optimal@optimalTx
  
  # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
  ### patients who are eligible are ones who have a complete case for x
  eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )
  
  # also update "A" column
  ### we select the eligible values which are in the current stage, and strata
  ### then we subset t only those eligible (have complete cases)
  long_data$A[long_data$stage == nDP & long_data[[paste0("strata", strata)]] == 1] [which(eligibility == 1)] <- optimal_treatments
  
  
  #### now, we want to extract the eligibility for the whole stage
  #### this means they have to BOTH have compelte cases, AND have strata1 == 1
  ### this should now be the same length as "eligibility"
  wholestage.eligibility <- as.numeric(
    ## within stage 10 variables, which are complete cases?
    ifelse(apply(get_all_vars(updated_formula, data), 1, function(x) all(!is.na(x))) &
             
             ## within stage 10, which patients are in strata 1?
             long_data %>% filter(stage == nDP)%>% pull(!!sym(paste0("strata", strata))) == 1 ,
           ## if both conditions are true, put 1; otherwise, put 0
           1, 0) )
  
  ### next, create a matrix with a column for EACH patient, and merge these into the pooled matrix
  #### for pooled matrix, should have dim of wholestage.eligibility x nDP
  ### this way, we can sequence through these probs and each patient will have same # of columns
  
  
  
  ##########################################
  ## putting optimal predictions together ##
  ##########################################
  
  # Initialize the shifted probability matrix. each patient should have one column (instead of just eligible stages)
  ## we overwrite the eligible patients
  
  shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(wholestage.eligibility))
  
  ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
  ## do this only for eligible patients within the stage who were both in strata 1 and have complete cases
  shiftedprob1[, which(wholestage.eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)
  
  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprob_done <- shiftedprob1 - rbind(shiftedprob1[-1L,], 0.0)
  
  ## sets very small values in pr to 0
  
  shiftedprob_done[abs(shiftedprob_done) < 1e-8] <- 0.0
  
  
  ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage
  ## add the shifted optimals into the pooled optimals together
  #### each patient will have the same number of stages by design, but ones that aren't included in the strata will just have 0s
  
  # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled <- matrix(0, nrow = nTimes, ncol = ncol(shiftedprob_done)*nDP)
  
  
  # Define the index to insert columns from the previous optimal
  ## this should be for the last stage
  insert_index <- seq(from = nDP, to = ncol(pr_pooled), by = nDP)
  
  # Initialize a counter for columns from append1_pr
  col_counter <- 1
  
  # Loop over each index to insert columns from append1_pr into combined_results
  for (w in seq_along(insert_index)) {
    # Insert columns from append1_pr into combined_results
    pr_pooled[,insert_index[w]] <-
      shiftedprob_done[, w]
    
    # Increment the counter for columns from append1_pr
    col_counter <- col_counter + 1
  }
  
  
  for (i in (nDP - 1):(1)) {
    
    ## start at stage nDP:
    message("Convergence Strata ", strata, " starting prediction from stage ", i)
    
    ####
    #### run for strata > 1
    ####
    
    ####### if the strata isn't 1, we perform this step for patients whose next stage (i + 1) is still in strata 2
    # otherwise, if their next stage is in strata 1, we skip this step for now
    # instead, we will append using their output probabilities
    
    ## create eligibility vector for if a patient's next strata is in 1
    # if i's strata == strata, and i + 1's strata is in == strata - 1, and the value of T isn't missing
    # we indicate a T. Otherwise, we indicate a F
    
    
    # Filter for the current stage
    current_stage <- long_data %>%
      filter(stage == i) %>%
      select(subj.id, paste0("strata", strata), T)
    
    # Filter for the next stage (i + 1)
    next_stage <- long_data %>%
      filter(stage == i + 1) %>%
      select(subj.id, paste0("strata", strata - 1), T)
    
    # Rename columns for clarity in join
    colnames(current_stage)[2] <- "current_strata"
    colnames(current_stage)[3] <- "current_T"
    colnames(next_stage)[2] <- "next_strata"
    colnames(next_stage)[3] <- "next_T"
    
    # Join the current stage and next stage data on subj.id
    merged_data <- left_join(current_stage, next_stage, by = "subj.id")
    
    # Create the logical vector based on the conditions
    nextstrat_dif <- with(merged_data, current_strata == 1 & next_strata == 1 & !is.na(current_T) & !is.na(next_T))
    
    # we also create a logical vector if the current strata == 1 and the next strata == 0 AKA both in current strata
    # Create the logical vector based on the conditions
    nextstrat_same <- with(merged_data, current_strata == 1 & next_strata == 0 & !is.na(current_T) & !is.na(next_T))
    
    
    ### NOTE: we only predict for those patients who have the same next strata
    ###
    x = get_all_vars(updated_formula, data)[which(nextstrat_same == 1), ]
    
    ####
    #### end for strata >1
    ####
    
    
    if (dim(x)[1] != 0){
      
      ## remove stage suffixto use in prediction
      new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
      colnames(x) <- new_col_names
      last.stage.pred <- .Predict(object = prev.iteration@FinalForest,
                                  newdata = x,
                                  params = params,
                                  findOptimal = T)
      
      ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
      # Iterate over each stage of interest
      
      # Extract optimal treatments for the current stage
      optimal_treatments <- last.stage.pred$optimal@optimalTx
    
    
      ####
      #### run for strata > 1
      ####
      
      # also update "A" column
      long_data$A[long_data$stage == (i+1)][which(nextstrat_same == 1)] <- optimal_treatments
      
      
      ####
      #### end for strata > 1
      ####
      
    
      
      ##########################################
      ## putting optimal predictions together ##
      ##########################################
    
    
      ###
      ### for  strata > 1
      ###
      
      # Initialize the shifted probability matrix
      ## we overwrite the eligible patients
      
      shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(nextstrat_same))
      
      ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
      ## do this only for eligible patients
      shiftedprob1[, which(nextstrat_same == 1) ] <- t(last.stage.pred$optimal@optimalY)
      
      ###
      ### end strata > 1
      ###
      
      ###################################################
      ################################################# NEXT START HERE
      ###############################################################
    
    
    
    }
    
    
    
    
    