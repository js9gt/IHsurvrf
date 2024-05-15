

setwd("~/survrf/Scripts/IHsurvrf")
source("R/IH.dtrSurv.R")



######## this is a function to assess convergence using results from IH.dtrSurv.R which fits the initial forest


IHdtrConv <- function(data,

                      ## the input forest results from the previous iteration
                      prev.iteration = NULL, ...){

  ## using the long data that was input
  data <- data

  prev.iteration <- prev.iteration


##### for all the previous setups for parameters and such, we want to use the values from the current environment

  # Now, we start the second iteration
  ### we use the final forest fit from IH.dtrSurv.R in res@FinalForest
  ### then, using the data at the FINAL STAGE, we feed this through the forest to obtain the predicted optimal survival curves for those patients
  ### then, using the second to last stage, we start the appending & pooling process



  message("Algorithm Iteration 2: starting from stage ", nDP)

  ## first, using the data at the last stage, predict the outcome using the final forest

  ### we need to re-construct the data to be in a predictable format
  # Extract the response variable from the formula
  response_var <- as.character(formula(prev.iteration@FinalForest@model)[[2]][-1])

  # Append the current stage number to each term
  response_with_stage <- paste0(response_var, "_", nDP)

  # Extract the terms of the formula excluding the response variable
  terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")

  # Append the current stage number to each term
  terms_with_stage <- paste0(terms, "_", nDP)

  # Reconstruct the formula
  updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))

  x = get_all_vars(updated_formula, data)

  ### for the new data, we need to get rid of the stage labels since our model doesn't have any stage labels:
  # Remove "_stage" suffix from all column names
  new_col_names <- gsub(paste0("_", nDP, "$"), "", colnames(x))
  colnames(x) <- new_col_names

  ## the arguments for predicting include only the data from the last stage, and the information from the final forest
  args <- list(prev.iteration, newdata = x)

  ## feeds this into predict() function which is defined in class_DTRSurv.R
  ## acts on objects of class DTRSurv
  ## the new data we are feeding through the forest is the covariates at the current stage


  ## this then subsets to objects of class DTRSurvStep, and calls .Predict() in class_IH.DTRSurv.R
  ## .Predict() on objects of DTRSurvStep is defined in class_IH.DTRSurvStep.R
  ## this predicts based on the results of the final forest
  ## this subsets to the "SurvRF" object to act on which is defined in class_IH.SurvRF.R

  ## follow the predicted optimal policy based on the input policy
  last.stage.pred <- do.call(predict, args)


  ### then, w/ the predicted survival curve, we would like to bring in a next new stage & append their survival curve, using the pooling


  ### now we need to pool the data together from the results. Inputs:
  ## model
  ## data: changed with whole data with all patients and all stages reformatted
  ## priorStep shifted probabilities: does each patient need to have their own matrix (nTime rows, 1 column)
  ## same params
  ## txName will just be "A"-- as we no longer separate things by stage
  ## mTry: calculated based on covariates
  ## sampleSIze == 1 to use all the pooled patients

  # Extract variable names from the data frame and remove "subj.id" and "rep.id" if they exist in the setdiff()
  ## then, extract all the unique portions before the underscore to be used as the common variable name
  variables <- unique(sub("_(\\d+)", "", setdiff(colnames(data), c("subj.id", "rep.id"))))


  ## reformatting the data in long form
  # Convert from wide to long format
  ## pivot_longer is function from tidyr


  long_data <-  data %>%
    pivot_longer(cols = starts_with(variables),
                 names_sep = "_",
                 names_to = c(".value", "stage"),
                 names_prefix = "",
                 values_to = "value",
                 values_drop_na = FALSE) %>%
    mutate(stage = as.integer(stage)) %>%
    mutate(A.iter1 = A) %>%
    arrange(subj.id, stage) %>%

    ## create a column called "A.opt.HC" to hold the optimal actions calculated from the last stage prediction using the overall forest trained in iteration 1
    mutate(A.opt.iter2 = NA)


  ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
  # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage-- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # Update A.opt.iter2 column for patients in the current stage for tracking
    long_data$A.opt.iter2[long_data$stage == nDP][which(eligibility == 1)] <- optimal_treatments

    # also update "A" column
    long_data$A[long_data$stage == nDP][which(eligibility == 1)] <- optimal_treatments

    ##########################################
    ## putting optimal predictions together ##
    ##########################################


  # Initialize the shifted probability matrix
  ## we overwrite the eligible patients

  shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

  ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
  ## do this only for eligible patients
  shiftedprob1[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)



  ##### NOTE: this step is equivalent to shifting by 0, essentially we just transform the optimal survival probabilities into probability mass vectors
  ## this step is essential to match with the probabilities that we've just computed
  ####################
  ####################
  ####################

  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprob1 <- shiftedprob1 - rbind(shiftedprob1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprob1[abs(shiftedprob1) < 1e-8] <- 0.0

  ### now, we need to create 0s for all other columns where patients were not eligible
  ## only overwrite the columns that had the optimal probabilities with the columns where the patients were eligible
  prev.op1  <- matrix(0, nrow = nTimes, ncol = length(eligibility))
  prev.op1 <- shiftedprob1


  ####################
  ####################
  ####################
  ####################

  ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage


  ## add the shifted optimals into the pooled optimals together

  # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled <- matrix(0, nrow = nTimes, ncol = ncol(prev.op1)*nDP)


  # Define the index to insert columns from the previous optimal
  ## this should be for the last stage
  insert_index <- seq(from = nDP, to = ncol(pr_pooled), by = nDP)

  # Initialize a counter for columns from append1_pr
  col_counter <- 1

  # Loop over each index to insert columns from append1_pr into combined_results
  for (w in seq_along(insert_index)) {
    # Insert columns from append1_pr into combined_results
    pr_pooled[,insert_index[w]] <-
      prev.op1[, w]

    # Increment the counter for columns from append1_pr
    col_counter <- col_counter + 1
  }


###### looping through the rest of the stages:

  for (i in (nDP - 1):1) {
    response_var <- as.character(formula(prev.iteration@FinalForest@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", i)
    terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", i)
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
    x = get_all_vars(updated_formula, data)
    new_col_names <- gsub(paste0("_", i, "$"), "", colnames(x))
    colnames(x) <- new_col_names
    args <- list(prev.iteration, newdata = x)
    last.stage.pred <- do.call(predict, args)

    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # Update A.opt.iter2 column for patients in the current stage for tracking
    long_data$A.opt.iter2[long_data$stage == i][which(eligibility == 1)] <- optimal_treatments

    # also update "A" column
    long_data$A[long_data$stage == i][which(eligibility == 1)] <- optimal_treatments


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients
    shiftedprob1[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)



    ##### NOTE: this step is equivalent to shifting by 0, essentially we just transform the optimal survival probabilities into probability mass vectors
    ## this step is essential to match with the probabilities that we've just computed
    ####################
    ####################
    ####################

    ## calculate the change in survival probability at each consecutive pair of time points
    ## subtracting each row of survShifted from the row above it
    ## append a row of 0s at the end to align with matrix dimensions
    ## probability mass vector representing the change in survival probabilities at each time point
    shiftedprob1 <- shiftedprob1 - rbind(shiftedprob1[-1L,], 0.0)

    ## sets very small values in pr to 0

    shiftedprob1[abs(shiftedprob1) < 1e-8] <- 0.0

    ### now, we need to create 0s for all other columns where patients were not eligible
    ## only overwrite the columns that had the optimal probabilities with the columns where the patients were eligible
    prev.op1  <- matrix(0, nrow = nTimes, ncol = length(eligibility))
    prev.op1 <- shiftedprob1


    ####################
    ####################
    ####################
    ####################

    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage


    ## add the shifted optimals into the pooled optimals together



    # Define the index to insert columns from the previous optimal
    ## this should be for the last stage
    insert_index <- seq(from = i, to = ncol(pr_pooled), by = nDP)

    # Initialize a counter for columns from append1_pr
    col_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled[,insert_index[j]] <-
        prev.op1[, j]

      # Increment the counter for columns from append1_pr
      col_counter <- col_counter + 1
    }
  }

  ## end of the for loop

  ### now, we use the full data as well as the input probabilities to fit a new random forest
  ## we also use the new optimal actions the patient received at each stage

  ## after appending, we fit another pooled random forest, including the data from nDP - 3

  conv1 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## convergence forest uses all the data
    data = long_data,
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input: these mTry values are all the same so I just use an arbitrary one
    mTry = mTry[[nDP - 3]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## this is used as TRUE for processing the treatment levels
    pool1 = TRUE,
    ## as TRUE allows for inputting the probability matrix to be used instead of creating one separately in the function
    appendstep1 = TRUE,
    inputpr = pr_pooled
  )

  # Set the column name in long_data
  long_data$A.final <- NA

  # Extract eligibility for the current forest
  eligibility_final <- conv1@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column
  long_data$A.final[eligibility_final] <- conv1@optimal@optimalTx
  long_data$A[eligibility_final] <- conv1@optimal@optimalTx

  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_final) )
  shiftedprobfinal[,eligibility_final] <-t(conv1@optimal@optimalY)



  ## get the final stage's area under the curve: we only look at the columns of the matrix in the last stage
  finalstagepr <- shiftedprobfinal[, seq(from = 1, to = ncol(shiftedprobfinal), by = nDP)]


  ## the result after appyling the area function to each column of the matrix of survival probabilities
  areas <- apply(finalstagepr, 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })

  ############# mean value results
  ## mean of maximum expected survival times across treatment levels & mean of maximum expected survival probabilities
  ## these are returned (not calculated) based on the maximum value for whatever the best treatment is
  ## ex) mean 1 (for trt1), mean2 (for trt2); we identify which the maximum is
  # value obtained from the first stage analysis

  ## after backwards recursion is complete, calculate the estimated value from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  valueTrain <- .meanValue(object = conv1)

  ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
  }

  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  conv_forest <- new(
    Class = "DTRSurv",
    "stageResults" = list(),
    "IHstageResults" = list(),
    "FinalForest" = conv1,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = areas
  )

  return(conv_forest)

}
