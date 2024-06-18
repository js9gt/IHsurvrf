

#setwd("~/survrf/Scripts/IHsurvrf")
#source("R/IH.dtrSurv.R")



######## this is a function to assess convergence using results from IH.dtrSurv.R which fits the initial forest


IHdtrConv <- function(data,

                      ## the input forest results from the previous iteration
                      prev.iteration = NULL, nDP, params, nTimes, models, mTry, ...){

  ## using the long data that was input
  data <- data

  prev.iteration <- prev.iteration

  nDP <- nDP

  params <- params

  nTimes <- nTimes

  models <- models

  mTry <- mTry

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

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3

#    xaxis <- params@timePoints
#    ## nDP
#    y1 <- last.stage.pred$optimal@optimalY[1,]
#
#    #eligibility
#
#    plot(xaxis, y1)
#    title("Convergence 3: stage nDP prediction")


    #############################
    #############################
    #############################
    #############################

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
  shiftedprob_done <- shiftedprob1 - rbind(shiftedprob1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprob_done[abs(shiftedprob_done) < 1e-8] <- 0.0

  ### now, we need to create 0s for all other columns where patients were not eligible
  ## only overwrite the columns that had the optimal probabilities with the columns where the patients were eligible
  prev.op1  <- matrix(0, nrow = nTimes, ncol = length(eligibility))
  prev.op1 <- shiftedprob_done


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


  ############## now, for stage 24, we work on appending
  ### to do: only do this for eligible patients for the new stage nDP - 3, and then perform shifting

  # identify individuals with complete data in new stage

  ## prepares model frame using the formula (input) and the data (input)-- NOTE we use the original model input
  ## models = Surv(Y,D)~X + A
  ## missing values are not specifically handled (not omitted)
  ## only include the variables that are used in the model formula mod
  ## If there are variables in your data that are not used in the formula, they won't be included in the resulting x.

  x_append1 <- stats::model.frame(formula = models,
                                  ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                  data = long_data %>% filter(stage == (nDP -1)) %>% dplyr::select(-matches("^A\\.")),
                                  na.action = na.pass)


  ## we want to exclude the A.opt.HC column and A.pool1 column
  elig_append1 <- stats::complete.cases(x_append1)

  # extract response and delta from model frame

  ## extract survival response
  response_append1 <- stats::model.response(data = x_append1)

  ## extract censoring indicator (delta) from the second column of the "response" data
  ## "L" is used to indicate that 2 is an integer
  delta_append1 <- response_append1[, 2L]

  ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
  response_append1 <- response_append1[, 1L]

  # remove response from x

  ## if first column of the model frame (x) is the response variable, remove this column
  ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
  if (attr(x = terms(x = models), which = "response") == 1L) {
    x_append1 <- x_append1[,-1L, drop = FALSE]
  }

  ## marks zeroed survival times and updates eligibility

  # responses that are zero (effectively) indicate censored at a previous stage

  ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
  zeroed <- abs(x = response_append1) < 1e-8

  ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

  elig_append1 <- elig_append1 & !zeroed

  ## if there are no eligible cases, output an error message

  if (sum(elig_append1) == 0L)
    stop("no cases have complete data", call. = FALSE)

  ## displays message indicating the number of cases that are still eligible for analysis at this stage

  message("cases in stage: ", sum(elig_append1))


  ## retrieve the number of timepoints from the params object
  ## defined in class_TimeInfo.R

  nTimes <- .NTimes(object = params)

  # create an empty matrix for all previously eligible cases

  ## initializes a survival matrix called "survMatrix" with 0's

  survMatrix <- matrix(data = 0.0,
                       nrow = nTimes,
                       ncol = nrow(x = x_append1))

  ## sets the first row to 1.0

  survMatrix[1L,] <- 1.0

  ## updates survMatrix with survival functions estimated from the previous step
  ## priorStep: A DTRSurvStep object. The analysis from a previous step

  # retrieve estimated OPTIMAL survival function from previous step
  ## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
  ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
  ## transpose the output of .OptimalY to align with the structure
  ## .OptimalY defined in class_Optimal.R
  ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
  ## optimal slot is  object of class optimal
  ## .OptimalY acts to return the optimalY slot



  ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
  ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
  ## with columns with values that are non-0
  survMatrix[, which(eligibility == 1)] <-  shiftedprob1[, colSums(shiftedprob1, na.rm = T) != 0]


  # shift the survival function down in time (T_i - Tq) and
  # transform to a probability mass vector for only those
  # eligible for this stage
  # .shiftMat is an internal function defined in shiftMat.R

  ## transforms the survival functions in survMatrix to probability mass format for the current stage
  ## shifts the survival function based on the observed survival times
  ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

  append1_pr_1 <- .shiftMat(
    timePoints = .TimePoints(object = params),

    ## extracts columns from survMatrix corresponding to cases that are eligible
    ## this is a matrix matrix where each column represents survival function for an individual
    survMatrix = survMatrix[, elig_append1, drop = FALSE],

    ## extracts survival times corresponding to eligible cases
    ## this is how much to shift survival function for each individual
    shiftVector = response_append1[elig_append1],

    ## probably transforming survival times into probabilities?
    surv2prob = TRUE
  )

  #############################
  #############################
  #############################
  ####### VISUALIZATION #######

  ## shifted probability output after adding new point for nDP - 3


  ## nDP
# y1 <- .shiftMat(
#   timePoints = .TimePoints(object = params),

#   ## extracts columns from survMatrix corresponding to cases that are eligible
#   ## this is a matrix matrix where each column represents survival function for an individual
#   survMatrix = survMatrix[, elig_append1, drop = FALSE],

#   ## extracts survival times corresponding to eligible cases
#   ## this is how much to shift survival function for each individual
#   shiftVector = response_append1[elig_append1],

#   ## probably transforming survival times into probabilities?
#   surv2prob = FALSE
# )[, 1]



# plot(xaxis, y1)
# title("Convergence 3: stage nDP - 1 appending to nDP prediction")

  #############################
  #############################
  #############################
  #############################

  ## sets very small values in pr to 0

  append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

  ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

  ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
  ## AKA we overwrite these shifted probabilities for eligible patients
  ## number of columns is the number of patients for this stage (nDP - 2)

  append1_pr <- matrix(0, nrow = nTimes, ncol = length(elig_append1))
  ## we overwrite the patients who are eligible in new stage (nDP - 3)
  append1_pr[,elig_append1] <-append1_pr_1

  ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations

  # Define the index to insert columns from the previous optimal
  ## this should be for the last stage
  insert_index <- seq(from = (nDP-1), to = ncol(pr_pooled), by = nDP)

  # Initialize a counter for columns from append1_pr
  col_counter <- 1

  # Loop over each index to insert columns from append1_pr into combined_results
  for (w in seq_along(insert_index)) {
    # Insert columns from append1_pr into combined_results
    pr_pooled[,insert_index[w]] <-
      append1_pr[, w]

    # Increment the counter for columns from append1_pr
    col_counter <- col_counter + 1
  }

###### looping through the rest of the stages:
  ## starting with stage 23, we append stage 23 survival to optimized survival at stage 24

  for (i in (nDP - 2):(nDP-4)) {

    ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
    ### if the current iteration is nDP - 2, we want to predict for nDP - 1
    ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

    ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
    response_var <- as.character(formula(prev.iteration@FinalForest@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", (i+1))
    terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", (i+1))
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
    x = get_all_vars(updated_formula, data)
    ## remove stage suffixto use in prediction
    new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
    colnames(x) <- new_col_names
    args <- list(prev.iteration, newdata = x)
    last.stage.pred <- do.call(predict, args)

    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # also update "A" column
    long_data$A[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3

#    xaxis <- params@timePoints
#    ## nDP
#    y1 <- last.stage.pred$optimal@optimalY[1,]
#
#
#    plot(xaxis, y1)
#    title("Convergence 3: stage nDP - 1 prediction")


    #############################
    #############################
    #############################
    #############################


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients
    shiftedprob1[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)




    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage

    ####################
    #################### appending the previous one
    ####################
    ####################

    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)


    ## we want to exclude the A.opt.HC column and A.pool1 column
    elig_append1 <- stats::complete.cases(x_append1)

    # extract response and delta from model frame

    ## extract survival response
    response_append1 <- stats::model.response(data = x_append1)

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]

    # remove response from x

    ## if first column of the model frame (x) is the response variable, remove this column
    ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
    if (attr(x = terms(x = models), which = "response") == 1L) {
      x_append1 <- x_append1[,-1L, drop = FALSE]
    }

    ## marks zeroed survival times and updates eligibility

    # responses that are zero (effectively) indicate censored at a previous stage

    ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
    zeroed <- abs(x = response_append1) < 1e-8

    ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

    elig_append1 <- elig_append1 & !zeroed

    ## if there are no eligible cases, output an error message

    if (sum(elig_append1) == 0L)
      stop("no cases have complete data", call. = FALSE)

    ## displays message indicating the number of cases that are still eligible for analysis at this stage

    message("cases in stage: ", sum(elig_append1))


    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases

    ## initializes a survival matrix called "survMatrix" with 0's

    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = nrow(x = x_append1))

    ## sets the first row to 1.0

    survMatrix[1L,] <- 1.0

    ## updates survMatrix with survival functions estimated from the previous step
    ## priorStep: A DTRSurvStep object. The analysis from a previous step

    # retrieve estimated OPTIMAL survival function from previous step
    ## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
    ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
    ## transpose the output of .OptimalY to align with the structure
    ## .OptimalY defined in class_Optimal.R
    ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
    ## optimal slot is  object of class optimal
    ## .OptimalY acts to return the optimalY slot



    ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
    ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
    ## with columns with values that are non-0

    survMatrix[, which(eligibility == 1)] <-  shiftedprob1[, colSums(shiftedprob1, na.rm = T) != 0]


    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

    append1_pr_1 <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig_append1, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response_append1[elig_append1],

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3


#    ## nDP
#    y1 <- .shiftMat(
#      timePoints = .TimePoints(object = params),
#
#      ## extracts columns from survMatrix corresponding to cases that are eligible
#      ## this is a matrix matrix where each column represents survival function for an individual
#      survMatrix = survMatrix[, elig_append1, drop = FALSE],
#
#      ## extracts survival times corresponding to eligible cases
#      ## this is how much to shift survival function for each individual
#      shiftVector = response_append1[elig_append1],
#
#      ## probably transforming survival times into probabilities?
#      surv2prob = FALSE
#    )[, 1]
#
#
#
#    plot(xaxis, y1)
#    title("Convergence 2: stage 1 appending to stage 2 prediction")
#

    #############################
    #############################
    #############################
    #############################

    ## sets very small values in pr to 0

    append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

    ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

    ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
    ## AKA we overwrite these shifted probabilities for eligible patients
    ## number of columns is the number of patients for this stage (nDP - 2)

    append1_pr <- matrix(0, nrow = nTimes, ncol = length(elig_append1))
    ## we overwrite the patients who are eligible in new stage (nDP - 3)
    append1_pr[,elig_append1] <-append1_pr_1



    # Define the index to insert columns from the previous optimal
    ## this should be for the last stage
    insert_index <- seq(from = i, to = ncol(pr_pooled), by = nDP)

    # Initialize a counter for columns from append1_pr
    col_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled[,insert_index[j]] <-
        append1_pr[, j]

      # Increment the counter for columns from append1_pr
      col_counter <- col_counter + 1
    }
  }

  ### now we fit a random forest for stages 21-25 with the data pooled together based on the output stubs and double stubs
  ## after appending, we fit another pooled random forest, including the data from nDP - 3

  s21.25 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## convergence forest uses all the data
    data = long_data[long_data$stage >= (nDP - 4) & long_data$stage <= nDP, ],
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
    ## the input probability should be only the last 5 of each stage
    inputpr = pr_pooled[, unlist(sapply(seq(21, 2500, by = 25), function(x) x:(x + 4)))]
  )


  # Set the column name in long_dai = ta
  long_data$A.pool5.1 <- NA

  # Extract eligibility for the current forest
  eligibility_s21.25 <- s21.25@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column
  long_data$A.pool5.1[which(eligibility_s21.25 == 1)] <- s21.25@optimal@optimalTx
  long_data$A[which(eligibility_s21.25 == 1)] <- s21.25@optimal@optimalTx

  # Initialize output probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobs21.25_1 <- matrix(0, nrow = nTimes, ncol = length(eligibility_s21.25) )
  shiftedprobs21.25_1[,which(eligibility_s21.25 == 1)] <-t(s21.25@optimal@optimalY)

  ## we need to convert these shifted probabilities into probability differences
  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprobs21.25 <- shiftedprobs21.25_1 - rbind(shiftedprobs21.25_1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprobs21.25[abs(shiftedprobs21.25) < 1e-8] <- 0.0


  ### now for stages 19, 18, ...15, we repeat the creation of the double stubs by feeding the predictions through this forest we just fit
  ### then, we will fit another forest with double stubs from stage 15-19 as well as the estimated optimal curves output from stage 20-25

  ###### looping through the rest of the stages:
  ## starting with stage 23, we append stage 23 survival to optimized survival at stage 24

  ## initialize matrix of probs to hold all 5 stages for all patients

  # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled5 <- matrix(0, nrow = nTimes, ncol = ncol(prev.op1)*5)


  colcounter <-5
  for (i in (nDP - 5):(nDP-9)) {

    ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
    ### if the current iteration is nDP - 2, we want to predict for nDP - 1
    ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

    ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
    response_var <- as.character(formula(s21.25@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", (i+1))
    terms <- attr(terms(s21.25@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", (i+1))
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
    x = get_all_vars(updated_formula, data)
    ## remove stage suffixto use in prediction
    new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
    colnames(x) <- new_col_names
    #args <- list(s20.25, newdata = x)
    last.stage.pred <- .Predict(object = s21.25,
                                newdata = x,
                                params = params,
                                findOptimal = T)

    ##########
    ########## TO DO: extract the optimal action and make a new column in the data for this; then work on the appending aspect to create new double stubs
    ##########


    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # Update A.opt.iter2 column for patients in the current stage for tracking
    #long_data$A.pool5.1[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    # also update "A" column
    long_data$A[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3

    #    xaxis <- params@timePoints
    #    ## nDP
    #    y1 <- last.stage.pred$optimal@optimalY[1,]
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 3: stage nDP - 1 prediction")


    #############################
    #############################
    #############################
    #############################


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob16.20 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients for a single stage which is stage 19, so every patient has 1 column
    shiftedprob16.20[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)




    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage

    ####################
    #################### appending the previous one
    ####################
    ####################

    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)


    ## we want to exclude the A.opt.HC column and A.pool1 column
    elig_append1 <- stats::complete.cases(x_append1)

    # extract response and delta from model frame

    ## extract survival response
    response_append1 <- stats::model.response(data = x_append1)

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]

    # remove response from x

    ## if first column of the model frame (x) is the response variable, remove this column
    ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
    if (attr(x = terms(x = models), which = "response") == 1L) {
      x_append1 <- x_append1[,-1L, drop = FALSE]
    }

    ## marks zeroed survival times and updates eligibility

    # responses that are zero (effectively) indicate censored at a previous stage

    ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
    zeroed <- abs(x = response_append1) < 1e-8

    ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

    elig_append1 <- elig_append1 & !zeroed

    ## if there are no eligible cases, output an error message

    if (sum(elig_append1) == 0L)
      stop("no cases have complete data", call. = FALSE)

    ## displays message indicating the number of cases that are still eligible for analysis at this stage

    message("cases in stage: ", sum(elig_append1))


    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases

    ## initializes a survival matrix called "survMatrix" with 0's

    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = nrow(x = x_append1))

    ## sets the first row to 1.0

    survMatrix[1L,] <- 1.0

    ## updates survMatrix with survival functions estimated from the previous step
    ## priorStep: A DTRSurvStep object. The analysis from a previous step

    # retrieve estimated OPTIMAL survival function from previous step
    ## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
    ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
    ## transpose the output of .OptimalY to align with the structure
    ## .OptimalY defined in class_Optimal.R
    ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
    ## optimal slot is  object of class optimal
    ## .OptimalY acts to return the optimalY slot



    ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
    ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
    ## with columns with values that are non-0

    survMatrix[, which(eligibility == 1)] <-  shiftedprob16.20[, colSums(shiftedprob16.20, na.rm = T) != 0]


    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

    append1_pr_1 <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig_append1, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response_append1[elig_append1],

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3


    #    ## nDP
    #    y1 <- .shiftMat(
    #      timePoints = .TimePoints(object = params),
    #
    #      ## extracts columns from survMatrix corresponding to cases that are eligible
    #      ## this is a matrix matrix where each column represents survival function for an individual
    #      survMatrix = survMatrix[, elig_append1, drop = FALSE],
    #
    #      ## extracts survival times corresponding to eligible cases
    #      ## this is how much to shift survival function for each individual
    #      shiftVector = response_append1[elig_append1],
    #
    #      ## probably transforming survival times into probabilities?
    #      surv2prob = FALSE
    #    )[, 1]
    #
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 2: stage 1 appending to stage 2 prediction")
    #

    #############################
    #############################
    #############################
    #############################

    ## sets very small values in pr to 0

    append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

    ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

    ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
    ## AKA we overwrite these shifted probabilities for eligible patients
    ## number of columns is the number of patients for this stage (nDP - 2)

    append1_pr <- matrix(0, nrow = nTimes, ncol = length(elig_append1))
    ## we overwrite the patients who are eligible in new stage (nDP - 3)
    append1_pr[,elig_append1] <-append1_pr_1



    # Define the index to insert columns from the previous optimal
    ## if i = nDP - 6 we should insert at the column 5
    ## if i = nDP - 7, we should insert at column 4
    insert_index <- seq(from = colcounter, to = ncol(pr_pooled5), by = 5)

    colcounter <- colcounter - 1

    # Initialize a counter for columns from append1_pr
    insert_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled5[,insert_index[j]] <-
        append1_pr[, j]

      # Increment the counter for columns from append1_pr
      insert_counter <- insert_counter + 1
    }
  }

  ### now, we construct another random forest, with stage 16-25 together, using the predicted probabilities frm 16-20
  ### we also use the output probabilities from 21-25 together

  # Initialize an empty matrix to store the combined results
  combined_matrix <- matrix(NA, nrow = nrow(pr_pooled5), ncol = 2 * ncol(pr_pooled5))

  # Define the number of columns to interleave at a time
  block_size <- 5

  # Interleave the columns in blocks of 5
  for (i in seq(0, ncol(pr_pooled5) - 1, by = block_size)) {
    combined_matrix[, (2 * i + 1):(2 * i + block_size)] <- pr_pooled5[, (i + 1):(i + block_size)]
    combined_matrix[, (2 * i + block_size + 1):(2 * i + 2 * block_size)] <- shiftedprobs21.25[, (i + 1):(i + block_size)]
  }



  ### now, we fit another pooled forest with stages 15-25

  ### now we fit a random forest for stages 20-25 with the data pooled together based on the output stubs and double stubs
  ## after appending, we fit another pooled random forest, including the data from nDP - 3

  s16.25 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## convergence forest uses all the data
    data = long_data[long_data$stage >= (nDP - 9) & long_data$stage <= nDP, ],
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
    ## the input probability should be only the last 5 of each stage
    inputpr = combined_matrix
  )

  ########
  ######## need to update the optimal actions, then we need to predict for stages 15, 14, 13, 12, 11
  ########

  # Set the column name in long_dai = ta
  long_data$A.pool5.2 <- NA

  # Extract eligibility for the current forest
  eligibility_s16.25 <- s16.25@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column
  long_data$A.pool5.2[which(eligibility_s16.25 == 1)] <- s16.25@optimal@optimalTx
  long_data$A[which(eligibility_s16.25 == 1)] <- s16.25@optimal@optimalTx

  # Initialize output probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobs16.25_1 <- matrix(0, nrow = nTimes, ncol = length(eligibility_s16.25) )
  shiftedprobs16.25_1[,which(eligibility_s16.25 == 1)] <-t(s16.25@optimal@optimalY)

  ## we need to convert these shifted probabilities into probability differences
  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprobs16.25 <- shiftedprobs16.25_1 - rbind(shiftedprobs16.25_1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprobs16.25[abs(shiftedprobs16.25) < 1e-8] <- 0.0


  ### now for stages 19, 18, ...15, we repeat the creation of the double stubs by feeding the predictions through this forest we just fit
  ### then, we will fit another forest with double stubs from stage 15-19 as well as the estimated optimal curves output from stage 20-25


  ## initialize matrix of probs to hold all 5 stages for all patients

  # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled5 <- matrix(0, nrow = nTimes, ncol = ncol(prev.op1)*5)

  colcounter <-5
  for (i in (nDP - 10):(nDP-14)) {

    ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
    ### if the current iteration is nDP - 2, we want to predict for nDP - 1
    ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

    ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
    response_var <- as.character(formula(s16.25@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", (i+1))
    terms <- attr(terms(s16.25@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", (i+1))
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
    x = get_all_vars(updated_formula, data)
    ## remove stage suffixto use in prediction
    new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
    colnames(x) <- new_col_names
    #args <- list(s20.25, newdata = x)
    last.stage.pred <- .Predict(object = s16.25,
                                newdata = x,
                                params = params,
                                findOptimal = T)

    ##########
    ########## TO DO: extract the optimal action and make a new column in the data for this; then work on the appending aspect to create new double stubs
    ##########


    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # Update A.opt.iter2 column for patients in the current stage for tracking
    #long_data$A.pool5.1[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    # also update "A" column
    long_data$A[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3

    #    xaxis <- params@timePoints
    #    ## nDP
    #    y1 <- last.stage.pred$optimal@optimalY[1,]
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 3: stage nDP - 1 prediction")


    #############################
    #############################
    #############################
    #############################


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob11.15 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients for a single stage which is stage 19, so every patient has 1 column
    shiftedprob11.15[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)




    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage

    ####################
    #################### appending the previous one
    ####################
    ####################

    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)


    ## we want to exclude the A.opt.HC column and A.pool1 column
    elig_append1 <- stats::complete.cases(x_append1)

    # extract response and delta from model frame

    ## extract survival response
    response_append1 <- stats::model.response(data = x_append1)

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]

    # remove response from x

    ## if first column of the model frame (x) is the response variable, remove this column
    ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
    if (attr(x = terms(x = models), which = "response") == 1L) {
      x_append1 <- x_append1[,-1L, drop = FALSE]
    }

    ## marks zeroed survival times and updates eligibility

    # responses that are zero (effectively) indicate censored at a previous stage

    ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
    zeroed <- abs(x = response_append1) < 1e-8

    ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

    elig_append1 <- elig_append1 & !zeroed

    ## if there are no eligible cases, output an error message

    if (sum(elig_append1) == 0L)
      stop("no cases have complete data", call. = FALSE)

    ## displays message indicating the number of cases that are still eligible for analysis at this stage

    message("cases in stage: ", sum(elig_append1))


    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases

    ## initializes a survival matrix called "survMatrix" with 0's

    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = nrow(x = x_append1))

    ## sets the first row to 1.0

    survMatrix[1L,] <- 1.0

    ## updates survMatrix with survival functions estimated from the previous step
    ## priorStep: A DTRSurvStep object. The analysis from a previous step

    # retrieve estimated OPTIMAL survival function from previous step
    ## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
    ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
    ## transpose the output of .OptimalY to align with the structure
    ## .OptimalY defined in class_Optimal.R
    ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
    ## optimal slot is  object of class optimal
    ## .OptimalY acts to return the optimalY slot



    ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
    ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
    ## with columns with values that are non-0

    survMatrix[, which(eligibility == 1)] <-  shiftedprob11.15[, colSums(shiftedprob11.15, na.rm = T) != 0]


    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

    append1_pr_1 <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig_append1, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response_append1[elig_append1],

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3


    #    ## nDP
    #    y1 <- .shiftMat(
    #      timePoints = .TimePoints(object = params),
    #
    #      ## extracts columns from survMatrix corresponding to cases that are eligible
    #      ## this is a matrix matrix where each column represents survival function for an individual
    #      survMatrix = survMatrix[, elig_append1, drop = FALSE],
    #
    #      ## extracts survival times corresponding to eligible cases
    #      ## this is how much to shift survival function for each individual
    #      shiftVector = response_append1[elig_append1],
    #
    #      ## probably transforming survival times into probabilities?
    #      surv2prob = FALSE
    #    )[, 1]
    #
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 2: stage 1 appending to stage 2 prediction")
    #

    #############################
    #############################
    #############################
    #############################

    ## sets very small values in pr to 0

    append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

    ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

    ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
    ## AKA we overwrite these shifted probabilities for eligible patients
    ## number of columns is the number of patients for this stage (nDP - 2)

    append1_pr <- matrix(0, nrow = nTimes, ncol = length(elig_append1))
    ## we overwrite the patients who are eligible in new stage (nDP - 3)
    append1_pr[,elig_append1] <-append1_pr_1



    # Define the index to insert columns from the previous optimal
    ## if i = nDP - 6 we should insert at the column 5
    ## if i = nDP - 7, we should insert at column 4
    insert_index <- seq(from = colcounter, to = ncol(pr_pooled5), by = 5)

    colcounter <- colcounter - 1

    # Initialize a counter for columns from append1_pr
    insert_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled5[,insert_index[j]] <-
        append1_pr[, j]

      # Increment the counter for columns from append1_pr
      insert_counter <- insert_counter + 1
    }
  }

  #### now we need to combine the predictions from stages 11-15 with the output from the forest for 16-25
  ## using pr_pooled5, we first put 5, then 10 from the column for 16-25
  # Initialize an empty matrix to hold the combined results
  combined_matrix <- matrix(ncol = ncol(pr_pooled5) + ncol(shiftedprobs16.25), nrow = nrow(pr_pooled5))

  # Define the column index counters for pr_pooled5 and shiftedprobs16.25
  pr_idx <- 1
  shifted_idx <- 1

  # Define the column index counter for the combined_matrix
  combined_idx <- 1

  # Number of iterations (blocks) needed
  num_iterations <- ncol(pr_pooled5)/5

  # Interleave columns
  for (i in 1:num_iterations) {
    # Get 5 columns from pr_pooled5
    combined_matrix[, combined_idx:(combined_idx + 4)] <- pr_pooled5[, pr_idx:(pr_idx + 4)]
    combined_idx <- combined_idx + 5
    pr_idx <- pr_idx + 5

    # Get 10 columns from shiftedprobs16.25
    combined_matrix[, combined_idx:(combined_idx + 9)] <- shiftedprobs16.25[, shifted_idx:(shifted_idx + 9)]
    combined_idx <- combined_idx + 10
    shifted_idx <- shifted_idx + 10
  }

  ## now, we fit another pooled forest with all 15 stages together

  ### now, we fit another pooled forest with stages 15-25

  ### now we fit a random forest for stages 20-25 with the data pooled together based on the output stubs and double stubs
  ## after appending, we fit another pooled random forest, including the data from nDP - 3

  s11.25 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## convergence forest uses all the data
    data = long_data[long_data$stage >= (nDP - 14) & long_data$stage <= nDP, ],
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
    ## the input probability should be only the last 5 of each stage
    inputpr = combined_matrix
  )

  ########
  ######## need to update the optimal actions, then we need to predict for stages 15, 14, 13, 12, 11
  ########

  # Set the column name in long_dai = ta
  long_data$A.pool5.3 <- NA

  # Extract eligibility for the current forest
  eligibility_s11.25 <- s11.25@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column
  long_data$A.pool5.3[which(eligibility_s11.25 == 1)] <- s11.25@optimal@optimalTx
  long_data$A[which(eligibility_s11.25 == 1)] <- s11.25@optimal@optimalTx

  # Initialize output probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobs11.25_1 <- matrix(0, nrow = nTimes, ncol = length(eligibility_s11.25) )
  shiftedprobs11.25_1[,which(eligibility_s11.25 == 1)] <-t(s11.25@optimal@optimalY)

  ## we need to convert these shifted probabilities into probability differences
  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprobs11.25 <- shiftedprobs11.25_1 - rbind(shiftedprobs11.25_1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprobs11.25[abs(shiftedprobs11.25) < 1e-8] <- 0.0


  ### now for stages 10, 9, 8, 7, 6, we repeat the creation of the double stubs by feeding the predictions through this forest we just fit
  ### then, we will fit another forest with double stubs from stage 15-19 as well as the estimated optimal curves output from stage 20-25

  ## initialize matrix of probs to hold all 5 stages for all patients

  # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled5 <- matrix(0, nrow = nTimes, ncol = ncol(prev.op1)*5)

  colcounter <-5
  for (i in (nDP - 15):(nDP-19)) {

    ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
    ### if the current iteration is nDP - 2, we want to predict for nDP - 1
    ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

    ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
    response_var <- as.character(formula(s11.25@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", (i+1))
    terms <- attr(terms(s11.25@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", (i+1))
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
    x = get_all_vars(updated_formula, data)
    ## remove stage suffixto use in prediction
    new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
    colnames(x) <- new_col_names
    #args <- list(s20.25, newdata = x)
    last.stage.pred <- .Predict(object = s11.25,
                                newdata = x,
                                params = params,
                                findOptimal = T)

    ##########
    ########## TO DO: extract the optimal action and make a new column in the data for this; then work on the appending aspect to create new double stubs
    ##########


    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # Update A.opt.iter2 column for patients in the current stage for tracking
    #long_data$A.pool5.1[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    # also update "A" column
    long_data$A[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3

    #    xaxis <- params@timePoints
    #    ## nDP
    #    y1 <- last.stage.pred$optimal@optimalY[1,]
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 3: stage nDP - 1 prediction")


    #############################
    #############################
    #############################
    #############################


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob6.10 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients for a single stage which is stage 19, so every patient has 1 column
    shiftedprob6.10[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)




    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage

    ####################
    #################### appending the previous one
    ####################
    ####################

    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)


    ## we want to exclude the A.opt.HC column and A.pool1 column
    elig_append1 <- stats::complete.cases(x_append1)

    # extract response and delta from model frame

    ## extract survival response
    response_append1 <- stats::model.response(data = x_append1)

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]

    # remove response from x

    ## if first column of the model frame (x) is the response variable, remove this column
    ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
    if (attr(x = terms(x = models), which = "response") == 1L) {
      x_append1 <- x_append1[,-1L, drop = FALSE]
    }

    ## marks zeroed survival times and updates eligibility

    # responses that are zero (effectively) indicate censored at a previous stage

    ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
    zeroed <- abs(x = response_append1) < 1e-8

    ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

    elig_append1 <- elig_append1 & !zeroed

    ## if there are no eligible cases, output an error message

    if (sum(elig_append1) == 0L)
      stop("no cases have complete data", call. = FALSE)

    ## displays message indicating the number of cases that are still eligible for analysis at this stage

    message("cases in stage: ", sum(elig_append1))


    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases

    ## initializes a survival matrix called "survMatrix" with 0's

    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = nrow(x = x_append1))

    ## sets the first row to 1.0

    survMatrix[1L,] <- 1.0

    ## updates survMatrix with survival functions estimated from the previous step
    ## priorStep: A DTRSurvStep object. The analysis from a previous step

    # retrieve estimated OPTIMAL survival function from previous step
    ## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
    ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
    ## transpose the output of .OptimalY to align with the structure
    ## .OptimalY defined in class_Optimal.R
    ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
    ## optimal slot is  object of class optimal
    ## .OptimalY acts to return the optimalY slot



    ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
    ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
    ## with columns with values that are non-0

    survMatrix[, which(eligibility == 1)] <-  shiftedprob6.10[, colSums(shiftedprob6.10, na.rm = T) != 0]


    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

    append1_pr_1 <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig_append1, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response_append1[elig_append1],

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3


    #    ## nDP
    #    y1 <- .shiftMat(
    #      timePoints = .TimePoints(object = params),
    #
    #      ## extracts columns from survMatrix corresponding to cases that are eligible
    #      ## this is a matrix matrix where each column represents survival function for an individual
    #      survMatrix = survMatrix[, elig_append1, drop = FALSE],
    #
    #      ## extracts survival times corresponding to eligible cases
    #      ## this is how much to shift survival function for each individual
    #      shiftVector = response_append1[elig_append1],
    #
    #      ## probably transforming survival times into probabilities?
    #      surv2prob = FALSE
    #    )[, 1]
    #
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 2: stage 1 appending to stage 2 prediction")
    #

    #############################
    #############################
    #############################
    #############################

    ## sets very small values in pr to 0

    append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

    ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

    ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
    ## AKA we overwrite these shifted probabilities for eligible patients
    ## number of columns is the number of patients for this stage (nDP - 2)

    append1_pr <- matrix(0, nrow = nTimes, ncol = length(elig_append1))
    ## we overwrite the patients who are eligible in new stage (nDP - 3)
    append1_pr[,elig_append1] <-append1_pr_1



    # Define the index to insert columns from the previous optimal
    ## if i = nDP - 6 we should insert at the column 5
    ## if i = nDP - 7, we should insert at column 4
    insert_index <- seq(from = colcounter, to = ncol(pr_pooled5), by = 5)

    colcounter <- colcounter - 1

    # Initialize a counter for columns from append1_pr
    insert_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled5[,insert_index[j]] <-
        append1_pr[, j]

      # Increment the counter for columns from append1_pr
      insert_counter <- insert_counter + 1
    }
  }

  #### now, we need to take 5 columns from pr_pooled5 and 15 columns from shiftedprobs11.25
  combined_matrix <- matrix(ncol = ncol(pr_pooled5) + ncol(shiftedprobs11.25), nrow = nrow(pr_pooled5))

  # Define the column index counters for pr_pooled5 and shiftedprobs16.25
  pr_idx <- 1
  shifted_idx <- 1

  # Define the column index counter for the combined_matrix
  combined_idx <- 1

  # Number of iterations (blocks) needed
  num_iterations <- ncol(pr_pooled5)/5

  for (i in 1:num_iterations) {
    # Get 5 columns from pr_pooled5
    combined_matrix[, combined_idx:(combined_idx + 4)] <- pr_pooled5[, pr_idx:(pr_idx + 4)]
    combined_idx <- combined_idx + 5
    pr_idx <- pr_idx + 5

    # Get 15 columns from shiftedprobs11.25
    combined_matrix[, combined_idx:(combined_idx + 14)] <- shiftedprobs11.25[, shifted_idx:(shifted_idx + 14)]
    combined_idx <- combined_idx + 15
    shifted_idx <- shifted_idx + 15
  }


  ## fit another random forest for stages 6-25
  s6.25 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## convergence forest uses all the data
    data = long_data[long_data$stage >= (nDP - 19) & long_data$stage <= nDP, ],
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
    ## the input probability should be only the last 5 of each stage
    inputpr = combined_matrix
  )


  ########
  ######## need to update the optimal actions, then we need to predict for stages 5, 4, 3, 2, 1
  ########

  # Set the column name in long_dai = ta
  long_data$A.pool5.4 <- NA

  # Extract eligibility for the current forest
  eligibility_s6.25 <- s6.25@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column
  long_data$A.pool5.4[which(eligibility_s6.25 == 1)] <- s6.25@optimal@optimalTx
  long_data$A[which(eligibility_s6.25 == 1)] <- s6.25@optimal@optimalTx

  # Initialize output probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobs6.25_1 <- matrix(0, nrow = nTimes, ncol = length(eligibility_s6.25) )
  shiftedprobs6.25_1[,which(eligibility_s6.25 == 1)] <-t(s6.25@optimal@optimalY)

  ## we need to convert these shifted probabilities into probability differences
  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprobs6.25 <- shiftedprobs6.25_1 - rbind(shiftedprobs6.25_1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprobs6.25[abs(shiftedprobs6.25) < 1e-8] <- 0.0

  # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled5 <- matrix(0, nrow = nTimes, ncol = ncol(prev.op1)*5)

  colcounter <-5
  for (i in (nDP - 20):(nDP-24)){

    ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
    ### if the current iteration is nDP - 2, we want to predict for nDP - 1
    ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

    ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
    response_var <- as.character(formula(s6.25@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", (i+1))
    terms <- attr(terms(s6.25@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", (i+1))
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = ", "), ") ~ ", paste(terms_with_stage, collapse = " + "))
    x = get_all_vars(updated_formula, data)
    ## remove stage suffixto use in prediction
    new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
    colnames(x) <- new_col_names
    #args <- list(s20.25, newdata = x)
    last.stage.pred <- .Predict(object = s6.25,
                                newdata = x,
                                params = params,
                                findOptimal = T)

    ##########
    ########## TO DO: extract the optimal action and make a new column in the data for this; then work on the appending aspect to create new double stubs
    ##########


    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # Update A.opt.iter2 column for patients in the current stage for tracking
    #long_data$A.pool5.1[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    # also update "A" column
    long_data$A[long_data$stage == (i+1)][which(eligibility == 1)] <- optimal_treatments

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3

    #    xaxis <- params@timePoints
    #    ## nDP
    #    y1 <- last.stage.pred$optimal@optimalY[1,]
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 3: stage nDP - 1 prediction")


    #############################
    #############################
    #############################
    #############################


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob1.5 <- matrix(0, nrow = nTimes, ncol = length(eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients for a single stage which is stage 19, so every patient has 1 column
    shiftedprob1.5[, which(eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)




    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage

    ####################
    #################### appending the previous one
    ####################
    ####################

    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)


    ## we want to exclude the A.opt.HC column and A.pool1 column
    elig_append1 <- stats::complete.cases(x_append1)

    # extract response and delta from model frame

    ## extract survival response
    response_append1 <- stats::model.response(data = x_append1)

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]

    # remove response from x

    ## if first column of the model frame (x) is the response variable, remove this column
    ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
    if (attr(x = terms(x = models), which = "response") == 1L) {
      x_append1 <- x_append1[,-1L, drop = FALSE]
    }

    ## marks zeroed survival times and updates eligibility

    # responses that are zero (effectively) indicate censored at a previous stage

    ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
    zeroed <- abs(x = response_append1) < 1e-8

    ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

    elig_append1 <- elig_append1 & !zeroed

    ## if there are no eligible cases, output an error message

    if (sum(elig_append1) == 0L)
      stop("no cases have complete data", call. = FALSE)

    ## displays message indicating the number of cases that are still eligible for analysis at this stage

    message("cases in stage: ", sum(elig_append1))


    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases

    ## initializes a survival matrix called "survMatrix" with 0's

    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = nrow(x = x_append1))

    ## sets the first row to 1.0

    survMatrix[1L,] <- 1.0

    ## updates survMatrix with survival functions estimated from the previous step
    ## priorStep: A DTRSurvStep object. The analysis from a previous step

    # retrieve estimated OPTIMAL survival function from previous step
    ## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
    ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
    ## transpose the output of .OptimalY to align with the structure
    ## .OptimalY defined in class_Optimal.R
    ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
    ## optimal slot is  object of class optimal
    ## .OptimalY acts to return the optimalY slot



    ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
    ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
    ## with columns with values that are non-0

    survMatrix[, which(eligibility == 1)] <-  shiftedprob1.5[, colSums(shiftedprob1.5, na.rm = T) != 0]


    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

    append1_pr_1 <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig_append1, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response_append1[elig_append1],

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )

    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3


    #    ## nDP
    #    y1 <- .shiftMat(
    #      timePoints = .TimePoints(object = params),
    #
    #      ## extracts columns from survMatrix corresponding to cases that are eligible
    #      ## this is a matrix matrix where each column represents survival function for an individual
    #      survMatrix = survMatrix[, elig_append1, drop = FALSE],
    #
    #      ## extracts survival times corresponding to eligible cases
    #      ## this is how much to shift survival function for each individual
    #      shiftVector = response_append1[elig_append1],
    #
    #      ## probably transforming survival times into probabilities?
    #      surv2prob = FALSE
    #    )[, 1]
    #
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 2: stage 1 appending to stage 2 prediction")
    #

    #############################
    #############################
    #############################
    #############################

    ## sets very small values in pr to 0

    append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

    ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

    ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
    ## AKA we overwrite these shifted probabilities for eligible patients
    ## number of columns is the number of patients for this stage (nDP - 2)

    append1_pr <- matrix(0, nrow = nTimes, ncol = length(elig_append1))
    ## we overwrite the patients who are eligible in new stage (nDP - 3)
    append1_pr[,elig_append1] <-append1_pr_1



    # Define the index to insert columns from the previous optimal
    ## if i = nDP - 6 we should insert at the column 5
    ## if i = nDP - 7, we should insert at column 4
    insert_index <- seq(from = colcounter, to = ncol(pr_pooled5), by = 5)

    colcounter <- colcounter - 1

    # Initialize a counter for columns from append1_pr
    insert_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled5[,insert_index[j]] <-
        append1_pr[, j]

      # Increment the counter for columns from append1_pr
      insert_counter <- insert_counter + 1
    }
  }


  #### now we need to combine the predictions from stages 1-5 with the output from the forest for 6-25
  ## using pr_pooled5, we first put 5, then 20 from the column for 6-25
  # Initialize an empty matrix to hold the combined results
  combined_matrix <- matrix(ncol = ncol(pr_pooled5) + ncol(shiftedprobs6.25), nrow = nrow(pr_pooled5))

  # Define the column index counters for pr_pooled5 and shiftedprobs16.25
  pr_idx <- 1
  shifted_idx <- 1

  # Define the column index counter for the combined_matrix
  combined_idx <- 1

  # Number of iterations (blocks) needed
  num_iterations <- ncol(pr_pooled5)/5

  # Interleave columns
  for (i in 1:num_iterations) {
    # Get 5 columns from pr_pooled5
    combined_matrix[, combined_idx:(combined_idx + 4)] <- pr_pooled5[, pr_idx:(pr_idx + 4)]
    combined_idx <- combined_idx + 5
    pr_idx <- pr_idx + 5

    # Get 20 columns from shiftedprobs6.25
    combined_matrix[, combined_idx:(combined_idx + 19)] <- shiftedprobs6.25[, shifted_idx:(shifted_idx + 19)]
    combined_idx <- combined_idx + 20
    shifted_idx <- shifted_idx + 20
  }



  ### now, we use the full data as well as the input probabilities to fit a new random forest
  ## we also use the new optimal actions the patient received at each stage

  ## we now fit a pooled forest with all stages together

  s1.25 <- .dtrSurvStep(
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
    inputpr = combined_matrix
  )

  # Set the column name in long_data
  long_data$A.final <- NA

  # Extract eligibility for the current forest
  eligibility_final <- s1.25@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column
  long_data$A.final[which(eligibility_final == 1)] <- s1.25@optimal@optimalTx
  long_data$A[which(eligibility_final == 1)] <- s1.25@optimal@optimalTx

  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(0, nrow = nTimes, ncol = length(eligibility_final) )
  shiftedprobfinal[,which(eligibility_final == 1)] <-t(s1.25@optimal@optimalY)

### once the forest has been trained with all the stages,
#  ## we take all the predicted optimals to train another forest
#
#  ### now, we use the full data as well as the input probabilities to fit a new random forest
#  ## we also use the new optimal actions the patient received at each stage
#
#  ## we now fit a pooled forest with all stages together
#
#  conv1 <- .dtrSurvStep(
#    ## use the model for pooled data
#    model = models,
#    ## convergence forest uses all the data
#    data = long_data,
#    priorStep = NULL,
#    params = params,
#    txName = "A",
#    ## set to the same as the input: these mTry values are all the same so I just use an arbitrary one
#    mTry = mTry[[nDP - 3]],
#    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
#    sampleSize = 1,
#    ## this is used as TRUE for processing the treatment levels
#    pool1 = TRUE,
#    ## as TRUE allows for inputting the probability matrix to be used instead of creating one separately in the function
#    appendstep1 = TRUE,
#    inputpr = shiftedprobpoolfinal
#  )
#
#  # Set the column name in long_data
#  long_data$A.final <- NA
#
#  # Extract eligibility for the current forest
#  eligibility_final <- conv1@eligibility
#
#  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
#  ## also update the "A" column
#  long_data$A.final[which(eligibility_final == 1)] <- conv1@optimal@optimalTx
#  long_data$A[which(eligibility_final == 1)] <- conv1@optimal@optimalTx
#
#  # Initialize the shifted probability matrix
#  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
#  ## we overwrite the eligible patients
#  shiftedprobfinal <- matrix(0, nrow = nTimes, ncol = length(eligibility_final) )
#  shiftedprobfinal[,which(eligibility_final == 1)] <-t(conv1@optimal@optimalY) %>% View()


  #############################
  #############################
  #############################
  ####### VISUALIZATION #######

     ## shifted probability output after adding new point for nDP - 3

#     xaxis <- params@timePoints
#     ## nDP
#     y1 <- shiftedprobfinal[,50]
#
#
#
#     plot(xaxis, y1)
#     title("Convergence 3: stage nDP output")


  #############################
  #############################
  #############################
  #############################


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

  ## after backwards recursion is complete, calculate the estimated value from the first stage: we only want to subset to eligible patients from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  ## for the mean value, we want to only subset for the patients from the last stage, not include the value across all stage

  mv_results <- matrix(0, nrow = length(eligibility_final), ncol = 2)
  mv_results[which(eligibility_final == 1), ] <- s1.25@valueAllTx$mean

  final_stagemv <- mv_results[seq(from = 1, to = nrow(mv_results), by = nDP), ]



  valueTrain <- mean(x = apply(

    ## access "mean" component of the valueAllTx slot which stores value of each tree for each treatment level
    X = final_stagemv,

    ##apply the "max" function across the margin of the mean matrix/array to apply across rows
    ## AKA for each row (individual),take maximum value across all columns (treatment options), then calculate the mean
    MARGIN = 1L,
    FUN = max))

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
    "FinalForest" = s1.25,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = areas,
    "n_it" = NA,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list()
  )

  return(conv_forest)

}
