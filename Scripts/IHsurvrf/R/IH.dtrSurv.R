
setwd("~/survrf/Scripts/IHsurvrf")
source("R/class_IH.DTRSurv.R")

source("R/IH.VerifyData.R")
source("R/IH.VerifyTxName.R")
source("R/IH.VerifyModels.R")
source("R/class_IH.CriticalValue.R")
source("R/class_IH.TimeInfo.R")
source("R/class_IH.TreeConditions.R")
source("R/IH.VerifySampleSize.R")
source("R/class_IH.Parameters.R")
source("R/area_function.R")
source("R/IH.dtrSurvConverge.R")
source("R/IH.dtrSurvConverge_otherstrat.R")
source("R/class_IH.DTRSurvRes.R")

library(tidyr)
library(dplyr)
library(survival)
library(IHsurvrf)

dyn.load("src/IHsurvrf.so")

IHdtrSurv <- function(data,
                         txName,
                         models,
                         ...,
                         timePoints = "uni",
                         nTimes = 200L,
                         tau = NULL,
                         criticalValue = "mean",
                         evalTime = NULL,
                         splitRule = NULL,
                         ERT = TRUE,
                         uniformSplit = NULL,
                         sampleSize = NULL,
                         replace = NULL,

                         ## fed into fortran via setUpBasics as "rs" (global param)
                         randomSplit = 0.2,
                         tieMethod = "random",
                         minEvent = 3L,
                         nodeSize = 6L,
                         nTree = 10L,
                         mTry = NULL,
                         stratifiedSplit = NULL,
                         stageLabel = ".",
                         ## number of strata we want to split into which fit different trees
                         nstrata = 2,

                         ## the number of days we want to go back by from tau to help define strata
                         ## the smaller, the more accurate the division of the percentage of events in each strata
                         windowsize = 50
                         ) {


  ## data verification

  # ensure that 'data' is provided as a data.frame or a matrix and does not
  # contain NaN values. If 'data' is appropriate, method returns a data.frame
  # object

  ## used in VerifyData.R script


  data <- .VerifyData(data = data)



  # total number of individuals in dataset
  nSamples <- nrow(x = data)

  # ensure that 'txName' is provided as a character or character vector and
  # that the provided names are present in 'data'. This input defines the
  # number of decision points for the analysis. If 'txName' is appropriate,
  # the object returned is the original input without modification.

  ## defined in VerifyTxName.R Script

  txName <- .VerifyTxName(txName = txName, data = data)

  ## The treatment variable name for
  #   each decision point. Each element corresponds to the respective decision point
  # element 1 = 1st decision; element 2 = 2nd decision, etc..

  # number of decision points in the analysis
  nDP <- length(x = txName)


  ## used in script VerifyModels.R
  ## this can be input as a list ofmodels, or a single formula (in which we probably assume common formula across all stages)

  stagemodels <- .VerifyModels(
    models = models,
    nDP = nDP,
    data = data,
    txName = txName,
    stageLabel = stageLabel
  )



  ## extracts the "response variable": matrix of survival response variables
  ## these were output from the .VerifyModels object
  response <- stagemodels$response

  ## extracts the "delta" variable
  del <- stagemodels$delta

  ## only including the models themselves & discarding "response" and "delta" components
  stagemodels <- stagemodels$models

  ## if models isn't a list already, wraps "models" in a list
  ## ensures uniform handling of "models" whether there is one or multiple decision points in an analysis

  # convert to list for single decision analyses for convenience
  if (!is.list(x = stagemodels))
    stagemodels <- list(stagemodels)

  params <- .parameters(
    timePoints = timePoints,
    tau = tau,
    nTimes = nTimes,

    ## response is created (not input)
    response = response,
    nTree = nTree,
    ERT = ERT,
    uniformSplit = uniformSplit,
    randomSplit = randomSplit,
    splitRule = splitRule,
    replace = replace,
    nodeSize = nodeSize,
    minEvent = minEvent,
    tieMethod = tieMethod,
    criticalValue = criticalValue,
    survivalTime = evalTime,

    ## nSamples is created (not input)
    nSamples = nSamples,
    stratifiedSplit = stratifiedSplit
  )


  # store basic information on fortran side

  ## .CriticalValueCriterion:
  ## if the object type is "mean", then return "surv.mean" if object type is "prob", return the "surv.prob"


  # retrieve index and fraction if survival type value estimator
  # set as 0's if mean estimator

  crit <- .CriticalValueCriterion(params)
  if (crit == "mean") {
    ## set index to 0 (integer)
    ind = 0L

    ## set fraction to 0.0
    frac = 0.0

    ## meaning, for mean-based criteria, specific index and fraction aren't applicable or needed

  } else if (crit == "surv.mean" || crit == "surv.prob") {
    ## retrieve the "Index" and the"sFraction" from the params object
    ## these are output depending on the input "survivalTime" parameter
    ind = params@sIndex
    frac = params@sFraction
  }

  # set basic parameter values in Fortran
  ## call a Fortran subroutine named "setUpBasics"
  ## in dtrSurv.f90 code
  ## basically just outputs the input info

  res = .Fortran(
    "setUpBasics",

    ## retrieve the number of timepoints
    t_nt = as.integer(x = .NTimes(object = params)),

    ## retrieve the difference in timepoints
    t_dt = as.double(x = .TimeDiff(object = params)),

    ## retrieves the "randomSplit" (input) parameter as a double
    ## this is then fed into fortran as the probability of random split (rs)
    t_rs = as.double(x = params@randomSplit),

    ## retrieves logical of whether extremeley randomized trees *ERT* are used (as integer)
    t_ERT = as.integer(x = params@ERT),

    ## retrieves logical of whether uniform splits are used (as integer)
    t_uniformSplit = as.integer(x = params@uniformSplit),

    ## min number of samples required at a tree node to consider a split (as integer)
    ## defined in class_IH.TreeConditions
    t_nodeSize = as.integer(x = .NodeSize(object = params)),

    ## min number of events present in node (as integer)
    ## defined in class_TreeConditions.R
    ## return the minEvent slot (integer): minimum number of events allowed in a node

    t_minEvent = as.integer(x = .MinEvent(object = params)),

    ## logical comparison of whether the "logrank" split rule is used or not
    ## if not, the truncated mean is used
    t_rule = as.integer(x = params@splitRule == 'logrank'),

    ## index (as determined by type of critical value estimator)
    t_sIndex = as.integer(x = ind),

    ## fraction (as determined by type of critical value estimator)
    t_sFraction = as.double(x = frac),

    ## A numeric object. The stratified random split coefficient (as double)
    t_stratifiedSplit = as.double(x = params@stratifiedSplit),

    ## whether or not samplingw/ replacement is used (as an integer)
    ## logical object. TRUE = sample with replacemen
    t_replace = as.integer(params@replace),
    PACKAGE = "IHsurvrf"
  )

  ## sample size verification

  # ensure that if given, sampleSize is 0 < sampleSize <= 1 and that
  # a value is provided for each decision point. If only 1 value is given,
  # it is assumed to be used for all decision points

  ## used in VerifySampleSize.R Script
  ## 'sampleSize' is numeric, numeric vector, or NULL. If NULL, set
  ##   to a default value based on ERT selection.

  ## if input is NULL & ERT is used, defaults sampleSize to 1.0 (using all data)
  ## if input is NULL & ERT is not used, uses a default fraction of 0.632 of the data used

  sampleSize <- .VerifySampleSize(sampleSize = sampleSize,
                                  ERT = params@ERT,
                                  nDP = nDP)

  ## checks if mTry is a scalar (has length 1). If so,

  # ensure that mTry is provided as a vector. At this point, there is
  # no verification of an appropriate value
  if (length(x = mTry) == 1L) {
    ## replicate this to match the number of decision points
    mTry <- rep(x = mTry, times = nDP)

    ## if mTry isn't NULL, and the input lenth if it doesn't match the number of decision points (and it's a vector), output error
  } else if (!is.null(x = mTry) && {
    length(x = mTry) != nDP
  }) {
    stop("if provided as vector, mTry must be provided for each dp",
         call. = FALSE)
  }

  ### now, we turn the data into long form so that we can operate by stage

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
    mutate(A.original = A) %>%
    arrange(subj.id, stage)

  ## we want to define the cumulative percentage of events that have occurred (delta = 0)
  ## when we go back a time window, we want to include in the upper strata the stages for patients who have cumulative time higher


  # Initialize starting point
  starting_thresh <- tau - windowsize

  # Initialize strata columns
  ## we want to do this automatically, and initialize a column for each strata
# Loop through the number of strata and add columns to long_data
for (i in 1:nstrata) {
  column_name <- paste0("strata", i)
  long_data[[column_name]] <- 0

}



  ############ semi-manual assigning of strata...
  ############ if strata = 1, we assign all membership in strata 1
  ############ if strata = 2, we assign the original way with a 70%/30% split
  ############ if strata = 5, we assign each strata with an equal split
  ############

  ###### run this if strata == 1
  if (nstrata == 1) {
    long_data$strata1 <- 1

  } else if (nstrata == 2) {

  ###### run this if strata == 2
  # Iterate over the data using a for loop
  for (i in seq(from = starting_thresh, to = 0, by = -windowsize)) {
    # Select data where cumulative time is greater than or equal to the current threshold
    selected_data <- long_data %>% filter(cumulative.time * tau >= i)

    # Calculate the percentage of events in current cutoff divided by total number of events
    event_percentage <- sum(selected_data$delta == 0, na.rm = TRUE) / sum(long_data$delta == 0, na.rm = TRUE)

    # Check if the event percentage meets the threshold
    ## NOTE: that the dimensions of these will not match the original data just because the original data contains rows of NAs
    ## this is because we simulate every patient to have 10 stages but if they have event/censored earlier they just have NA
    if (event_percentage >= 0.3) {
      # Assign strata membership
      long_data$strata1[ (long_data$cumulative.time * tau) >= i] <- 1
      long_data$strata2[ (long_data$cumulative.time * tau) < i] <- 1

      ## subset the data
      strata1 <- selected_data
      strata2 <- long_data %>% filter(cumulative.time*tau < i)

      selected_thresh <- i

      ### print the cutoff time for the strata
      message("Cumulative cutoff:", selected_thresh, "; longest length:", tau)


      break


    }
  }
} else if (nstrata == 5) {

  ###### run this if strata == 5

  nstrata <- 5  # Number of strata
  thresholds <- c(0.65, 0.07, 0.07, 0.07)  # Event cutoffs for the first n-1 strata

  # Initialize variables
  remaining_data <- long_data
  previous_thresh <- starting_thresh

  # Initialize a vector to store the values of i
  selected_thresh <- c()

  # Loop through the first n-1 strata
  for (k in 1:(nstrata - 1)) {
    current_threshold <- thresholds[k]

    # Loop over time windows starting from the previous strata's ending time
    for (i in seq(from = previous_thresh, to = 0, by = -windowsize)) {
      selected_data <- remaining_data %>% filter(cumulative.time * tau >= i)

      # Calculate the percentage of events in the current window
      event_percentage <- sum(selected_data$delta == 0, na.rm = TRUE) / sum(long_data$delta == 0, na.rm = TRUE)

      # Check if the event percentage meets or exceeds the threshold for the current strata
      if (event_percentage >= current_threshold) {
        strata_col <- paste0("strata", k)

        # Assign to the current strata only within the remaining data
        long_data[[strata_col]] <- ifelse(long_data[[strata_col]] == 0 & (long_data$cumulative.time * tau) >= i, 1, long_data[[strata_col]])

        # Subset the remaining data for the next strata
        remaining_data <- remaining_data %>% filter(cumulative.time * tau < i)

        # Update the previous threshold for the next strata
        previous_thresh <- i

        # Store the used value of i
        selected_thresh <- c(selected_thresh, i)

        # Print the cutoff time for the strata
        message("Cumulative cutoff for strata ", k, ":", i, "; event percentage:", event_percentage)

        # Exit the loop once the current strata is assigned
        break
      }
    }
  }

  # Ensure mutually exclusive membership based on the smallest strata
  for (i in 1:nrow(long_data)) {
    # Find the smallest strata with a value of 1
    assigned_strata <- which(long_data[i, paste0("strata", 1:nstrata)] == 1)

    if (length(assigned_strata) > 0) {
      # Keep only the smallest strata, set the rest to 0
      smallest_strata <- min(assigned_strata)
      long_data[i, paste0("strata", 1:nstrata)] <- 0
      long_data[i, paste0("strata", smallest_strata)] <- 1
    }
  }

  # Assign the remaining data to the nth strata
  # Assign the remaining data to the last strata using an ifelse statement
  long_data$strata5 <- ifelse(
    long_data$strata1 == 0 & long_data$strata2 == 0 & long_data$strata3 == 0 & long_data$strata4 == 0,
    1,
    0
  )

}



  #############
  #############

## commented out: plotting visualizations
#####
##### plotting to see how many stages a patient has
#####

#  #### see how many stages patients have
#strata1_stages <- long_data %>%
#  filter(strata1 == 1 & !is.na(T)) %>%
#  group_by(subj.id) %>%
#  summarise(stage_count = n()) %>%
#  ungroup() %>%
#  count(stage_count)
#
#strata2_stages <- long_data %>%
#  filter(strata2 == 1 & !is.na(T)) %>%
#  group_by(subj.id) %>%
#  summarise(stage_count = n()) %>%
#  ungroup() %>%
#  count(stage_count)
#
## Create the bar plot
#ggplot(strata1_stages, aes(x = stage_count, y = n)) +
#  geom_bar(stat = "identity") +
#  labs(title = "Number of Stages Patients Have in Strata 1,
#       with 10% split in events",
#       x = "Number of Stages",
#       y = "Number of Patients") +
#  theme_minimal()

#####################################
#####################################
#####################################


  # Step 1: Dynamically get column names starting with "strata"
  strata_columns <- names(long_data)[grepl("^strata", names(long_data))]


  # Step 2: Convert back to wide format
  data <- long_data %>%
    pivot_wider(
      id_cols = subj.id,
      names_from = stage,
      values_from = c(variables,  all_of(strata_columns)),
      names_sep = "_"
    )

  # Get the column name as a string
  respname <- as.character(attr(terms(models), "variables")[[2]][[2]])

  ### now, for strata1, we want to pool the observations together, treating each pts stage separately
  ### fit a forest using this

  #### NOTE: no matter how many strata we have, we will always use this step initially
  ## if we only have 1 strata, all the pts are captured in this
  ## if we have 2 strata, this is our first, and our second uses strata2
  ## if we have 5 strata, this is our first, our second will use strata2, third strata3, etc...

  rm(variables)
  rm(starting_thresh)
  gc()

  s1.strata1 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata with complete cases (don't include if they've failed already)
    data = long_data %>%
      filter(strata1 == 1 & !is.na(.data[[respname]])),
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
    mTry = mTry[[1]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## we are in the first step of pooling patient data after running HC code
    pool1 = F,

    ## we're not inputting any survival probabilities
    appendstep1 = F
  )


  ## now, for strata 1, we want to track the optimal and eligibility
  ## first initialize a column name
  long_data$A.s1.strata1 <- NA

  eligibility_s1.strata1 <- s1.strata1@eligibility





  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE
  ### the eligibility is only for patients in strata 1 AND not NA, so we have to filter for that too
  long_data$A.s1.strata1[long_data$strata1 == 1 & !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata1 == 1)] <-
    s1.strata1@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data$strata1 == 1 & !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata1 == 1)] <- s1.strata1@optimal@optimalTx


  # Initialize the probability matrix that's output
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata1) )
  shiftedprobfinal[,eligibility_s1.strata1] <-t(s1.strata1@optimal@optimalY)

## note: we don't know the final stage's probability since each patient's final stage may be different




  ### ------------ SCRATCH FOLLOWING A PATIENT AND TRACKING SURVIVAL CURVE ---------------- ###

#  conv_1_curve <- res.strata1.1@prev_probs[, 1]
#  conv_2_curve <- res.strata1.1@prev_probs[, 1]
#  conv_3_curve <- res.strata1.1@prev_probs[, 1]
#  conv_4_curve <- res.strata1.1@prev_probs[, 1]
#  conv_5_curve <- res.strata1.1@prev_probs[, 1]
#  conv_6_curve <- res.strata1.1@prev_probs[, 1]
#  conv_7_curve <- res.strata1.1@prev_probs[, 1]
#  conv_8_curve <- res.strata1.1@prev_probs[, 1]
#  conv_9_curve <- res.strata1.1@prev_probs[, 1]
#  conv_23_curve <- res.strata1.1@prev_probs[, 1]
#
#  pool1_curve <- cbind(params@timePoints, shiftedprobfinal[, 1], conv_1_curve, conv_2_curve, conv_3_curve, conv_4_curve, conv_5_curve,
#                       conv_6_curve, conv_7_curve, conv_8_curve, conv_9_curve, conv_10_curve,
#                       conv_11_curve, conv_12_curve, conv_13_curve, conv_14_curve, conv_15_curve,
#                       conv_16_curve, conv_17_curve, conv_19_curve, conv_20_curve, conv_21_curve,
#                       conv_22_curve, conv_23_curve)
#

  ### ------------------------------------------------------------------------------------- ###




  ## the result after appyling the area function to each column of the matrix of survival probabilities
## since each patient has a different number of visits in the strata, we look at the area of all stages included
  areas <- apply(shiftedprobfinal, 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })

  ## after backwards recursion is complete, calculate the estimated value from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  valueTrain <- .meanValue(object = s1.strata1)


  ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
  }


  # store values in call structure for returned object

  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  # Initialize a flag to indicate whether to continue iterations
  continue_iterations <- TRUE

  conv_iterations <- 1

  ## now we want to create a matrix for areas where each row is one iteration of the forest
  ## but, we want to zero out any values where the T was too close to 0 so we can compare survival curves
  area_mat <- areas
  ## we add this part to areas if action > 2
  #[-which(long_data$strata1 == 1 & long_data$T < 1e-8)]

  # Initialize vectors to store avg_diff and res@valueTrain values
  avg_diff_values <- c()
  valueTrain_values <- c(valueTrain)

  ## commented out: visualizations
  ##############################
  ############################## VISUALIZATION
  ##############################

#  xaxis <- params@timePoints
#  y1 <- t(s1.strata1@optimal@optimalY)[,33]
#  plot(xaxis, y1)
#
  ##############################
  ##############################
  ##############################

  ### now, we need to work on convergence
  ### we start at stage nDP (10), see if any patients are in strata1
  #### if so, we predict a stub using s1.strata1 forest
  ### go back to stage nDP - 1 (9), see if any patients are in strata1
  ### if they are, and they were eligible in previous stage, append to get double stub
  ### if they are, and this is their first stage, predict using s1.strata1 forest
  #### keep going back until stage 1
  ### once we make predictions for all these stages, use these predicted stubs/double stubs to fit another forest s2.strata1
  ### we repeat until forest convergence

  # Construct an object of class DTRSurv for the output
  ## this needs to be initialized outside the loop
  res.strata1.1 <- new(
    Class = "DTRSurv",

    #####
    ###### delete these?
    ######
    "stageResults" = list(seq(1:nDP)),

    #####
    ###### delete these?
    ######

    "IHstageResults" = list(),
    "FinalForest" = s1.strata1,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = area_mat,
    "n_it" = conv_iterations,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list(),
    "long_data" = long_data,

    ## empty matrix since we don't use previous probabilities
    "prev_probs" = matrix(nrow = 2, ncol = 2)
  )

  ## free up some memory

  rm(s1.strata1)
  gc()

  max_iter <- long_data %>%
    filter(strata1 == 1 & !is.na(.data[[respname]])) %>%
    group_by(subj.id) %>%           # Group by unique patient ID
    summarise(num_stages = n_distinct(stage)) %>%  # Count distinct stages for each patient
    summarise(max_stages = max(num_stages)) %>%    # Find the maximum number of stages
    pull(max_stages)

  #########
  ## now, for each patient, we find their max # of stages to know how many iterations we need to run before freezing their values:
  #### count backwards where the final stage is 1, ...their last stage, put NA in all the other entries (# pts * nDP) size

  ## create an eligibility matrix for whether we are looking at a patient's conv_iteration stage; whether the current iteration counter matches the counter above ^^

  while(continue_iterations && conv_iterations <= max_iter){

    message("Convergence Re-fitting Iteration:", conv_iterations)

    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv(data = data,
                                 prev.iteration = res.strata1.1, nDP = nDP, params = params, nTimes = nTimes,
                                 models = models, mTry = mTry, strata = 1, long_data = res.strata1.1@long_data,
                                 prev_probs = res.strata1.1@prev_probs)

    ######## for the pt that was zeroed out, we want to get rid of their integral from the original res.strata1

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <- rbind(res.strata1.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata1.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata1.1@prev_probs <- convergence_res@prev_probs

    ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
    last_two_rows_diff <- (abs(diff(area_mat[ (nrow(area_mat)-1):nrow(area_mat), ]))/area_mat[ (nrow(area_mat)-1), ]) * 100

    avg_diff <- mean(last_two_rows_diff, na.rm = T)

    # Store avg_diff value to track each iteration
    avg_diff_values <- c(avg_diff_values, avg_diff)

    print(avg_diff)

    ## we update the final forest used with the most recent forest estimated in the convergence step
    res.strata1.1@FinalForest <- convergence_res@FinalForest

    ## we update the value with the most recent estimated value in the convergence step
    res.strata1.1@value <- convergence_res@value

    ## we also update the long_data
    res.strata1.1@long_data <- convergence_res@long_data

    # Store res@valueTrain value to track each iteration
    valueTrain_values <- c(valueTrain_values, res.strata1.1@value)

    ## track these values in the forest output
    res.strata1.1@n_it <- conv_iterations
    res.strata1.1@avgKM_diff <- as.matrix(avg_diff_values)
    res.strata1.1@valueTrain_list <- valueTrain_values

    ## Increment the iteration counter
    conv_iterations <- conv_iterations + 1

    ## Check if the maximum number of iterations is reached
    if(conv_iterations > max_iter) {
      continue_iterations <- FALSE  # Stop the loop when max_iter is reached
    }

    # Free memory by removing convergence_res and triggering garbage collection after each iteration
    rm(convergence_res)
    gc()
  }



## we only run these steps if nstrata > 1
  ## we want to dynamically create this forest + the convergence loop depending on the number of strata
  #### we loop through each of the strata from 2:nstrata
  if (nstrata == 2) {

    # Get the variable name as a character
    resp_name <- deparse(attr(terms(models), "variables")[[2]][[2]])

  ## now, for strata2, we want to pool the observations together, treating each pts stage separately
  s1.strata2 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata that's complete
    data = long_data %>% filter(!!sym(paste0("strata",2)) == 1 & !is.na(!!sym(resp_name))),
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
    mTry = mTry[[1]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## we are in the first step of pooling patient data after running HC code
    pool1 = F,

    ## this is true if we are inputting survival probabilities
    appendstep1 = F
  )


  ## now, for strata 2, we want to track the optimal and eligibility
  ## first initialize a column name
  long_data$A.s1.strata2 <- NA

  eligibility_s1.strata2 <- s1.strata2@eligibility

  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE


  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE
  ### the eligibility is only for patients in strata 1 AND not NA, so we have to filter for that too
  long_data$A.s1.strata2[long_data[[paste0("strata", 2)]] == 1 & !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata2 == 1)] <- s1.strata2@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data[[paste0("strata", 2)]] == 1 & !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata2 == 1)] <- s1.strata2@optimal@optimalTx



  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata2) )
  shiftedprobfinal[,eligibility_s1.strata2] <-t(s1.strata2@optimal@optimalY)



  ## the result after appyling the area function to each column of the matrix of survival probabilities
  ##### NOTE: the issue is that each patient has a different number of visits belonging to this strata
  ##### meaning, for visits in strata 2, we must input 0 probability
  areas <- apply(shiftedprobfinal, 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })

  ## now we want to create a matrix for areas where each row is one iteration of the forest
  area_mat <- areas


  # Construct an object of class DTRSurv for the output
  ## this needs to be initialized outside the loop
  res.strata2.1 <- new(
    Class = "DTRSurv",

    #####
    ###### delete these?
    ######
    "stageResults" = list(seq(1:nDP)),

    #####
    ###### delete these?
    ######

    "IHstageResults" = list(),
    "FinalForest" = s1.strata2,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = area_mat,
    "n_it" = conv_iterations,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list(),
    "long_data" = long_data,

    ## NOTE: for the previous probs, we want to access res.strata(i - 1).1
    ############# THIS IS ANOTHER CHANGE TO MAKE


    "prev_probs" = res.strata1.1@prev_probs
  )


  ## after backwards recursion is complete, calculate the estimated value from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  valueTrain <- .meanValue(object = s1.strata2)


  ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
  }


  # store values in call structure for returned object

  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  # Initialize a flag to indicate whether to continue iterations
  continue_iterations <- TRUE

  conv_iterations <- 1

  # Initialize vectors to store avg_diff and res@valueTrain values
  avg_diff_values <- c()
  valueTrain_values <- c(valueTrain)

  ## free up memory
  rm(s1.strata2)
  gc()


  max_iter <- long_data %>%
    filter(strata2 == 1 & !is.na(.data[[respname]])) %>%
    group_by(subj.id) %>%           # Group by unique patient ID
    summarise(num_stages = n_distinct(stage)) %>%  # Count distinct stages for each patient
    summarise(max_stages = max(num_stages)) %>%    # Find the maximum number of stages
    pull(max_stages)




  while(continue_iterations && conv_iterations <= max_iter){

    message("Convergence Re-fitting Iteration:", conv_iterations)

    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv_otherstrata(data = data,
                                             prev.iteration = res.strata2.1, nDP = nDP, params = params, nTimes = nTimes,
                                             models = models, mTry = mTry, strata = 2,
                                             # use the most recent long_data
                                             long_data = res.strata2.1@long_data,
                                             prev_probs = res.strata2.1@prev_probs)
    ######## for the pt that was zeroed out, we want to get rid of their integral from the original res.strata1

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <- rbind(res.strata2.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata2.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata2.1@prev_probs <- convergence_res@prev_probs

    ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
    last_two_rows_diff <- (abs(diff(area_mat[ (nrow(area_mat)-1):nrow(area_mat), ]))/area_mat[ (nrow(area_mat)-1), ]) * 100

    avg_diff <- mean(last_two_rows_diff, na.rm = T)

    # Store avg_diff value to track each iteration
    avg_diff_values <- c(avg_diff_values, avg_diff)

    print(avg_diff)

    ## we update the final forest used with the most recent forest estimated in the convergence step
    res.strata2.1@FinalForest <- convergence_res@FinalForest

    ## we update the value with the most recent estimated value in the convergence step
    res.strata2.1@value <- convergence_res@value

    ## we also update the long_data
    res.strata2.1@long_data <- convergence_res@long_data

    # Store res@valueTrain value to track each iteration
    valueTrain_values <- c(valueTrain_values, res.strata2.1@value)

    ## track these values in the forest output
    res.strata2.1@n_it <- conv_iterations
    res.strata2.1@avgKM_diff <- as.matrix(avg_diff_values)
    res.strata2.1@valueTrain_list <- valueTrain_values

    ## Increment the iteration counter
    conv_iterations <- conv_iterations + 1

    ## Check if the maximum number of iterations is reached
    if(conv_iterations > max_iter) {
      continue_iterations <- FALSE  # Stop the loop when max_iter is reached
    }

    # Free memory by removing convergence_res and triggering garbage collection after each iteration
    rm(convergence_res)
    gc()
  }


  ## we create this one if our nstrata == 2
  ###### NOTE: need to create a different result if nstrata == 1 and if nstrata == 5
  res <- new(
    Class = "DTRSurvRes",

    ## strata 1 forest

    "Forest1" = res.strata1.1,

    ## strata 2 forest
    "Forest2" = res.strata2.1,

    "Forest3" = NA,
    "Forest4" = NA,
    "Forest5" = NA,

    "call" = cl,
    "params" = params,
    "long_data" = res.strata2.1@long_data,
    "prev_probs" = res.strata2.1@prev_probs,
    "n_it" = NA,
    "avgKM_diff" = matrix(0),
    "cutoff" = selected_thresh
  )



   } else if (nstrata == 5) {

       ## now, for strata2, we want to pool the observations together, treating each pts stage separately
       s1.strata2 <- .dtrSurvStep(
         ## use the model for pooled data
         model = models,
         ## only include the data in the first strata that's complete
         data = long_data %>% filter(!!sym(paste0("strata",2)) == 1 & !is.na(as.character(attr(terms(models), "variables")[[2]][[2]]))),
         priorStep = NULL,
         params = params,
         txName = "A",
         ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
         mTry = mTry[[1]],
         ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
         sampleSize = 1,
         ## we are in the first step of pooling patient data after running HC code
         pool1 = F,

         ## this is true if we are inputting survival probabilities
         appendstep1 = F
         #    inputpr = input.strata2[, colSums(input.strata2) != 0]
       )

       ## now, for strata 2, we want to track the optimal and eligibility
       ## first initialize a column name
       long_data$A.s1.strata2 <- NA

       eligibility_s1.strata2 <- s1.strata2@eligibility

       # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE


       # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE
       ### the eligibility is only for patients in strata 1 AND not NA, so we have to filter for that too
       long_data$A.s1.strata2[long_data[[paste0("strata", 2)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata2 == 1)] <- s1.strata2@optimal@optimalTx
       ## also update the "A" column to be used for the convergence aspect
       long_data$A[long_data[[paste0("strata", 2)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata2 == 1)] <- s1.strata2@optimal@optimalTx



       # Initialize the shifted probability matrix
       ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
       ## we overwrite the eligible patients
       shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata2) )
       shiftedprobfinal[,eligibility_s1.strata2] <-t(s1.strata2@optimal@optimalY)



       ## the result after appyling the area function to each column of the matrix of survival probabilities
       ##### NOTE: the issue is that each patient has a different number of visits belonging to this strata
       ##### meaning, for visits in strata 2, we must input 0 probability
       areas <- apply(shiftedprobfinal, 2, function(surv_prob_col) {
         area_under_curve(surv_prob_col, params@timePoints)
       })

       ## now we want to create a matrix for areas where each row is one iteration of the forest
       area_mat <- areas


       # Construct an object of class DTRSurv for the output
       ## this needs to be initialized outside the loop
       res.strata2.1 <- new(
         Class = "DTRSurv",

         #####
         ###### delete these?
         ######
         "stageResults" = list(seq(1:nDP)),

         #####
         ###### delete these?
         ######

         "IHstageResults" = list(),
         "FinalForest" = s1.strata2,
         "value" = valueTrain,
         "call" = cl,
         "params" = params,
         "integral_KM" = area_mat,
         "n_it" = conv_iterations,
         "avgKM_diff" = matrix(nrow = 2, ncol = 2),
         "valueTrain_list" = list(),
         "long_data" = long_data,

         ## NOTE: for the previous probs, we want to access res.strata(i - 1).1
         ############# THIS IS ANOTHER CHANGE TO MAKE


         "prev_probs" = res.strata1.1@prev_probs
       )

       ## after backwards recursion is complete, calculate the estimated value from the first stage
       ## calculates mean values of expected survival times and survival probabilities
       ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
       ## .meanValue() function defined in class_DTRSurvStep.R

       valueTrain <- .meanValue(object = s1.strata2)


       ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

       message("Estimated Value:", appendLF = FALSE)
       for (i in 1L:length(valueTrain)) {
         message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
       }


       # store values in call structure for returned object

       ## captures current function call, including function name and all arguments passed to it
       cl <- match.call()

       ## ensures name of called function is set to "dtrSurv"
       cl[[1L]] <- as.name("IHdtrSurv")

       # Initialize a flag to indicate whether to continue iterations
       continue_iterations <- TRUE

       conv_iterations <- 1

       # Initialize vectors to store avg_diff and res@valueTrain values
       avg_diff_values <- c()
       valueTrain_values <- c(valueTrain)


       while(continue_iterations){

         message("Convergence Re-fitting Iteration:", conv_iterations)


         ####### Here, we will construct a new iteration of forest training to check for convergence
         convergence_res <- IHdtrConv_otherstrata(data = data,
                                                  prev.iteration = res.strata2.1, nDP = nDP, params = params, nTimes = nTimes,
                                                  models = models, mTry = mTry, strata = 2,
                                                  # use the most recent long_data
                                                  long_data = res.strata2.1@long_data,
                                                  prev_probs = res.strata2.1@prev_probs)

         ## now we want to create a matrix for areas where each row is one iteration of the forest
         area_mat <- rbind(res.strata2.1@integral_KM, convergence_res@integral_KM)

         # Update res@integral_KM with area_mat
         res.strata2.1@integral_KM <- area_mat

         # update the optimal output probabilities
         res.strata2.1@prev_probs <- convergence_res@prev_probs


         ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
         #### meaning, for the same patient, if the difference between their estimated survival curves is large, we go through another iteration of training
         ## we update res@FinalForest with the forest in convergence_res@FinalForest
         ## we update the matrix of areas

         # Check the absolute value of the percent difference
         last_two_rows_diff <- ( abs(diff(area_mat[ (nrow(area_mat)-1):nrow(area_mat), ]))/area_mat[ (nrow(area_mat)-1), ])*100

         avg_diff <- mean(last_two_rows_diff)

         # Store avg_diff value to track each iteration
         avg_diff_values <- c(avg_diff_values, avg_diff)

         ## we update the final forest used with the most recent forest estimated in the convergence step
         # Update res@FinalForest with the forest in convergence_res@FinalForest
         res.strata2.1@FinalForest <- convergence_res@FinalForest

         ## we update the value with the most recent estimated value in the convergence steo
         res.strata2.1@value <- convergence_res@value

         ## we also update the long_data
         res.strata2.1@long_data <- convergence_res@long_data

         # Store res@valueTrain value to track each iteration
         valueTrain_values <- c(valueTrain_values, res.strata2.1@value)

         ## track these values in the forest output
         res.strata2.1@n_it <- conv_iterations
         res.strata2.1@avgKM_diff <- as.matrix(avg_diff_values)
         res.strata2.1@valueTrain_list <- valueTrain_values



         ### if the absolute change starts fluctuating and the last 10 iteration don't change by more than 0.5, we stop
         ## OR, if the forest itself has converging curves
         change_last <- mean(abs(diff(tail(as.matrix(avg_diff_values), 3))))


         ## wait until the absolute change in survival curves is less than something or the average change gets small

         if (avg_diff > 0.1) {
           if (is.na(change_last) || is.nan(change_last) || abs(change_last) > 1) {
             # If change_last5 is NaN or abs(change_last5) > 0.5, continue the loop
             continue_iterations <- TRUE

             cat("Difference in survival curve:", avg_diff, "\n")
             cat("Last 3 stages change:", change_last, "\n")

             # Increment the iteration counter
             conv_iterations <- conv_iterations + 1
           } else {
             # If avg_diff > 0.1 but abs(change_last5) <= 0.5, stop the loop
             continue_iterations <- FALSE

             cat("Difference in survival curve:", avg_diff, "\n")
             cat("Last 5 stages change:", change_last, "\n")

             # Update res@integral_KM with area_mat
             res.strata2.1@integral_KM <- area_mat

             # Update res@FinalForest with the forest in convergence_res@FinalForest
             res.strata2.1@FinalForest <- convergence_res@FinalForest

             # Update the value with the most recent estimated value in the convergence step
             res.strata2.1@value <- convergence_res@value

             # Update the long_data
             res.strata2.1@long_data <- convergence_res@long_data

             # Update the optimal output probabilities
             res.strata2.1@prev_probs <- convergence_res@prev_probs

             # Track these values in the forest output
             res.strata2.1@n_it <- conv_iterations
             res.strata2.1@avgKM_diff <- as.matrix(avg_diff_values)
             res.strata2.1@valueTrain_list <- valueTrain_values
           }
         } else {
           # If avg_diff <= 0.1, stop the loop
           continue_iterations <- FALSE

           cat("Difference in survival curve:", avg_diff, "\n")
           cat("Last 5 stages change:", change_last, "\n")

           # Update res@integral_KM with area_mat
           res.strata2.1@integral_KM <- area_mat

           # Update res@FinalForest with the forest in convergence_res@FinalForest
           res.strata2.1@FinalForest <- convergence_res@FinalForest

           # Update the value with the most recent estimated value in the convergence step
           res.strata2.1@value <- convergence_res@value

           # Update the long_data
           res.strata2.1@long_data <- convergence_res@long_data

           # Update the optimal output probabilities
           res.strata2.1@prev_probs <- convergence_res@prev_probs

           # Track these values in the forest output
           res.strata2.1@n_it <- conv_iterations
           res.strata2.1@avgKM_diff <- as.matrix(avg_diff_values)
           res.strata2.1@valueTrain_list <- valueTrain_values
         }
       }

  ########## for strata 3-5, we only run this portion if nstrata == 5
  ########## ADD THIS CONDITIONAL LOOP IN

  ## now, for strata3, we want to pool the observations together, treating each pts stage separately
  s1.strata3 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata that's complete
    data = long_data %>% filter(!!sym(paste0("strata",3)) == 1 & !is.na(as.character(attr(terms(models), "variables")[[2]][[2]]))),
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
    mTry = mTry[[1]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## we are in the first step of pooling patient data after running HC code
    pool1 = F,

    ## this is true if we are inputting survival probabilities
    appendstep1 = F
    #    inputpr = input.strata2[, colSums(input.strata2) != 0]
  )

  ## now, for strata 2, we want to track the optimal and eligibility
  ## first initialize a column name
  long_data$A.s1.strata3 <- NA

  eligibility_s1.strata3 <- s1.strata3@eligibility

  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE


  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE
  ### the eligibility is only for patients in strata 1 AND not NA, so we have to filter for that too
  long_data$A.s1.strata3[long_data[[paste0("strata", 3)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata3 == 1)] <- s1.strata3@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data[[paste0("strata", 3)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata3 == 1)] <- s1.strata3@optimal@optimalTx



  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata3) )
  shiftedprobfinal[,eligibility_s1.strata3] <-t(s1.strata3@optimal@optimalY)



  ## the result after appyling the area function to each column of the matrix of survival probabilities
  ##### NOTE: the issue is that each patient has a different number of visits belonging to this strata
  ##### meaning, for visits in strata 2, we must input 0 probability
  areas <- apply(shiftedprobfinal, 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })

  ## now we want to create a matrix for areas where each row is one iteration of the forest
  area_mat <- areas


  # Construct an object of class DTRSurv for the output
  ## this needs to be initialized outside the loop
  res.strata3.1 <- new(
    Class = "DTRSurv",

    #####
    ###### delete these?
    ######
    "stageResults" = list(seq(1:nDP)),

    #####
    ###### delete these?
    ######

    "IHstageResults" = list(),
    "FinalForest" = s1.strata3,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = area_mat,
    "n_it" = conv_iterations,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list(),
    "long_data" = long_data,

    ## NOTE: for the previous probs, we want to access res.strata(i - 1).1
    ############# THIS IS ANOTHER CHANGE TO MAKE


    "prev_probs" = res.strata2.1@prev_probs
  )

  ## after backwards recursion is complete, calculate the estimated value from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  valueTrain <- .meanValue(object = s1.strata3)


  ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
  }


  # store values in call structure for returned object

  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  # Initialize a flag to indicate whether to continue iterations
  continue_iterations <- TRUE

  conv_iterations <- 1

  # Initialize vectors to store avg_diff and res@valueTrain values
  avg_diff_values <- c()
  valueTrain_values <- c(valueTrain)


  while(continue_iterations){

    message("Convergence Re-fitting Iteration:", conv_iterations)


    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv_otherstrata(data = data,
                                             prev.iteration = res.strata3.1, nDP = nDP, params = params, nTimes = nTimes,
                                             models = models, mTry = mTry, strata = 3,
                                             # use the most recent long_data
                                             long_data = res.strata3.1@long_data,
                                             prev_probs = res.strata3.1@prev_probs)

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <- rbind(res.strata3.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata3.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata3.1@prev_probs <- convergence_res@prev_probs


    ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
    #### meaning, for the same patient, if the difference between their estimated survival curves is large, we go through another iteration of training
    ## we update res@FinalForest with the forest in convergence_res@FinalForest
    ## we update the matrix of areas

    # Check the absolute value of the percent difference
    last_two_rows_diff <- ( abs(diff(area_mat[ (nrow(area_mat)-1):nrow(area_mat), ]))/area_mat[ (nrow(area_mat)-1), ])*100

    avg_diff <- mean(last_two_rows_diff)

    # Store avg_diff value to track each iteration
    avg_diff_values <- c(avg_diff_values, avg_diff)

    ## we update the final forest used with the most recent forest estimated in the convergence step
    # Update res@FinalForest with the forest in convergence_res@FinalForest
    res.strata3.1@FinalForest <- convergence_res@FinalForest

    ## we update the value with the most recent estimated value in the convergence steo
    res.strata3.1@value <- convergence_res@value

    ## we also update the long_data
    res.strata3.1@long_data <- convergence_res@long_data

    # Store res@valueTrain value to track each iteration
    valueTrain_values <- c(valueTrain_values, res.strata3.1@value)

    ## track these values in the forest output
    res.strata3.1@n_it <- conv_iterations
    res.strata3.1@avgKM_diff <- as.matrix(avg_diff_values)
    res.strata3.1@valueTrain_list <- valueTrain_values



    ### if the absolute change starts fluctuating and the last 10 iteration don't change by more than 0.5, we stop
    ## OR, if the forest itself has converging curves
    change_last <- mean(abs(diff(tail(as.matrix(avg_diff_values), 3))))


    ## wait until the absolute change in survival curves is less than something or the average change gets small

    if (avg_diff > 0.1) {
      if (is.na(change_last) || is.nan(change_last) || abs(change_last) > 1) {
        # If change_last5 is NaN or abs(change_last5) > 0.5, continue the loop
        continue_iterations <- TRUE

        cat("Difference in survival curve:", avg_diff, "\n")
        cat("Last 3 stages change:", change_last, "\n")

        # Increment the iteration counter
        conv_iterations <- conv_iterations + 1
      } else {
        # If avg_diff > 0.1 but abs(change_last5) <= 0.5, stop the loop
        continue_iterations <- FALSE

        cat("Difference in survival curve:", avg_diff, "\n")
        cat("Last 5 stages change:", change_last, "\n")

        # Update res@integral_KM with area_mat
        res.strata3.1@integral_KM <- area_mat

        # Update res@FinalForest with the forest in convergence_res@FinalForest
        res.strata3.1@FinalForest <- convergence_res@FinalForest

        # Update the value with the most recent estimated value in the convergence step
        res.strata3.1@value <- convergence_res@value

        # Update the long_data
        res.strata3.1@long_data <- convergence_res@long_data

        # Update the optimal output probabilities
        res.strata3.1@prev_probs <- convergence_res@prev_probs

        # Track these values in the forest output
        res.strata3.1@n_it <- conv_iterations
        res.strata3.1@avgKM_diff <- as.matrix(avg_diff_values)
        res.strata3.1@valueTrain_list <- valueTrain_values
      }
    } else {
      # If avg_diff <= 0.1, stop the loop
      continue_iterations <- FALSE

      cat("Difference in survival curve:", avg_diff, "\n")
      cat("Last 5 stages change:", change_last, "\n")

      # Update res@integral_KM with area_mat
      res.strata3.1@integral_KM <- area_mat

      # Update res@FinalForest with the forest in convergence_res@FinalForest
      res.strata3.1@FinalForest <- convergence_res@FinalForest

      # Update the value with the most recent estimated value in the convergence step
      res.strata3.1@value <- convergence_res@value

      # Update the long_data
      res.strata3.1@long_data <- convergence_res@long_data

      # Update the optimal output probabilities
      res.strata3.1@prev_probs <- convergence_res@prev_probs

      # Track these values in the forest output
      res.strata3.1@n_it <- conv_iterations
      res.strata3.1@avgKM_diff <- as.matrix(avg_diff_values)
      res.strata3.1@valueTrain_list <- valueTrain_values
    }
  }


  ### now, we calculate info for the 4th strata

  ## now, for strata4, we want to pool the observations together, treating each pts stage separately
  s1.strata4 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata that's complete
    data = long_data %>% filter(!!sym(paste0("strata",4)) == 1 & !is.na(as.character(attr(terms(models), "variables")[[2]][[2]]))),
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
    mTry = mTry[[1]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## we are in the first step of pooling patient data after running HC code
    pool1 = F,

    ## this is true if we are inputting survival probabilities
    appendstep1 = F
    #    inputpr = input.strata2[, colSums(input.strata2) != 0]
  )

  ## now, for strata 2, we want to track the optimal and eligibility
  ## first initialize a column name
  long_data$A.s1.strata4 <- NA

  eligibility_s1.strata4 <- s1.strata4@eligibility

  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE


  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE
  ### the eligibility is only for patients in strata 1 AND not NA, so we have to filter for that too
  long_data$A.s1.strata4[long_data[[paste0("strata", 4)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata4 == 1)] <- s1.strata4@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data[[paste0("strata", 4)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata4 == 1)] <- s1.strata4@optimal@optimalTx



  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata4) )
  shiftedprobfinal[,eligibility_s1.strata4] <-t(s1.strata4@optimal@optimalY)



  ## the result after appyling the area function to each column of the matrix of survival probabilities
  ##### NOTE: the issue is that each patient has a different number of visits belonging to this strata
  ##### meaning, for visits in strata 2, we must input 0 probability
  areas <- apply(shiftedprobfinal, 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })

  ## now we want to create a matrix for areas where each row is one iteration of the forest
  area_mat <- areas


  # Construct an object of class DTRSurv for the output
  ## this needs to be initialized outside the loop
  res.strata4.1 <- new(
    Class = "DTRSurv",

    #####
    ###### delete these?
    ######
    "stageResults" = list(seq(1:nDP)),

    #####
    ###### delete these?
    ######

    "IHstageResults" = list(),
    "FinalForest" = s1.strata4,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = area_mat,
    "n_it" = conv_iterations,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list(),
    "long_data" = long_data,

    ## NOTE: for the previous probs, we want to access res.strata(i - 1).1
    ############# THIS IS ANOTHER CHANGE TO MAKE


    "prev_probs" = res.strata3.1@prev_probs
  )

  ## after backwards recursion is complete, calculate the estimated value from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  valueTrain <- .meanValue(object = s1.strata4)


  ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
  }


  # store values in call structure for returned object

  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  # Initialize a flag to indicate whether to continue iterations
  continue_iterations <- TRUE

  conv_iterations <- 1

  # Initialize vectors to store avg_diff and res@valueTrain values
  avg_diff_values <- c()
  valueTrain_values <- c(valueTrain)


  while(continue_iterations){

    message("Convergence Re-fitting Iteration:", conv_iterations)


    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv_otherstrata(data = data,
                                             prev.iteration = res.strata4.1, nDP = nDP, params = params, nTimes = nTimes,
                                             models = models, mTry = mTry, strata = 4,
                                             # use the most recent long_data
                                             long_data = res.strata4.1@long_data,
                                             prev_probs = res.strata4.1@prev_probs)

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <- rbind(res.strata4.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata4.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata4.1@prev_probs <- convergence_res@prev_probs


    ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
    #### meaning, for the same patient, if the difference between their estimated survival curves is large, we go through another iteration of training
    ## we update res@FinalForest with the forest in convergence_res@FinalForest
    ## we update the matrix of areas

    # Check the absolute value of the percent difference
    last_two_rows_diff <- ( abs(diff(area_mat[ (nrow(area_mat)-1):nrow(area_mat), ]))/area_mat[ (nrow(area_mat)-1), ])*100

    avg_diff <- mean(last_two_rows_diff)

    # Store avg_diff value to track each iteration
    avg_diff_values <- c(avg_diff_values, avg_diff)

    ## we update the final forest used with the most recent forest estimated in the convergence step
    # Update res@FinalForest with the forest in convergence_res@FinalForest
    res.strata4.1@FinalForest <- convergence_res@FinalForest

    ## we update the value with the most recent estimated value in the convergence steo
    res.strata4.1@value <- convergence_res@value

    ## we also update the long_data
    res.strata4.1@long_data <- convergence_res@long_data

    # Store res@valueTrain value to track each iteration
    valueTrain_values <- c(valueTrain_values, res.strata4.1@value)

    ## track these values in the forest output
    res.strata4.1@n_it <- conv_iterations
    res.strata4.1@avgKM_diff <- as.matrix(avg_diff_values)
    res.strata4.1@valueTrain_list <- valueTrain_values



    ### if the absolute change starts fluctuating and the last 10 iteration don't change by more than 0.5, we stop
    ## OR, if the forest itself has converging curves
    change_last <- mean(abs(diff(tail(as.matrix(avg_diff_values), 3))))


    ## wait until the absolute change in survival curves is less than something or the average change gets small

    if (avg_diff > 0.1) {
      if (is.na(change_last) || is.nan(change_last) || abs(change_last) > 1) {
        # If change_last5 is NaN or abs(change_last5) > 0.5, continue the loop
        continue_iterations <- TRUE

        cat("Difference in survival curve:", avg_diff, "\n")
        cat("Last 3 stages change:", change_last, "\n")

        # Increment the iteration counter
        conv_iterations <- conv_iterations + 1
      } else {
        # If avg_diff > 0.1 but abs(change_last5) <= 0.5, stop the loop
        continue_iterations <- FALSE

        cat("Difference in survival curve:", avg_diff, "\n")
        cat("Last 5 stages change:", change_last, "\n")

        # Update res@integral_KM with area_mat
        res.strata4.1@integral_KM <- area_mat

        # Update res@FinalForest with the forest in convergence_res@FinalForest
        res.strata4.1@FinalForest <- convergence_res@FinalForest

        # Update the value with the most recent estimated value in the convergence step
        res.strata4.1@value <- convergence_res@value

        # Update the long_data
        res.strata4.1@long_data <- convergence_res@long_data

        # Update the optimal output probabilities
        res.strata4.1@prev_probs <- convergence_res@prev_probs

        # Track these values in the forest output
        res.strata4.1@n_it <- conv_iterations
        res.strata4.1@avgKM_diff <- as.matrix(avg_diff_values)
        res.strata4.1@valueTrain_list <- valueTrain_values
      }
    } else {
      # If avg_diff <= 0.1, stop the loop
      continue_iterations <- FALSE

      cat("Difference in survival curve:", avg_diff, "\n")
      cat("Last 5 stages change:", change_last, "\n")

      # Update res@integral_KM with area_mat
      res.strata4.1@integral_KM <- area_mat

      # Update res@FinalForest with the forest in convergence_res@FinalForest
      res.strata4.1@FinalForest <- convergence_res@FinalForest

      # Update the value with the most recent estimated value in the convergence step
      res.strata4.1@value <- convergence_res@value

      # Update the long_data
      res.strata4.1@long_data <- convergence_res@long_data

      # Update the optimal output probabilities
      res.strata4.1@prev_probs <- convergence_res@prev_probs

      # Track these values in the forest output
      res.strata4.1@n_it <- conv_iterations
      res.strata4.1@avgKM_diff <- as.matrix(avg_diff_values)
      res.strata4.1@valueTrain_list <- valueTrain_values
    }
  }

  ### end of 4th strata forest

  ### start of 5th strata forest

  ## now, for strata3, we want to pool the observations together, treating each pts stage separately
  s1.strata5 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata that's complete
    data = long_data %>% filter(!!sym(paste0("strata",5)) == 1 & !is.na(as.character(attr(terms(models), "variables")[[2]][[2]]))),
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
    mTry = mTry[[1]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## we are in the first step of pooling patient data after running HC code
    pool1 = F,

    ## this is true if we are inputting survival probabilities
    appendstep1 = F
    #    inputpr = input.strata2[, colSums(input.strata2) != 0]
  )

  ## now, for strata 2, we want to track the optimal and eligibility
  ## first initialize a column name
  long_data$A.s1.strata5 <- NA

  eligibility_s1.strata5 <- s1.strata5@eligibility

  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE


  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE
  ### the eligibility is only for patients in strata 1 AND not NA, so we have to filter for that too
  long_data$A.s1.strata5[long_data[[paste0("strata", 5)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata5 == 1)] <- s1.strata5@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data[[paste0("strata", 5)]] == 1 & !is.na(long_data$T)][which(eligibility_s1.strata5 == 1)] <- s1.strata5@optimal@optimalTx



  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata5) )
  shiftedprobfinal[,eligibility_s1.strata5] <-t(s1.strata5@optimal@optimalY)



  ## the result after appyling the area function to each column of the matrix of survival probabilities
  ##### NOTE: the issue is that each patient has a different number of visits belonging to this strata
  ##### meaning, for visits in strata 2, we must input 0 probability
  areas <- apply(shiftedprobfinal, 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })

  ## now we want to create a matrix for areas where each row is one iteration of the forest
  area_mat <- areas


  # Construct an object of class DTRSurv for the output
  ## this needs to be initialized outside the loop
  res.strata5.1 <- new(
    Class = "DTRSurv",

    #####
    ###### delete these?
    ######
    "stageResults" = list(seq(1:nDP)),

    #####
    ###### delete these?
    ######

    "IHstageResults" = list(),
    "FinalForest" = s1.strata5,
    "value" = valueTrain,
    "call" = cl,
    "params" = params,
    "integral_KM" = area_mat,
    "n_it" = conv_iterations,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list(),
    "long_data" = long_data,

    ## NOTE: for the previous probs, we want to access res.strata(i - 1).1
    ############# THIS IS ANOTHER CHANGE TO MAKE


    "prev_probs" = res.strata4.1@prev_probs
  )

  ## after backwards recursion is complete, calculate the estimated value from the first stage
  ## calculates mean values of expected survival times and survival probabilities
  ## this is calculated across all PTS, so, once all pts have received their estimated optimal treatment --> what's the mean of all their survival times
  ## .meanValue() function defined in class_DTRSurvStep.R

  valueTrain <- .meanValue(object = s1.strata5)


  ## display the estimated value calculated in the first stage, and iterates through each element and prints names and values

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[i]], appendLF = FALSE)
  }


  # store values in call structure for returned object

  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  # Initialize a flag to indicate whether to continue iterations
  continue_iterations <- TRUE

  conv_iterations <- 1

  # Initialize vectors to store avg_diff and res@valueTrain values
  avg_diff_values <- c()
  valueTrain_values <- c(valueTrain)


  while(continue_iterations){

    message("Convergence Re-fitting Iteration:", conv_iterations)


    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv_otherstrata(data = data,
                                             prev.iteration = res.strata5.1, nDP = nDP, params = params, nTimes = nTimes,
                                             models = models, mTry = mTry, strata = 5,
                                             # use the most recent long_data
                                             long_data = res.strata5.1@long_data,
                                             prev_probs = res.strata5.1@prev_probs)

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <- rbind(res.strata5.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata5.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata5.1@prev_probs <- convergence_res@prev_probs


    ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
    #### meaning, for the same patient, if the difference between their estimated survival curves is large, we go through another iteration of training
    ## we update res@FinalForest with the forest in convergence_res@FinalForest
    ## we update the matrix of areas

    # Check the absolute value of the percent difference
    last_two_rows_diff <- ( abs(diff(area_mat[ (nrow(area_mat)-1):nrow(area_mat), ]))/area_mat[ (nrow(area_mat)-1), ])*100

    avg_diff <- mean(last_two_rows_diff)

    # Store avg_diff value to track each iteration
    avg_diff_values <- c(avg_diff_values, avg_diff)

    ## we update the final forest used with the most recent forest estimated in the convergence step
    # Update res@FinalForest with the forest in convergence_res@FinalForest
    res.strata5.1@FinalForest <- convergence_res@FinalForest

    ## we update the value with the most recent estimated value in the convergence steo
    res.strata5.1@value <- convergence_res@value

    ## we also update the long_data
    res.strata5.1@long_data <- convergence_res@long_data

    # Store res@valueTrain value to track each iteration
    valueTrain_values <- c(valueTrain_values, res.strata5.1@value)

    ## track these values in the forest output
    res.strata5.1@n_it <- conv_iterations
    res.strata5.1@avgKM_diff <- as.matrix(avg_diff_values)
    res.strata5.1@valueTrain_list <- valueTrain_values



    ### if the absolute change starts fluctuating and the last 10 iteration don't change by more than 0.5, we stop
    ## OR, if the forest itself has converging curves
    change_last <- mean(abs(diff(tail(as.matrix(avg_diff_values), 3))))


    ## wait until the absolute change in survival curves is less than something or the average change gets small

    if (avg_diff > 0.1) {
      if (is.na(change_last) || is.nan(change_last) || abs(change_last) > 1) {
        # If change_last5 is NaN or abs(change_last5) > 0.5, continue the loop
        continue_iterations <- TRUE

        cat("Difference in survival curve:", avg_diff, "\n")
        cat("Last 3 stages change:", change_last, "\n")

        # Increment the iteration counter
        conv_iterations <- conv_iterations + 1
      } else {
        # If avg_diff > 0.1 but abs(change_last5) <= 0.5, stop the loop
        continue_iterations <- FALSE

        cat("Difference in survival curve:", avg_diff, "\n")
        cat("Last 5 stages change:", change_last, "\n")

        # Update res@integral_KM with area_mat
        res.strata5.1@integral_KM <- area_mat

        # Update res@FinalForest with the forest in convergence_res@FinalForest
        res.strata5.1@FinalForest <- convergence_res@FinalForest

        # Update the value with the most recent estimated value in the convergence step
        res.strata5.1@value <- convergence_res@value

        # Update the long_data
        res.strata5.1@long_data <- convergence_res@long_data

        # Update the optimal output probabilities
        res.strata5.1@prev_probs <- convergence_res@prev_probs

        # Track these values in the forest output
        res.strata5.1@n_it <- conv_iterations
        res.strata5.1@avgKM_diff <- as.matrix(avg_diff_values)
        res.strata5.1@valueTrain_list <- valueTrain_values
      }
    } else {
      # If avg_diff <= 0.1, stop the loop
      continue_iterations <- FALSE

      cat("Difference in survival curve:", avg_diff, "\n")
      cat("Last 5 stages change:", change_last, "\n")

      # Update res@integral_KM with area_mat
      res.strata5.1@integral_KM <- area_mat

      # Update res@FinalForest with the forest in convergence_res@FinalForest
      res.strata5.1@FinalForest <- convergence_res@FinalForest

      # Update the value with the most recent estimated value in the convergence step
      res.strata5.1@value <- convergence_res@value

      # Update the long_data
      res.strata5.1@long_data <- convergence_res@long_data

      # Update the optimal output probabilities
      res.strata5.1@prev_probs <- convergence_res@prev_probs

      # Track these values in the forest output
      res.strata5.1@n_it <- conv_iterations
      res.strata5.1@avgKM_diff <- as.matrix(avg_diff_values)
      res.strata5.1@valueTrain_list <- valueTrain_values
    }
  }

  ### end of 5th strata forest

  ## we create this one if our nstrata == 2
  ###### NOTE: need to create a different result if nstrata == 1 and if nstrata == 5
  res <- new(
    Class = "DTRSurvRes",

    ## strata 1 forest

    "Forest1" = res.strata1.1,

    ## strata 2 forest
    "Forest2" = res.strata2.1,

    "Forest3" = res.strata3.1,
    "Forest4" = res.strata4.1,
    "Forest5" = res.strata5.1,

    "call" = cl,
    "params" = params,
    "long_data" = res.strata5.1@long_data,
    "prev_probs" = res.strata5.1@prev_probs,
    "n_it" = NA,
    "avgKM_diff" = matrix(0),
    "cutoff" = selected_thresh
  )

  # end of the if statement for nstrata == 5
  }else if (nstrata == 1){

  ## if nstrata == 1, we only use one forest
  ## for the other results, they are located depending on the value of nstrata
  #### if nstrata == 2, we can use the original two strata
  #### if nstrata == 5, we'll use all 5 strata

  ## now, we need to return both forests depending on the strata
  ### so, we create a new class called "DTRSurvRes"


  ## we create this one if our nstrata == 1
  ###### NOTE: need to create a different result if nstrata == 2 and if nstrata == 5
  res <- new(
    Class = "DTRSurvRes",

  ## strata 1 forest

    "Forest1" = res.strata1.1,

    "Forest2" = NA,
  "Forest3" = NA,
  "Forest4" = NA,
  "Forest5" = NA,

    "call" = cl,
    "params" = params,
    "long_data" = res.strata1.1@long_data,
    "prev_probs" = res.strata1.1@prev_probs,
    "n_it" = NA,
    "avgKM_diff" = matrix(0),

  ## there's no cutoff for the data split
    "cutoff" = 0
  )


}


  ## end of function
}
