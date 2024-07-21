
setwd("~/survrf/Scripts/IHsurvrf")
source("R/class_IH.DTRSurv.R")

source("R/IH.VerifyData.R")
source("R/IH.VerifyTxName.R")
source("R/IH.VerifyUsePrevTime.R")
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
                         usePrevTime = TRUE,
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
                         pooled = FALSE,
                         stratifiedSplit = NULL,
                         stageLabel = ".",
                         stage.start = NULL,
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

  if (!is.null(stage.start)) {
    nDP <- stage.start
  }

  # ensures that the usePrevTime input is of appropriate type (a logical)

  ## defined in VerifyUsePrevTime.R script

  usePrevTime <- .VerifyUsePrevTimes(usePrevTimes = usePrevTime)

  ## used in script VerifyModels.R
  ## this can be input as a list ofmodels, or a single formula (in which we probably assume common formula across all stages)

  stagemodels <- .VerifyModels(
    models = models,
    nDP = nDP,
    data = data,
    txName = txName,
    stageLabel = stageLabel,
    usePrevTime = usePrevTime
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
    pooled = pooled,
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
    arrange(subj.id, stage) %>%

    ## create a column called "A.opt.HC" to hold the optimal actions calculated from HC's method
    mutate(A.opt.HC = NA)

  ## we want to define the cumulative percentage of events that have occurred (delta = 0)
  ## when we go back a time window, we want to include in the upper strata the stages for patients who have cumulative time higher


  # Initialize starting point
  starting_thresh <- tau - windowsize

  # Initialize strata columns
  long_data$strata1 <- 0
  long_data$strata2 <- 0

  # Iterate over the data using a for loop
  for (i in seq(from = starting_thresh, to = 0, by = -windowsize)) {
    # Select data where cumulative time is greater than or equal to the current threshold
    selected_data <- long_data %>% filter(cumulative.time * tau >= i)

    # Calculate the percentage of events in current cutoff divided by total number of events
    event_percentage <- sum(selected_data$delta == 0, na.rm = TRUE) / sum(long_data$delta == 0, na.rm = TRUE)

    # Check if the event percentage meets the threshold
    ## NOTE: that the dimensions of these will not match the original data just because the original data contains rows of NAs
    ## this is because we simulate every patient to have 10 stages but if they have event/censored earlier they just have NA
    if (event_percentage >= 0.5) {
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


  # Step 2: Convert back to wide format
  data <- long_data %>%
    pivot_wider(
      id_cols = c(subj.id),
      names_from = stage,
      values_from = c(T, delta, A, baseline1, baseline2, state1, state2, prior.visit.length, cumulative.time, nstages, action.1.count, action.0.count, strata1, strata2),
      names_sep = "_"
    )


  ### now, for strata1, we want to pool the observations together, treating each pts stage separately
  ### fit a forest using this

  s1.strata1 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata with complete cases (don't include if they've failed already)
    data = long_data %>% filter(strata1 == 1 & !is.na(T)),
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
  long_data$A.s1.strata1[long_data$strata1 == 1 & !is.na(long_data$T)][which(eligibility_s1.strata1 == 1)] <- s1.strata1@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data$strata1 == 1 & !is.na(long_data$T)][which(eligibility_s1.strata1 == 1)] <- s1.strata1@optimal@optimalTx


  # Initialize the probability matrix that's output
  shiftedprobfinal <- matrix(NA, nrow = nTimes, ncol = length(eligibility_s1.strata1) )
  shiftedprobfinal[,eligibility_s1.strata1] <-t(s1.strata1@optimal@optimalY)

## note: we don't know the final stage's probability since each patient's final stage may be different


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
  area_mat <- areas

  # Initialize vectors to store avg_diff and res@valueTrain values
  avg_diff_values <- c()
  valueTrain_values <- c(valueTrain)

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
    "stageResults" = list(seq(1:10)),

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


  while(continue_iterations){

    message("Convergence Re-fitting Iteration:", conv_iterations)


    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv(data = data,
                                   prev.iteration = res.strata1.1, nDP = nDP, params = params, nTimes = nTimes,
                                   models = models, mTry = mTry, strata = 1, long_data = res.strata1.1@long_data,
                                   prev_probs = res.strata1.1@prev_probs)

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <- rbind(res.strata1.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata1.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata1.1@prev_probs <- convergence_res@prev_probs


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
    res.strata1.1@FinalForest <- convergence_res@FinalForest

    ## we update the value with the most recent estimated value in the convergence steo
    res.strata1.1@value <- convergence_res@value

    ## we also update the long_data
    res.strata1.1@long_data <- convergence_res@long_data

    # Store res@valueTrain value to track each iteration
    valueTrain_values <- c(valueTrain_values, res.strata1.1@value)

    ## track these values in the forest output
    res.strata1.1@n_it <- conv_iterations
    res.strata1.1@avgKM_diff <- as.matrix(avg_diff_values)
    res.strata1.1@valueTrain_list <- valueTrain_values




    ## wait until the absolute change (not avg is less than 0.01%)

    if(avg_diff > 0.1) {
      # If the condition is met, continue the loop
      continue_iterations <- TRUE

      print(avg_diff)



      ## increment the iteration counter

      conv_iterations <- conv_iterations + 1

    } else {
      # If the condition is not met, stop the loop
      continue_iterations <- FALSE


      # Update res@integral_KM with area_mat
      res.strata1.1@integral_KM <- area_mat


      ## we update the final forest used with the most recent forest estimated in the convergence step
      # Update res@FinalForest with the forest in convergence_res@FinalForest
      res.strata1.1@FinalForest <- convergence_res@FinalForest

      ## we update the value with the most recent estimated value in the convergence steo
      res.strata1.1@value <- convergence_res@value

      ## we also update the long_data
      res.strata1.1@long_data <- convergence_res@long_data

      # update the optimal output probabilities
      res.strata1.1@prev_probs <- convergence_res@prev_probs

      ## track these values in the forest output
      res.strata1.1@n_it <- conv_iterations
      res.strata1.1@avgKM_diff <- as.matrix(avg_diff_values)
      res.strata1.1@valueTrain_list <- valueTrain_values

    }

  }


  ###########
  ########### now, we use the previously optimized convergence probabilities to use as input for the new forest
  ## for the observations in strata 2, if their k+1 stage is in strata 1, we append the strata 2 observation
  ## we create a bunch of "double stubs" but not through prediction, through appending
  #### for patients whose next stage is NOT in strata 1, we just treat them like in HC code where the
  ## the first row of the survival probability is set to 1
  ## then, we convert these by subtracting from the row above so it's a difference in probs

  ## first, we initialize a matrix of probabilities, with ncol == each patient gets nDP stages

  input.strata2 <- matrix(0, nrow = nTimes, ncol = nrow(long_data %>% filter(stage == 1))*nDP )


  ## for stage 10, get the current eligibility; if the patient is in strata 2, give them first row == 1
  stage.elig <- long_data %>%
    # Filter the data to include only rows where stage is equal to i
    filter(stage == nDP) %>%

    # Mutate the data to add a new column 'eligibility'
    mutate(
      eligibility = as.numeric(
        # Use ifelse to check two conditions for assigning eligibility
        ifelse(
          # Condition 1: Check if strata 2 == 1
          !!sym(paste0("strata", 2)) == 1 &
            # Condition 2: Check if the row has complete cases excluding columns starting with "A"
            complete.cases(dplyr::select(., -matches("^A\\."))),
          # If both conditions are TRUE, assign 1
          1,
          # Otherwise, assign 0
          0
        )
      )
    ) %>%

    # Extract the 'eligibility' column as a vector
    pull(eligibility)


  ## now, for every 10th column (equivalent to pulling all pt stage 10), we overwrite with 1 in first row for eligible pt
  ## meaning, for all patients who are eligible in stage 10 (and have non-NA values, we know they survived so we put 1)
  input.strata2[, seq( (nDP),
                       ncol(input.strata2), by = (nDP))][, which(stage.elig == 1)][1L, ] <- 1.0

  ## now, we loop through each stage:

  for (i in (nDP-1):1){

    message("Prep for strata 2, appending for stage", i)

    ## for the current stage, we first get an eligibility for the current observations if theyy're in strata 2

    stage.elig <- long_data %>%
      # Filter the data to include only rows where stage is equal to i
      filter(stage == i) %>%

      # Mutate the data to add a new column 'eligibility'
      mutate(
        eligibility = as.numeric(
          # Use ifelse to check two conditions for assigning eligibility
          ifelse(
            # Condition 1: Check if the strata column (constructed dynamically) is equal to 1
            !!sym(paste0("strata", 2)) == 1 &
              # Condition 2: Check if the row has complete cases excluding columns starting with "A"
              complete.cases(dplyr::select(., -matches("^A\\."))),
            # If both conditions are TRUE, assign 1
            1,
            # Otherwise, assign 0
            0
          )
        )
      ) %>%

      # Extract the 'eligibility' column as a vector
      pull(eligibility)




    next.stage.elig <- long_data %>%
      # Filter the data to include only rows where stage is equal to i
      filter(stage == (i+1)) %>%

      # Mutate the data to add a new column 'eligibility'
      mutate(
        eligibility = as.numeric(
          # Use ifelse to check two conditions for assigning eligibility
          ifelse(
            # Condition 1: Check if the strata column (constructed dynamically) is equal to 1
            !!sym(paste0("strata", 1)) == 1 &
              # Condition 2: Check if the row has complete cases excluding columns starting with "A"
              complete.cases(dplyr::select(., -matches("^A\\."))),
            # If both conditions are TRUE, assign 1
            1,
            # Otherwise, assign 0
            0
          )
        )
      ) %>%

      # Extract the 'eligibility' column as a vector
      pull(eligibility)



    ## if both are true (both == 1), then we need to extract the survival probability, and append the observed
    ## meaning, the current stage has a patient in strata 2, next stage is in strata 1
    ## this then is the input survival probability
    app.elig <- ifelse(stage.elig == 1 & next.stage.elig == 1, TRUE, FALSE)

    ## first, we subset the stage i survival matrix, then overwrite it with the optimal for cases where this is true
    ## we overwrite it with the eligible
    input.strata2[, seq( i, ncol(input.strata2), by = nDP)][, which(app.elig == 1)] <-

      res.strata1.1@prev_probs[, seq( (i+1), ncol(input.strata2), by = nDP)][, which(app.elig == 1)]

    ## then, for these patients, we append, but don't transform into difference of probabilities (we do this manually after)
    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)

    # extract response and delta from model frame

    ## extract survival response for all 300 patients
    ### however, those that are not in elig_append1 = TRUE will receive an NA
    response_append1 <- stats::model.response(data = x_append1)
    ## if there are any that are not eligible for this stage, we turn into NA
    response_append1[!app.elig, ] <- NA

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]
    delta_append1[!app.elig] <- NA

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]
    response_append1[!app.elig] <- NA

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

    app.elig <- app.elig & !zeroed


    ### if sum(app.elig) is 0, we skip over this, we should just leave the probabilities as 0

    if (sum(app.elig) != 0){

    ## now, we overwrite these probabilities for eligible patients with appended probs
    input.strata2[, seq( i, ncol(input.strata2), by = nDP)][, which(app.elig == 1)] <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = input.strata2[, seq( i, ncol(input.strata2), by = nDP)][, app.elig, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response_append1[app.elig],

      ## probably transforming survival times into probabilities?
      surv2prob = FALSE
    )

    }


    #### we create another eligibility if the next stage is in the same strata
    next.stage.same.strata <- long_data %>%
      # Filter the data to include only rows where stage is equal to i
      filter(stage == (i+1)) %>%

      # Mutate the data to add a new column 'eligibility'
      mutate(
        eligibility = as.numeric(
          # Use ifelse to check two conditions for assigning eligibility
          ifelse(
            # Condition 1: Check if the strata column (constructed dynamically) is equal to 1
            !!sym(paste0("strata", 2)) == 1 &
              # Condition 2: Check if the row has complete cases excluding columns starting with "A"
              complete.cases(dplyr::select(., -matches("^A\\."))),
            # If both conditions are TRUE, assign 1
            1,
            # Otherwise, assign 0
            0
          )
        )
      ) %>%

      # Extract the 'eligibility' column as a vector
      pull(eligibility)




    ## otherwise, if the current stage is in strata 2, but the next stage isn't in strata 1, we just use vector (1, 0000)
    ## also, if the current stage is in strata 2, and there is no next stage
    no.app.1 <- ifelse(stage.elig == 1 & next.stage.elig == 0, TRUE, FALSE)

    ## now, for these patients, we just change their first row to 1
    ## now, for every 10th column (equivalent to pulling all pt stage 10), we overwrite with 1 in first row for eligible pt
    input.strata2[, seq( (i),
                         ncol(input.strata2), by = (nDP))][, which(no.app.1 == 1)][1L, ] <- 1.0


    ## note, that patients who are not in the strata will eventually be filtered out as we will only use cols where
    ## the sum of the column is not 0



  }

  ## now, we need to shift these probabilities to become differences
  input.strata2 <- input.strata2 - rbind( as.matrix(input.strata2[-1L,]), 0.0)

  ## sets very small values in pr to 0

  input.strata2[abs(input.strata2) < 1e-8] <- 0.0



  ## now, for strata2, we want to pool the observations together, treating each pts stage separately
  s1.strata2 <- .dtrSurvStep(
    ## use the model for pooled data
    model = models,
    ## only include the data in the first strata that's complete
    data = long_data %>% filter(strata2 == 1 & !is.na(T)),
    priorStep = NULL,
    params = params,
    txName = "A",
    ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
    mTry = mTry[[1]],
    ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
    sampleSize = 1,
    ## we are in the first step of pooling patient data after running HC code
    pool1 = F,
    appendstep1 = F
    #inputpr = input.strata2[, colSums(input.strata2) != 0]
  )

  ## now, for strata 2, we want to track the optimal and eligibility
  ## first initialize a column name
  long_data$A.s1.strata2 <- NA

  eligibility_s1.strata2 <- s1.strata2@eligibility

  # Update actions for A.s1.strata1 for rows where strata1 == 1 and eligibility_final is TRUE

  ###############
  ################ check to see if the eligibility matches with the dimensions... might need filter for NA (T)?
  ###############

  long_data$A.s1.strata2[which(long_data$strata2 == 1)][which(eligibility_s1.strata2 == 1)] <- s1.strata2@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[which(long_data$strata2 == 1)][which(eligibility_s1.strata2 == 1)] <- s1.strata2@optimal@optimalTx

  #################
  #################
  #################


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
    "stageResults" = list(seq(1:10)),

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




    ## wait until the absolute change (not avg is less than 0.01%)

    if(avg_diff > 2) {
      # If the condition is met, continue the loop
      continue_iterations <- TRUE


      print(avg_diff)


      ## increment the iteration counter

      conv_iterations <- conv_iterations + 1

    } else {
      # If the condition is not met, stop the loop
      continue_iterations <- FALSE


      # Update res@integral_KM with area_mat
      res.strata2.1@integral_KM <- area_mat


      ## we update the final forest used with the most recent forest estimated in the convergence step
      # Update res@FinalForest with the forest in convergence_res@FinalForest
      res.strata2.1@FinalForest <- convergence_res@FinalForest

      ## we update the value with the most recent estimated value in the convergence steo
      res.strata2.1@value <- convergence_res@value

      ## we also update the long_data
      res.strata2.1@long_data <- convergence_res@long_data

      # update the optimal output probabilities
      res.strata2.1@prev_probs <- convergence_res@prev_probs

      ## track these values in the forest output
      res.strata2.1@n_it <- conv_iterations
      res.strata2.1@avgKM_diff <- as.matrix(avg_diff_values)
      res.strata2.1@valueTrain_list <- valueTrain_values

    }

  }

  ## now, we need to return both forests depending on the strata
  ### so, we create a new class called "DTRSurvRes"

  # Construct an object of class DTRSurv for the output
  ## this needs to be initialized outside the loop
  res <- new(
    Class = "DTRSurvRes",

  ## strata 1 forest

    "Forest1" = res.strata1.1,

  ## strata 2 forest
    "Forest2" = res.strata2.1,
    "call" = cl,
    "params" = params,
    "long_data" = res.strata2.1@long_data,
    "prev_probs" = res.strata2.1@prev_probs,
    "n_it" = NA,
    "avgKM_diff" = matrix(0),
    "cutoff" = selected_thresh
  )






  ## end of function
}
