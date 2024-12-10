
#' Indefinite Horizon Dynamic Treatment Regime Estimation for Survival Outcomes
#'
#' Provides methods for estimating multi-stage (indefinite horizon) optimal dynamic treatment
#'   regimes for survival outcomes, specifically truncated mean survival time, with dependent censoring.
#'
#'----------------- CHANGE THIS DESCRIPTION -----------------
#' We use a common formula for all decision points, i.e.,
#'   'models' is a single formula object; therefore, your data must follow a specific
#'   format. Specifically, 'stageLabel' = "_" is a requirement, so covariates must be named as
#'   xxx_1 for the first decision point,
#'   xxx_2 for the second, xxx_3 for the third, etc. The exact structure of the
#'   'xxx' can be generally defined; however, it cannot contain the stageLabel.
#'   For example, if the column names are
#'   (Y_1, Y_2, delta_1, delta_2, A_1, A_2, X_1, X_2)
#'   'models' = Surv(Y,delta) ~ X + A would lead to
#'   Surv(Y_1, delta_1) ~ X_1 + A_1 as the first stage model; and
#'   Surv(Y_2, delta_2) ~ X_2 + A_2
#'   as the second stage.
#'   The censoring indicator column, if using more than 1 strata, must be named "delta",
#'   which will be in the Surv( xx, delta). Otherwise, it can be named in any manner.
#'
#'   Y_k is the length of Stage k so that (Y_1 + Y_2 + ... + Y_K) is the
#'   overall observed failure time, delta_k is the censoring status at Stage k,
#'   d_k = 0 if a subject was censored at Stage k, and 1 if he/she experienced
#'   failure during that stage or moved to Stage k+1.
#'   A_k is the treatment at Stage k, k=1,2,..., K. Note that every quantity
#'   here is stage-wide. In other words, Y_2 is the length of Stage 2
#'   and is not cumulative from the baseline.
#'
#'   When one experienced censoring or failure at Stage k, it should be that
#'   Y_j = 0 for all j > k and instantaneous failure (Y_k < 1e-8) is not allowed;
#'   E.g., when delta_(k-1) = 1 and Y_k = 0, the person is considered dead at
#'   Stage k-1, but when d_(k-1) = 1 and Y_k = 2, the person made it to Stage k
#'   and either experienced failure or censoring (depending on d_k) during Stage k.
#'
#'   Any subject with missing values at Stage k will be ignored.
#'
#' @param ... Ignored. Present only to require named inputs.
#'
#' @param data A data.frame object. The full dataset including patient ID, treatments
#'   received, all stage covariates, observed times, and censoring
#'   indicators.
#'   Can be provided as a matrix object if column headers are included.
#'   Can contain missing data coded as NA, but cannot contain NaN.
#'   If performing analysis with multiple strata, cumulative.time is also required for strata splitting.
#'
#' @param txName A character vector object. The treatment variable name for
#'   each decision point. Each element corresponds to the respective decision
#'   point (element 1 = 1st decision; element 2 = 2nd decision, etc.).
#'
#' @param models A single formula. Contains a formula defining the
#'   response as a Surv() object and the covariate structure of the model.
#'   Note that this model should not include any terms of order > 1. This
#'  is a common formula to be used across all decision
#'   points.
#'   Example: Surv(visit.length, delta) ~ A + prev.visit.length
#'   In this example, visit.length is the observed time, delta is the censoring indicator, A is treatment.
#'   Any number of covariates can be added in the model.
#'
#'
#' @param timePoints A character object or a numeric vector object. If a character
#'   object, must be one of \{"quad", "uni", "exp"\} indicating the distribution
#'   from which the time points are to be calculated. For character input,
#'   input 'nTimes' must also be provided. If a numeric vector, the
#'   time points to be used. If 0 is not the first value, it will be
#'   concatenated by the software.
#'
#' @param nTimes An integer object. The total number of time points to be
#'   generated and considered. Used in conjunction with input 'timePoints'
#'   when 'timePoints' is a character; ignored otherwise.
#'
#' @param tau A numeric object or NULL. The study length. If NULL, the
#'   maximum timePoint is used.
#'
#' @param criticalValue A character object. Must be \{"mean"\}, which is the
#'  estimator for the value of a treatment rule which is the mean survival time.
#'
#' @param splitRule A character object OR NULL.
#'   Must be \{"mean"\} indicating the test used to determine an optimal split.
#'   Indicates usage of the truncated mean test.
#'
#' @param ERT A logical object. If TRUE, the Extremely Randomized Trees
#'   algorithm is used to select the candidate variable.
#'
#' @param sampleSize A numeric object, numeric vector object, or NULL.
#'   The fraction (0 < sampleSize <= 1) of the data to be used for each
#'   tree in the forest. If only
#'   one value is given, it is assumed to be the fraction for all decision
#'   points. If a vector is given, the length must be equal to the total
#'   number of decision points and each element corresponds to its respective
#'   decision point. If NULL and 'ERT' is TRUE,
#'   sampleSize defaults to 1.0. If NULL and 'ERT'
#'   is FALSE, sampleSize defaults to 0.632.
#'
#' @param uniformSplit A logical object. If 'ERT' and 'uniformSplit' are TRUE,
#'   the random cutoff is sampled from a uniform distribution over the range
#'   of available covariate values. If 'ERT' is TRUE and 'uniformSplit' is
#'   FALSE, a case is randomly selected and the cutoff is taken to be the mean
#'   cutoff between it and the next largest covariate value. If 'ERT' is FALSE,
#'   input is ignored.
#'
#' @param replace A logical object or NULL. If TRUE, the sample drawn for each
#'   of the nTree trees may have duplicate records. If FALSE, no individual is
#'   present in the sample for than once. If NULL, 'replace' = !'ERT'.
#'
#' @param randomSplit A numeric object. The probability that a random split
#'   will occur. Must be 0 < randomSplit < 1.
#'
#' @param tieMethod A character object. Must be one of
#'   \{"first", "random"\}. If multiple splits lead to the same
#'   value, the method by which the tie is broken.
#'
#' @param minEvent An integer object. The minimum number of events that must be
#'   present in a node.
#'
#' @param nodeSize An integer object. The minimum number of individuals that
#'   must be present in a node.
#'
#' @param nTree An integer object. The number of trees to grow.
#'
#' @param mTry An integer or integer vector object. The maximum number of
#'   covariates to sample for each split. If a vector, each element
#'   corresponds to its respective decision point.
#'
#'
#' @param stratifiedSplit A numeric object. The stratified random split
#'    coefficient. Covariates for which the number of splits (s_i) is less
#'    than s*stratifiedSplit/d are explored preferentially
#     (s is the total number of splits, d is the
#'    total number of covariates under consideration).
#'
#' @param stageLabel A character object. If using a common formula, the
#'    character used to separate the covariate from the decision point label.
#'    Must be "_".
#'    Example: "_", if the decision points are labeled Y_1, Y_2, Y_2
#'
#' @references She, J., Egberg, M., and Kosorok, M.R.
#'  An optimal dynamic treatment regime estimator for indefinite-horizon survival outcomes.
#'  Submitted.
#'
#' @include IH.VerifyData.R IH.VerifyTxName.R IH.VerifyModels.R IH.VerifySampleSize.R
#' @include IH.dtrSurvConverge.R IH.dtrSurvConverge_otherstrat.R
#' @include class_IH.Parameters.R
#' @include class_IH.DTRSurvStep.R class_IH.DTRSurv.R class_IH.DTRSurvRes.R
#' @import methods
#' @export
#' @useDynLib IHsurvrf
#' @import survival
#' @import tidyr
#' @import dplyr
#'
#' @returns An S4 object of class DTRSurvRes containing the key results and
#'   input parameters of the analysis.
#'
#' @seealso \code{\link{predict}} for retrieving the optimal treatment
#'    and/or the optimal survival curves.
#'
#' @examples
#'
#'
#' dt <- data.frame("Y_1" = sample(1:100,100,TRUE), "Y_2" = sample(1:100,100,TRUE),
#'                  "D_1" = rbinom(100, 1, 0.9), "D_2" = rbinom(100,1,0.9),
#'                  "A_1" = rbinom(100, 1, 0.5), "A_2" = rbinom(100,1,0.5),
#'                  "X_1" = rnorm(100), "X_2" = rnorm(100), subj.id = 1:100)
#'
#' 2 stage analysis with 1 strata
#' IHdtrSurv(data = dt,
#'         txName = c("A_1", "A_2"),
#'         stageLabel = "_",
#'         models = Surv(Y,D)~X+A,
#'         nstrata = 1)
#'
#'
#'
#'  5 stage analysis with 2 strata. We note that cumulative.time is required for more than 1 strata.
#'  We also note that due to the use of more than 1 strata, we need to name the censoring indicator "delta"
#'
#' dt <- data.frame("Y_1" = sample(1:100,100,TRUE), "Y_2" = sample(1:100,100,TRUE),
#'                  "Y_3" = sample(1:100,100,TRUE), "Y_4" = sample(1:100,100,TRUE),
#'                  "Y_5" = sample(1:100,100,TRUE),
#'                   "delta_1" = rbinom(100, 1, 0.9), "delta_2" = rbinom(100,1,0.9),
#'                  "delta_3" = rbinom(100,1,0.9), "delta_4" = rbinom(100,1,0.9),
#'                  "delta_5" = rbinom(100,1,0.9),
#'                   "A_1" = rbinom(100, 1, 0.5), "A_2" = rbinom(100,1,0.5),
#'                  "A_3" = rbinom(100,1,0.5), "A_4" = rbinom(100,1,0.5),
#'                  "A_5" = rbinom(100,1,0.5),
#'                   "X_1" = rnorm(100), "X_2" = rnorm(100),
#'                  "X_3" = rnorm(100), "X_4" = rnorm(100),
#'                  "X_5" = rnorm(100), subj.id = 1:100)
#'
#' dt$cumulative.time_1 <- 0
#' dt$cumulative.time_2 <- dt$Y_1
#' dt$cumulative.time_3 <- dt$Y_1 + dt$Y_2
#' dt$cumulative.time_4 <- dt$Y_1 + dt$Y_2 + dt$Y_3
#' dt$cumulative.time_5 <- dt$Y_1 + dt$Y_2 + dt$Y_3 + dt$Y_4
#'
#' IHdtrSurv(data = dt,
#'         txName = c("A_1", "A_2", "A_3", "A_4", "A_5"),
#'         stageLabel = "_",
#'         models = Surv(Y,delta)~X+A,
#'         nstrata = 2, windowsize = 5)
####################################################################################



IHdtrSurv <- function(data,
                      txName,
                      models,
                      ...,
                      timePoints = "uni",
                      nTimes = 200L,
                      tau = NULL,
                      criticalValue = "mean",
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
                      nstrata = 1,

                      ## the number of days we want to go back by from tau to help define strata
                      ## the smaller, the more accurate the division of the percentage of events in each strata
                      windowsize = 10) {
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
  if (!is.list(x = stagemodels)){stagemodels <- list(stagemodels)}

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

    ## nSamples is created (not input)
    nSamples = nSamples,
    stratifiedSplit = stratifiedSplit
  )

  tau <- params@tau


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
  variables <-
    unique(sub("_(\\d+)", "", setdiff(colnames(data), c(
      "subj.id", "rep.id"
    ))))


  ## reformatting the data in long form
  # Convert from wide to long format
  ## pivot_longer is function from tidyr


  long_data <-  data %>%
    pivot_longer(
      cols = starts_with(variables),
      names_sep = "_",
      names_to = c(".value", "stage"),
      names_prefix = "",
      values_to = "value",
      values_drop_na = FALSE
    ) %>%
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
    for (i in seq(from = starting_thresh,
                  to = 0,
                  by = -windowsize)) {
      # Select data where cumulative time is greater than or equal to the current threshold
      selected_data <-
        long_data %>% filter(cumulative.time * tau >= i)

      # Calculate the percentage of events in current cutoff divided by total number of events
      event_percentage <-
        sum(selected_data$delta == 0, na.rm = TRUE) / sum(long_data$delta == 0, na.rm = TRUE)

      # Check if the event percentage meets the threshold
      ## NOTE: that the dimensions of these will not match the original data just because the original data contains rows of NAs
      ## this is because we simulate every patient to have 10 stages but if they have event/censored earlier they just have NA
      if (event_percentage >= 0.3) {
        # Assign strata membership
        long_data$strata1[(long_data$cumulative.time * tau) >= i] <-
          1
        long_data$strata2[(long_data$cumulative.time * tau) < i] <-
          1

        ## subset the data
        strata1 <- selected_data
        strata2 <- long_data %>% filter(cumulative.time * tau < i)

        selected_thresh <- i

        ### print the cutoff time for the strata
        message("Cumulative cutoff:",
                selected_thresh,
                "; longest length:",
                tau)


        break


      }
    }
  }


  # Step 1: Dynamically get column names starting with "strata"
  strata_columns <-
    names(long_data)[grepl("^strata", names(long_data))]


  # Step 2: Convert back to wide format
  data <- long_data %>%
    pivot_wider(
      id_cols = subj.id,
      names_from = stage,
      values_from = c(variables,  all_of(strata_columns)),
      names_sep = "_"
    )

  # Get the column name as a string
  respname <-
    as.character(attr(terms(models), "variables")[[2]][[2]])

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
  long_data$A.s1.strata1[long_data$strata1 == 1 &
                           !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata1 == 1)] <-
    s1.strata1@optimal@optimalTx
  ## also update the "A" column to be used for the convergence aspect
  long_data$A[long_data$strata1 == 1 &
                !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata1 == 1)] <-
    s1.strata1@optimal@optimalTx


  # Initialize the probability matrix that's output
  shiftedprobfinal <-
    matrix(NA,
           nrow = nTimes,
           ncol = length(eligibility_s1.strata1))
  shiftedprobfinal[, eligibility_s1.strata1] <-
    t(s1.strata1@optimal@optimalY)


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

  while (continue_iterations && conv_iterations <= max_iter) {
    message("Convergence Re-fitting Iteration:", conv_iterations)

    ####### Here, we will construct a new iteration of forest training to check for convergence
    convergence_res <- IHdtrConv(
      data = data,
      prev.iteration = res.strata1.1,
      nDP = nDP,
      params = params,
      nTimes = nTimes,
      models = models,
      mTry = mTry,
      strata = 1,
      long_data = res.strata1.1@long_data,
      prev_probs = res.strata1.1@prev_probs
    )

    ######## for the pt that was zeroed out, we want to get rid of their integral from the original res.strata1

    ## now we want to create a matrix for areas where each row is one iteration of the forest
    area_mat <-
      rbind(res.strata1.1@integral_KM, convergence_res@integral_KM)

    # Update res@integral_KM with area_mat
    res.strata1.1@integral_KM <- area_mat

    # update the optimal output probabilities
    res.strata1.1@prev_probs <- convergence_res@prev_probs

    ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
    last_two_rows_diff <-
      (abs(diff(area_mat[(nrow(area_mat) - 1):nrow(area_mat), ])) / area_mat[(nrow(area_mat) -
                                                                                1), ]) * 100

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
    if (conv_iterations > max_iter) {
      continue_iterations <-FALSE  # Stop the loop when max_iter is reached
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
      data = long_data %>% filter(!!sym(paste0("strata", 2)) == 1 &
                                    !is.na(!!sym(resp_name))),
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
    long_data$A.s1.strata2[long_data[[paste0("strata", 2)]] == 1 &
                             !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata2 == 1)] <-
      s1.strata2@optimal@optimalTx
    ## also update the "A" column to be used for the convergence aspect
    long_data$A[long_data[[paste0("strata", 2)]] == 1 &
                  !is.na(long_data[[as.character(attr(terms(models), "variables")[[2]][[2]])]])][which(eligibility_s1.strata2 == 1)] <-
      s1.strata2@optimal@optimalTx



    # Initialize the shifted probability matrix
    ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
    ## we overwrite the eligible patients
    shiftedprobfinal <-
      matrix(NA,
             nrow = nTimes,
             ncol = length(eligibility_s1.strata2))
    shiftedprobfinal[, eligibility_s1.strata2] <-
      t(s1.strata2@optimal@optimalY)



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


    while (continue_iterations && conv_iterations <= max_iter) {
      message("Convergence Re-fitting Iteration:", conv_iterations)

      ####### Here, we will construct a new iteration of forest training to check for convergence
      convergence_res <- IHdtrConv_otherstrata(
        data = data,
        prev.iteration = res.strata2.1,
        nDP = nDP,
        params = params,
        nTimes = nTimes,
        models = models,
        mTry = mTry,
        strata = 2,
        # use the most recent long_data
        long_data = res.strata2.1@long_data,
        prev_probs = res.strata2.1@prev_probs
      )
      ######## for the pt that was zeroed out, we want to get rid of their integral from the original res.strata1

      ## now we want to create a matrix for areas where each row is one iteration of the forest
      area_mat <-
        rbind(res.strata2.1@integral_KM,
              convergence_res@integral_KM)

      # Update res@integral_KM with area_mat
      res.strata2.1@integral_KM <- area_mat

      # update the optimal output probabilities
      res.strata2.1@prev_probs <- convergence_res@prev_probs

      ## if the difference between these rows on average is greater than 0.005, then we go through another iteration
      last_two_rows_diff <-
        (abs(diff(area_mat[(nrow(area_mat) - 1):nrow(area_mat), ])) / area_mat[(nrow(area_mat) -
                                                                                  1), ]) * 100

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
      if (conv_iterations > max_iter) {
        continue_iterations <-
          FALSE  # Stop the loop when max_iter is reached
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
  } else if (nstrata == 1) {
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
  # end of function
}

#' Hidden methods
#'
#' @name IHdtrSurv-internal-api
#' @keywords internal
#' @import methods
NULL
