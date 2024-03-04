


# Class for storing primary results from a single stage of the Q-learning
#   survival analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot txName A character object. The name of the treatment variable for the
#    Q-learning step
#
#  @slot txLevels A vector. The treatment options available.
#
#  @slot model A formula object. The model for extracting the covariates to be
#    considered for splitting
#
#  @slot survRF A SurvRF object. The primary results of the tree building
#    algorithm
#
#  @slot eligibility A logical vector object. TRUE indicates that the case was
#    eligible to be included in analysis
#
#  @slot valueAllTx A list object. The value of the tree for each tx level.
#
#  @slot optimal An Optimal object. The estimated optimal tx and optimal value.


library(stats)
source("R/IH.survRF.R")
source("R/class_IH.Optimal.R")
source("R/IH.shiftMat.R")

## defines a new S4 class called DTRSurvStep for storing a single stage of Q-learning survival analysis
setClass(
  Class = "DTRSurvStep",

  ## slots (contain attributes/properties) that objects of this class contains

  slots = c(
    ## txName: a character object- stores name of treatment variable
    "txName" = "character",

    ## vector that contains available treatment options
    "txLevels" = "vector",

    ## holds model formula for extracting covariates to be considered for RF splitting
    "model" = "formula",

    ## stores SurvRF object-- contains primary results of tree building algorithm
    "survRF" = "ANY",

    ## logical vector (T/F): T= case eligible to be included
    "eligibility" = "logical",

    ## stores value of tree for each treatment level as a list
    "valueAllTx" = "list",

    ## contains estimated optimal treatment and its value
    "optimal" = "Optimal"
  )
)


#-------------------------------------------------------------------------------
# Internal function to create a random forest for a single Q-learning step
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Internal function to create random forest
#
# Function is not exported
#
# @param model: A survival formula object, the rhs of which specifies the
#   covariates to be considered in the splitting algorithm
#
# @params data: A data.frame object containing covariate and treatment histories
#
#
# @params: priorStep A DTRSurvStep object. The analysis from a previous step
#
# @params: params A Parameters object.
#
# @params: txName A character object or a character vector object. The names
#   of the treatment variables in data
#
# @params: mTry An integer object or NULL. If integer, the maximum number of
#    covariates to consider for each split. If NULL, mTry is sqrt(np)
#
# @params: sampleSize An integer object. The number of samples to draw for each
#    tree


## specify that certain functions from a package should be imported into namespace of package
## several functions from the "stats" package are being imported



## TimeInfo script sourced: survRF which sources Parameters which sources TimeInfo


## define an internal function .dtrSurvStep
.dtrSurvStep <- function(...,
                         ## input arguments
                         model,
                         data,
                         priorStep,
                         params,
                         txName,
                         mTry,
                         sampleSize) {
  # model is input
  mod <- model

  ## extract first order terms (main effect) from the model for splitting in the random forest algorithm

  # identify order 1 terms in formula
  ## == 1L checks if terms are of order 1, result is a logical vector
  order1 <- attr(x = stats::terms(x = mod), which = "order") == 1L
  if (any(order1)) {
    ## extract the labels from the first order terms & stores them (names of the covariates used in the model)
    stageCov <-
      attr(x = stats::terms(x = mod), which = "term.labels")[order1]

    ## if there are no first order terms
  } else {
    ## stop with an error
    stop("problem in identifying covariates, verify formula\n",
         call. = FALSE)
  }

  ## identifies interaction terms, but are ignored as the message output suggests

  # warn about order > 1
  orderHigh <- attr(x = stats::terms(x = mod), which = "order") > 1L
  if (any(orderHigh))
    message("interaction terms are ignored")

  # extract model frame

  ## prepares model frame using the formula (input) and the data (input)
  ## missing values are not specifically handled (not omitted)

  x <-
    stats::model.frame(formula = mod,
                       data = data,
                       na.action = na.pass)

  ## generates message to display the model formula used in the analysis

  message("model ", appendLF = FALSE)
  tm <- as.character(mod)
  message(tm[2], " ~ ", tm[3])

  ## determines eligible cases based on complete cases in x (model frame)

  # identify individuals with complete data
  elig <- stats::complete.cases(x)

  # extract response and delta from model frame

  ## extract survival response
  response <- stats::model.response(data = x)

  ## extract censoring indicator (delta) from the second column of the "response" data
  ## "L" is used to indicate that 2 is an integer
  delta <- response[, 2L]

  ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
  response <- response[, 1L]

  # remove response from x

  ## if first column of the model frame (x) is the response variable, remove this column
  ## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
  if (attr(x = terms(x = mod), which = "response") == 1L) {
    x <- x[,-1L, drop = FALSE]
  }

  ## marks zeroed survival times and updates eligibility

  # responses that are zero (effectively) indicate censored at a previous stage

  ## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
  zeroed <- abs(x = response) < 1e-8

  ## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

  elig <- elig & !zeroed

  ## if there are no eligible cases, output an error message

  if (sum(elig) == 0L)
    stop("no cases have complete data", call. = FALSE)

  ## displays message indicating the number of cases that are still eligible for analysis at this stage

  message("cases in stage: ", sum(elig))


  ## checks if mTry (input) is specified (mTry is the argument for x)
  ## if it's NULL, it's set to the square root of the number of columns in x, rounded up
  ## this is common in RF modeling, where this is often the default value of mTry

  # maximum number of covariates to try
  if (is.null(x = mTry)) {
    mTry <- as.integer(x = ceiling(x = sqrt(x = ncol(x = x))))
    message("maximum # of covariates considered for splitting set to ", mTry)

    ## if specified mTry > number of covariates (columns in x)
    ## indicate that the value is too large, and resets it to square root of n columns in x

  } else if (mTry > ncol(x = x)) {
    message("mTry reset as it is larger than the # of available covariates")
    mTry <- as.integer(x = ceiling(x = sqrt(x = ncol(x = x))))
    message("maximum # of covariates considered for splitting set to ", mTry)
  } else {
    ## if mTry is specified and not greater than number of covariates, converted to an integer
    mTry <- as.integer(x = mTry)
  }



  ## if priorStep (input) is null, AKA there was no previous step
  ##       this means we are currently in the last stage

  if (is.null(x = priorStep)) {
    # priorStep is NULL for first step of the analysis.
    # transform the time variable to a probability mass vector

    ## applies function FUN to each element of X which checks if each survival time (s) is less than each time point
    ## .TimePoints() part specifies an additional argument to be passed to fun
    ## .TimePoints() defined in class_TimeInfo.R

    # identify time points <= response
    tSurv <- sapply(
      X = response[elig],

      # results in a logical integer (0 = FALSE)
      FUN = function(s, tp) {
        as.integer(x = {
          s < tp
        })
      },
      tp = .TimePoints(object = params)
    )

    # time point nearest the response without going over
    # {nTimes x nElig}

    ## pr: probability mass vector which binds transformed survival times (1, 0) (excluding first row) w values with a row of 1s
    ## for tSurv,each element is 1 if the survival time < corresponding time point. Otherwise, it's 0
    ## then subtract the transformed survival times to create a probability mass vector
    ## this calculates the difference between successive elements in each row of tSurv by translating binary indicators into p.survival until each time point
    ## there's now mass 1 at time point where each survival time is censored or an event, and 0 for all other time points

    pr <- {
      rbind(tSurv[-1L,], 1) - tSurv
    }

    ## if priorStep isn't null,

  } else {
    # priorStep is not NULL when q < Q
    ### priorStep is of class DTRSurvStep AKA the results from the prior step of the analysis
    #
    # the number of timepoints
    # .NTimes() is a getter method defined for Parameters objects

    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases

    ## initializes a survival matrix called "survMatrix" with 0's

    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = nrow(x = x))

    ## sets the first row to 1.0

    survMatrix[1L,] <- 1.0

    ## updates survMatrix with survival functions estimated from the previous step
    ## priorStep: A DTRSurvStep object. The analysis from a previous step

    # retrieve estimated OPTIMAL survival function from previous step
    ## then accesses the "eligibility" slot (logical) to select only the columns where the patient is still eligible
    ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
    ## transpose the output of .OptimalY to align with the structure
    ## .OptimalY defined in class_Optimal.R
    ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
    ## optimal slot is  object of class optimal
    ## .OptimalY acts to return the optimalY slot

    survMatrix[, priorStep@eligibility] <-
      t(.OptimalY(object = priorStep@optimal))

    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

    pr <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = response[elig],

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )

    ## sets very small values in pr to 0

    pr[abs(pr) < 1e-8] <- 0.0

  }


  ## if there are any NA values in the matrix, stop function execution and display error message
  if (any(is.na(x = pr)))
    stop("NA not permitted in pr -- contact maintainer",
         call. = FALSE)

  ## checks if any values in pr are outside the range [0,1]. If so, stops with error message
  ## pr is expected to be a matrix of probabilities and therefore between 0 and 1

  if (any(pr > 1.0) || any(pr < 0.0)) {
    stop("pr must obey 0 <= pr <= 1 -- contact maintainer", call. = FALSE)
  }

  ## if txName variable is a factor,

  # identify tx levels in limited data
  if (is.factor(x = data[, txName])) {
    ## extract levels of the factor treatment variable for eligible cases
    txLevels <- levels(x = factor(x = data[elig, txName]))
  } else {
    ## if txName is not a factor, fines unique values of the treatment for eligible cases & sort them to use as treatment levels
    txLevels <- sort(x = unique(x = data[elig, txName]))
  }

  ## checks if there's only one treatment level in the data
  ## if so, display a message to the user-- this is important bc analyses often require multiple treatment levels for comparison

  if (length(x = txLevels) == 1L) {
    message("***only one treatment level in data***")
  }

  ## if .Pooled (defined in class_TreeConditions.R) which returns a logical vector,
  ## pooled analysis applies survival RF model to the entire dataset

  if (.Pooled(object = params)) {
    ## displays message indicating pooled analysis is being conducted & listing the treatment levels involved

    message("pooled analysis; treatments ",
            paste(txLevels, collapse = " "))

    # this will be a SurvRF object

    ## applies the .survRF function for survival RF modeling
    ## coded in survRF.R script which generates a forest using the shifted times
    ## returns object of class "SurvRF"

    result <- .survRF(
      ## selects subset of predictor matrix for eligible cases & ensures result is a matrix
      x = x[elig, , drop = FALSE],
      y = response[elig],
      pr = pr,
      delta = delta[elig],
      params = params,
      mTry = mTry,
      txLevels = txLevels,
      model = mod,
      sampleSize = sampleSize
    )

    ## if not pooled,

  } else {
    ## a stratified analysis is conducted, so the analysis is performed separately for each treatment level

    message("stratified analysis")

    ## initialize an empty list "result" to store output from each treatment level

    # result will be a list of SurvRF objects
    result <- list()

    ## iterates over each treatment level "txLevels"

    for (i in 1L:length(x = txLevels)) {
      ## displays message indicating current treatment level being processed

      message("  treatment level ", txLevels[i])

      ## name of the current treatment level

      nms <- as.character(x = txLevels[i])

      ## logical vector indicating which cases the data corresponds to current treatment level

      di <- {
        data[elig, txName] == txLevels[i]
      }

      ## logical vector combining eligibility & current treatment level

      use <- elig & {
        data[, txName] == txLevels[i]
      }

      ## applies the .survRF function to the subset of data for the current treatment level

      result[[nms]] <- .survRF(
        x = x[use, , drop = FALSE],
        y = response[use],
        pr = pr[, di],
        delta = delta[use],
        params = params,
        mTry = mTry,
        txLevels = txLevels[i],
        model = mod,
        sampleSize = sampleSize
      )

    }

    ## converts list of Survrf objects into a single object of class "SurvRFStratified"
    ## this is defined in class_IH.SurvRF.R

    result <- new(Class = "SurvRFStratified", "strat" = result)

  }


  # calculate the estimated values for all treatment levels
  # .PredictAll() is a method; called here for objects of class SurvRF (or SurvRF stratified), which
  # is defined in file class_SurvRF.R line 318


  ## method with data.frame calculates value for new data from SurvRF object
  ## .PredictAll returns list object containing survFunc, mean, and? survProb

  resV <- .PredictAll(
    ## survival RF model stored in "result" is passed as an input
    ## has info containing list of individual trees (trees), aggregated results from forest (forest), variables,
    ## mTry,nCat,levels for each categorical variable (xLevels)
    object = result,

    ## newdata corresponds to eligible cases to be used for predictions
    ## data =  data.frame object containing covariate and treatment histories
    ## this is the same data that we grew the forest on

    newdata = data[elig,],
    params = params,

    ## using the input model, trtmt name, trtmtn levels avaiable at the current stage
    model = mod,
    txName = txName,
    txLevels = txLevels
  )


  ## this .PredictAll() sets all cases to receive first treatment, then uses input formula and newdata to make predictions
  ## uses .Predict()-- class_SurvRF.R line 176 to output: USING THE INPUT "newdata" of the eligible patients
  ## makes predictions based on first tree using .predictSurvTree():
  ## the tree has already been fit on the data, so we find the number of nodes in the tree
  ## traverse through the tree with the patients using the trained tree cutoffs to determined if they go left or right
  ## also, get the predicted mean and survival probability for each patient based on the terminal node they belong to
  ## predMean and predSurvProb
  ## then loops through the rest of the trees
  ## then for each tree, we add the predicted results form the survival function, mean, survival probability
  ## we divide by the total number of trees to get predicted values for each patient
  ## OUTPUT: forest-averaged predictions of survival function, mean survival, survival probability
  ## as part of output, get the survival function and convert that into a list
  ## iterate through rest of the treatment levels & set all cases to receive that treatment, use input formula and newdata to predic
  ## uses .Predict() to output:
  ## using the same method, OUTPUT: forest-averaged predictions of survival function, mean survival, survival probability
  ## as part of the output, get the survival function for the i-th treatment level
  ## also get the mean and survival probabilities by binding them in columns:
  ## so, we have mean from first treatment, mean frm second treatment, etc. in columns
  ## we also do this for survival probability

  ## next, we use .optimal() to calculate the optimal treatment based on the predictions (output = opt)
  ## depends on what the criteria for optimizing is (mean survival prob), mean survival prob and then survival probability, survival probability
  ## identify which element of what you're looking at (mean survival prob, surv prob) and return the prediction w/ maximum critieria

  ## we return a list of the predicted mean, survival probability, survival function, and the optimal treatment (output = predicted)


  ## create a new object of class"DTRSurvStep"

  result <- new(
    Class = "DTRSurvStep",
    "txName" = txName,
    "txLevels" = txLevels,
    "model" = mod,
    "survRF" = result,
    "eligibility" = elig,

    ## predicted values for all treatment levels
    "valueAllTx" = resV$predicted,

    ## optimal treatment recommendations extracted
    "optimal" = resV$optimal
  )

  ## return newly created "DTRSurvStep" object

  return(result)

}






















