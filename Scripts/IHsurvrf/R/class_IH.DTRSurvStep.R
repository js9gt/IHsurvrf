


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
#source("~/survrf/Scripts/IHsurvrf/R/IH.survRF.R")
source("R/class_IH.Optimal.R")
#source("~/survrf/Scripts/IHsurvrf/R/class_IH.Optimal.R")
source("R/IH.shiftMat.R")
#source("~/survrf/Scripts/IHsurvrf/R/IH.shiftMat.R")
## For policy evaluation:
#source("R/class_IH.SurvRF.R")

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
    "optimal" = "Optimal",

    ## stores each patient's shifted probabilities from appending
    "stageappend" = "matrix",

    ## a logical statement telling us we have enough patients in the current stage (at least 5)
    "sumElig" = "logical",

    ## TRUE means we only have 1 treatment level present-- we need at least 2
    "enoughtrt" = "logical"
  )
)

# generic defined in class_Optimal.R
setMethod(f = ".OptimalAsList",
          signature = c(object = "DTRSurvStep"),
          definition = function(object, ...) {
            return( .OptimalAsList(object = object@optimal) )
          })


#-------------------------------------------------------------------------------
# method to retrieve predicted values: used in policy evaluation step
#-------------------------------------------------------------------------------
# if findOptimal is TRUE, method stops with error
# if findOptimal is FALSE, method returns a Value object
#-------------------------------------------------------------------------------

## .Predict() on these object classes (survRF or survRFStratified) are defined in class_IHSurvRF.R

setMethod(f = ".Predict",
          signature = c(object = "DTRSurvStep",
                        newdata = NULL),
          definition = function(object, newdata, ...) {

            ## the code that's executed when the function is called
            ## method returns the contents of the "forest" slot from the SurvRF object

            return( .Predict(object = object@survRF, newdata = NULL, ...) )

          })

#-------------------------------------------------------------------------------
# method to predict value for new data
#-------------------------------------------------------------------------------
# if findOptimal is TRUE, method returns a list containing a Value object and
#   an Optimal object
# if findOptimal is FALSE, method returns a Value object
#-------------------------------------------------------------------------------


setMethod(f = ".Predict",
          signature = c(object = "DTRSurvStep",
                        newdata = "data.frame"),
          definition = function(object,
                                newdata,
                                ...,
                                params,
                                findOptimal) {

            # update model to remove response
            mod <- update(object@model, NULL ~ .)


          # ensure data contains all model covariates
          x <- tryCatch(expr = stats::model.frame(formula = mod,
                                                  data = newdata),
                        error = function(e) {
                          stop("variables in the training data missing in newdata",
                               call. = FALSE)
                        })

          # remove response from x
          # this should no longer happen, but keeping anyway
          if (attr(x = terms(x = mod), which = "response") == 1L) {
            x <- x[,-1L,drop = FALSE]
          }




          ## if findOptimal = TRUE, return the optimal predictions

           if (findOptimal) {
             # if optimal is requested make predictions for all possible
             # treatment options-- go through each treatment ant see what the patient's estimated survival would have been

             resV <- .PredictAll(object = object@survRF,
                                 newdata = newdata,
                                 params = params,
                                 model = mod,
                                 txName = object@txName,
                                 txLevels = object@txLevels)



             ## return a list containing the predictions for each treatment level and the optimal treatment decision

             return( list("value" = resV$predicted, "optimal" = resV$optimal) )

           } else {


             ## otherwise, calculate the estimatated survival values that the patient receives using the actions they already got
             return( .Predict(object = object@survRF,
                              newdata = x,
                              params = params, ...) )
           }          })


# Define a new function for DTRSurvStep
PredDTRSurvStep <- function(object, newdata, ..., params, findOptimal) {

  # Ensure object is of the correct class
  if (!inherits(object, "DTRSurvStep")) {
    stop("object must be of class 'DTRSurvStep'")
  }

  # Check if newdata is of the correct type
  if (!is.data.frame(newdata) && !is.null(newdata)) {
    stop("newdata must be a data.frame or NULL")
  }

  # Code for handling NULL newdata case
  if (is.null(newdata)) {
    return(.Predict(object = object@survRF, newdata = NULL, ...))
  }

  # Update model to remove response
  mod <- update(object@model, NULL ~ .)

  # Ensure data contains all model covariates
  x <- tryCatch(expr = stats::model.frame(formula = mod,
                                          data = newdata),
                error = function(e) {
                  stop("variables in the training data missing in newdata",
                       call. = FALSE)
                })

  # Remove response from x if necessary
  if (attr(x = terms(x = mod), which = "response") == 1L) {
    x <- x[,-1L,drop = FALSE]
  }

  # If findOptimal is TRUE, return the optimal predictions
  if (findOptimal) {
    resV <- .PredictAll(object = object@survRF,
                        newdata = newdata,
                        params = params,
                        model = mod,
                        txName = object@txName,
                        txLevels = object@txLevels)
    return(list("value" = resV$predicted, "optimal" = resV$optimal))
  } else {
    return(.Predict(object = object@survRF,
                    newdata = x,
                    params = params, ...))
  }
}



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
# @params: pool1: A logical which defaults to FALSE. FALSE meaning we are not in the first step of pooling data for the overall RF after fitting HC method


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
                         sampleSize,
                         pool1 = FALSE,
                         appendstep1 = FALSE,
                         inputpr = NULL,
                         input_prop
                         ) {
  # model is input
  mod <- model

  ## read in the propensity score
  prop <-input_prop

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
  ## only include the variables that are used in the model formula mod
  ## If there are variables in your data that are not used in the formula, they won't be included in the resulting x.

  x <- stats::model.frame(formula = mod,
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

  # Check if sum(elig) < 3
#  if (sum(elig) < 5) {
#    return(new(
#      Class = "DTRSurvStep",
#      "txName" = "A",
#      "txLevels" = c(0, 1),
#      "model" = mod,
#      "survRF" = NA,
#      "eligibility" = T,
#
#      ## predicted values for all treatment levels
#      "valueAllTx" = list(0, 0),
#
#      ## optimal treatment recommendations extracted
#      "optimal" = new("Optimal",
#                      optimalTx = c(0),  # Initial value for optimalTx slot
#                      optimalY = matrix(0, nrow = 2, ncol = 2),  # Initial value for optimalY slot
#                      type = "mean"),
#
#      ## adding a new slot to retrieve the appending probabilities for each stage
#      "stageappend" = matrix(0, nrow = 2, ncol = 2),
#
#      ## logical if we have enough patients in the stage to move on
#      ## TRUE means we don't
#      "sumElig" = TRUE
#
#    ))
#  }

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


  # If appendstep1 is FALSE, execute the code block
  if (!appendstep1) {
    if (is.null(x = priorStep)) {
      # PriorStep is NULL for the first step of the analysis.
      # Transform the time variable to a probability mass vector

      ## applies function FUN to each element of X which checks if each survival time (s) is less than each time point
      ## .TimePoints() part specifies an additional argument to be passed to fun
      ## .TimePoints() defined in class_TimeInfo.R

      # Identify time points <= response
      tSurv <- sapply(
        X = response[elig],
        # Results in a logical integer (0 = FALSE)
        FUN = function(s, tp) {
          as.integer(x = {
            s < tp
          })
        },
        tp = .TimePoints(object = params)
      )

      # Time point nearest the response without going over
      # {nTimes x nElig}
      pr <- rbind(tSurv[-1L,], 1) - tSurv

    } } else {

#      {
#      # PriorStep is not NULL when q < Q
#      ### priorStep is of class DTRSurvStep AKA the results from the prior step of the analysis
#      #
#      # the number of timepoints
#      # .NTimes() is a getter method defined for Parameters objects
#
#      ## retrieve the number of timepoints from the params object
#      ## defined in class_TimeInfo.R
#
#      nTimes <- .NTimes(object = params)
#
#      # create an empty matrix for all previously eligible cases
#
#      ## initializes a survival matrix called "survMatrix" with 0's
#
#      survMatrix <- matrix(data = 0.0,
#                           nrow = nTimes,
#                           ncol = nrow(x = x))
#
#      ## sets the first row to 1.0
#
#      survMatrix[1L,] <- 1.0
#
#      ## updates survMatrix with survival functions estimated from the previous step
#      ## priorStep: A DTRSurvStep object. The analysis from a previous step
#
#      # retrieve estimated OPTIMAL survival function from previous step
#      ## then accesses the "eligibility" slot (logical) to select only the columns where the patient is still eligible
#      ## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
#      ## transpose the output of .OptimalY to align with the structure: AKA number of rows = timepoints
#      ## when it’s retrieved with .OptimalY shows up as rows = pts, and columns = timepoints, so we just transpose it back to its normal shape
#      ## .OptimalY defined in class_Optimal.R
#      ## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
#      ## optimal slot is  object of class optimal
#      ## .OptimalY acts to return the optimalY slot
#
#      survMatrix[, priorStep@eligibility] <-
#        t(.OptimalY(object = priorStep@optimal))
#
#      # shift the survival function down in time (T_i - Tq) and
#      # transform to a probability mass vector for only those
#      # eligible for this stage
#      # .shiftMat is an internal function defined in shiftMat.R
#
#      ## transforms the survival functions in survMatrix to probability mass format for the current stage
#      ## shifts the survival function based on the observed survival times
#      ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage
#
#      pr <- .shiftMat(
#        timePoints = .TimePoints(object = params),
#
#        ## extracts columns from survMatrix corresponding to cases that are eligible
#        ## this is a matrix matrix where each column represents survival function for an individual
#        survMatrix = survMatrix[, elig, drop = FALSE],
#
#        ## extracts survival times corresponding to eligible cases
#        ## this is how much to shift survival function for each individual
#        shiftVector = response[elig],
#
#        ## probably transforming survival times into probabilities?
#        surv2prob = TRUE
#      )
#
#      ## sets very small values in pr to 0
#
#      pr[abs(pr) < 1e-8] <- 0.0
#
#    }
#  } else
    # Use the input value of "pr" for further calculations

    pr <- inputpr

    pr[abs(pr) < 1e-8] <- 0.0

    if (is.null(inputpr)) {
      stop("Input value for 'pr' is required when 'appendstep1' is TRUE.")
    }
  }



  ## if there are any NA values in the matrix, stop function execution and display error message
  if (any(is.na(x = pr)))
    stop("NA not permitted in pr -- contact maintainer",
         call. = FALSE)

  ## checks if any values in pr are outside the range [0,1]. If so, stops with error message
  ## pr is expected to be a matrix of probabilities and therefore between 0 and 1

  if (any(pr > 1.0) || any(pr < 0.0)) {
    # Print values that are out of bounds
    print(pr[pr > 1.0 | pr < 0.0])

    # Clip values outside the range [0, 1]
    pr[pr > 1.0] <- 1.0
    pr[pr < 0.0] <- 0.0

    # Optional: Print a message indicating values were clipped
    message("Some pr values were outside the range [0, 1] and have been clipped.")
  }


  ## if txName variable is a factor,

  # identify tx levels in limited data
  if (is.factor(x = data[, txName]) & !pool1) {
    ## extract levels of the factor treatment variable for eligible cases
    txLevels <- levels(x = factor(x = data[elig, txName]))
  } else if ( !is.factor(x = data[, txName]) & !pool1) {  # Add condition to execute else block only if pool1 is FALSE

    ## if txName is not a factor and pool1 is FALSE, find unique values of the treatment for eligible cases & sort them to use as treatment levels
    txLevels <- sort(x = as.numeric(as.character(unlist(unique(data[elig, txName])))))
  }else  if (pool1) {  # If pool1 is TRUE
    ## if pool1 is TRUE, sort unique values of A.opt.HC for eligible cases
    txLevels <- sort(x = as.numeric(unlist(unique(data[elig, txName]))))

  }

  ## checks if there's only one treatment level in the data
  ## if so, display a message to the user-- this is important bc analyses often require multiple treatment levels for comparison

  if (length(x = txLevels) == 1L) {
    message("***only one treatment level in data***")

    ## also return early,with indicator for txLevel
#    return(new(
#      Class = "DTRSurvStep",
#      "txName" = "A",
#      "txLevels" = c(0, 1),
#      "model" = mod,
#      "survRF" = NA,
#      "eligibility" = T,
#
#      ## predicted values for all treatment levels
#      "valueAllTx" = list(0, 0),
#
#      ## optimal treatment recommendations extracted
#      "optimal" = new("Optimal",
#                      optimalTx = c(0),  # Initial value for optimalTx slot
#                      optimalY = matrix(0, nrow = 2, ncol = 2),  # Initial value for optimalY slot
#                      type = "mean"),
#
#      ## adding a new slot to retrieve the appending probabilities for each stage
#      "stageappend" = matrix(0, nrow = 2, ncol = 2),
#
#      ## logical if we have enough patients in the stage to move on
#      ## TRUE means we don't
#      "sumElig" = FALSE,
#
#      ## TRUE means we only have 1 treatment level present-- we need at least 2
#      "enoughtrt" = TRUE
#
#    ))

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
      sampleSize = sampleSize,

      prop = prop
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
    "optimal" = resV$optimal,

    ## adding a new slot to retrieve the appending probabilities for each stage
    "stageappend" = pr,

    ## logical if we have enough patients in the stage to move on
    "sumElig" = FALSE,

    ## TRUE means we only have 1 treatment level present-- we need at least 2
    "enoughtrt" = FALSE

  )

  ## return newly created "DTRSurvStep" object

  return(result)

}

#-------------------------------------------------------------------------------
# Internal function to find the mean values of maximum expected survival times across treatment levels
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## defines an internal function .meanValue
## returns mean values of expected survival times and survival probabilities
## used in dtrSurv.R

.meanValue <- function(object, ...) {
  ## initializes empty list
  res <- list()

  ## find the mean of maximum expected survival times across treatment levels
  # returns the mean of the expected survival times & adds them to the result in a new element "Et"
  res[[ "Et" ]] <-  mean(x = apply(

    ## access "mean" component of the valueAllTx slot which stores value of each tree for each treatment level
    X = object@valueAllTx$mean,

    ##apply the "max" function across the margin of the mean matrix/array to apply across rows
    ## AKA for each row (individual),take maximum value across all columns (treatment options), then calculate the mean
    MARGIN = 1L,
    FUN = max))

  ## if survival probabilities are available,
  ## checks the "survProb" comonent of the valueAllTx slot

  if (!is.null(object@valueAllTx$survProb)) {

    ## also calculate their means
    res[[ "St" ]] <-  mean(x = apply(X = object@valueAllTx$survProb,
                                     MARGIN = 1L,
                                     FUN = max))
  }

  ## return results as a list
  return( res )
}

























