
# Virtual class to denote objects are arising from survRF step
#
# Methods
#   .Predict(object, newdata, ...) {new; not allowed}

source("R/IH.predictSurvTree.R")
#source("~/survrf/Scripts/IHsurvrf/R/IH.predictSurvTree.R")
source("R/class_IH.TimeInfo.R")
#source("~/survrf/Scripts/IHsurvrf/R/class_IH.TimeInfo.R")



## create a new, virtual class called "SurvRFObject"
## since it's virtual, you cannot create objects of this class, used as a base class from which other classes inherit
setClass(Class = "SurvRFObject",
         contains = c("VIRTUAL"))




#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values (in preparation for defining .Predict on dif object classes)
#-------------------------------------------------------------------------------

## creates a generic function called .Predict()

setGeneric(name = ".Predict",
           def = function(object, newdata, ...) { standardGeneric(".Predict") })


## creates a method for .Predict() for objects of ANY class and new data of ANY class
## stops execution w/ an error to prevent function from being used without specific class implementation

setMethod(f = ".Predict",
          signature = c(object = "ANY",
                        newdata = "ANY"),
          definition = function(object, newdata, ...) { stop("not allowed") })

#-------------------------------------------------------------------------------
# method to make predictions for new data at each tx level (in preparation for defining .PredictAll on dif object classes)
#-------------------------------------------------------------------------------

## create a genetic function called .PredictAll()

setGeneric(name = ".PredictAll",
           def = function(object, ...) {
             standardGeneric(".PredictAll")
           })

## executes method for objects of ANY class
## stops execution w/ an error to prevent function from being used without specific class implementation

setMethod(f = ".PredictAll",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

          ## ------------------------------------------------ ##
          ##   For SurvRF objects: .Predict and .PredictAll   ##
          ## ------------------------------------------------ ##

## defines a new S4 class called "SurvRF"

setClass(Class = "SurvRF",
         slots = c(

           ## A list object. The results of the tree building algorithm for
           ##    each tree in the forest
           "trees" = "list",

           ## A list object. The values averaged across all trees in the
           ##    forest
           "forest" = "list",

           ## A character vector. The variables considered in the
           ##    analysis
           "variables" = "character",

           ## An integer. The maximum number of covariates considered for
           #    splitting
           "mTry" = "integer",

           ## An integer vector. The number of categories for each covariate
           #    considered. >=2 unordered factor, 1 ordered factor, 0 continuous
           "nCat" = "integer",

           ## A list object. The categories in each covariate considered.
           "xLevels" = "list"),

         ## SurvRF will inherit from SurvRFObject (virtual class)

         contains = c("SurvRFObject"))


#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values: for survRF object
#-------------------------------------------------------------------------------
# method with NULL retrieves fitted values from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing survFunc, mean, and? survProb
#-------------------------------------------------------------------------------

## creates a new method for the .Predict() function for the SurvRF class when the newdata argument is NULL
## meant to retrieve fitted values from a SurvRF object instead of making predictions on new data

setMethod(f = ".Predict",
          signature = c(object = "SurvRF",
                        newdata = NULL),

          ## the code that's executed when the function is called
          ## method returns the contents of the "forest" slot from the SurvRF object
          ## this object contains list object containing survFunc, mean, and? survProb

          definition = function(object, newdata, ...) {
            return( object@forest )
          })


#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculates value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing survFunc, mean, and? survProb
#-------------------------------------------------------------------------------


## creates a new method for the SurfRF class when the newdata argument is a dataframe

setMethod(f = ".Predict",
          signature = c(object = "SurvRF",
                        newdata = "data.frame"),
          definition = function(object,
                                newdata,
                                ...,
                                params) {


            # verify that there are not new levels in the the data
            # this assumes that newdata has been passed in with
            # covariates in the order used in the analysis. This is
            # guaranteed if the data.frame is created from the model.

            ## extract levels of all factors in "newdata" to make sure the categorical variables match those used for model training
            xLevels <- lapply(X = newdata, FUN = levels)

            #browser()

            ## iterate over each variable

            for (i in length(x = xLevels)) {

              ## if both the newdata level for ith variable and model level for ith variable are null, skip current iteration

              if (is.null(x = xLevels[[ i ]]) &&
                  is.null(x = object@xLevels[[ i ]])) next

              ## if data introduces new factor levels not in model, stop execution and present error message

              if (any(! {xLevels[[ i ]] %in% object@xLevels[[ i ]]})) {
                stop("new factor levels not present in the training data",
                     call. = FALSE)

                #browser()
              }
            }

            # verify that type of data is the same as the training data
            # type means numeric (nCat = 0), ordered factor (nCat = 1), or
            # factor (nCat = length(levels(x)))

            ## determines number of categories for each variable in newdata

            nCat <- sapply(X = xLevels, FUN = length)

            ## if ncat is ordered factors, adjust this value to 1

            nCat <- ifelse(test = sapply(X = newdata, FUN = is.ordered),
                           yes = 1L,
                           no = nCat)

            ## if the datatypes in model don't match the data types in nCat,
            # object@nCat is from the input survRF object

            if (any(unlist(x = object@nCat) != unlist(x = nCat))) {

              ## stop running and output warning message

              stop("type of predictors in newdata do not match the training data",
                   call. = FALSE)
            }

            ## retrieves number of trees in the forest from the input survRF object

            nTree <- length(x = object@trees)

            ## predict value for the first tree

            # predict for first tree
            # .predictSurvTree() is an internal function defined in predictSurvTree.R

            #browser()

            newResult <- .predictSurvTree(x = newdata,
                                          params = params,
                                          nCat = nCat,

                                          ## this is how you know the prediction is for the first tree
                                          ## the trees slot contains a list of all the trees in the forest


                                          nodes = object@trees[[ 1L ]])

            #browser()

            ## create counter to set tree index to 2
            ## loop through the rest of the trees

            iTree <- 2L
            while (iTree <= nTree) {

              ## for each tree, predict for each tree using .predictSurvTree


              # predict for tree iTree; sum result
              # .predictSurvTree() is an internal function defined in predictSurvTree.R

              tmp <- .predictSurvTree(x = newdata,
                                      params = params,
                                      nCat = nCat,
                                      nodes = object@trees[[ iTree ]])

              ## ensemble averaging

              ## add results from the first tree with the predicted value from the ith tree
              ## AKA across all tree, take the cumulative sum of these:
              ##                 survFunc
              ##                 mean
              ##                 survProb (if not null)

              newResult$survFunc <- newResult$survFunc + tmp$survFunc
              newResult$mean <- newResult$mean + tmp$mean

              ## if the survProb isn't null, also take the cumulative sum


              if (!is.null(x = newResult$survProb)) {
                newResult$survProb <- newResult$survProb + tmp$survProb
              }

              ## increase tree counter by 1

              iTree <- iTree + 1L
            }

            ## averages aggregated predictions by the total number of trees to get final prediction values

            newResult$survFunc <- newResult$survFunc / nTree
            newResult$mean <- newResult$mean / nTree
            if (!is.null(x = newResult$survProb)) {
              newResult$survProb <- newResult$survProb / nTree
            }

            ## return averaged predictions of survFunc, mean, and survProb

            return( newResult )

          })


#-------------------------------------------------------------------------------
# method to make predictions for new data for SurvRF object (pooled analysis) from class_DTRSurvStep.R
#-------------------------------------------------------------------------------
# method with data.frame calculates value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing survFunc, mean, and? survProb
#-------------------------------------------------------------------------------

## this code defines a method for the .PredictAll() function

setMethod(f = ".PredictAll",

          ## used for object class SurvRF

          signature = c(object = "SurvRF"),
          definition = function(object, ..., newdata, model, txLevels, txName, params) {

            ## initially sets all cases in newdata to receive the first treatment level

            # set all cases to receive first treatment
            newdata[[ txName ]][] <- txLevels[1L]

            ## creates a model frame "x" that uses the input formula in model & the new data
            ## this will be used to make predictions

            # extract new model frame
            x <- stats::model.frame(formula = model, data = newdata)

            ## if the model formula includes a response variable, remove this and only leave the predictors

            # remove response from x
            if (attr(x = terms(x = model), which = "response") == 1L) {
              x <- x[,-1L,drop = FALSE]
            }

            ## return averaged predictions of survFunc, mean, and survProb over all trees

            # calculate the estimated values for this treatment level (first)
            # .Predict() is a method; called here for objects of class SurvRF, which
            # is defined in this file-- uses file .predictSurvTree.R (and fortran)
            ## this basically returns tree-averaged results (survival func, mean survival time, surv prob
            res <- .Predict(object = object, newdata = x, params = params, ...)

            ## convert the output "survFunc" into a list

            res$survFunc <- list(res$survFunc)

            ## set index of treatment level to 2
            ## iterate through the rest of the treatment levels

            i <- 2L

            # repeat this process for each tx level
            while (i <= length(x = txLevels)) {

              ## update newdata to set the treatment to the current level

              # set all cases to receive the ith treatment
              newdata[[ txName ]][] <- txLevels[i]

              # extract new model frame
              x <- stats::model.frame(formula = model, data = newdata)

              # remove response from x
              if (attr(x = terms(x = model), which = "response") == 1L) {
                x <- x[,-1L,drop=FALSE]
              }

              # calculate the estimated values for this treatment level
              # .Predict() is a method; called here for objects of class SurvRF, which
              # is defined in file class_SurvRF.R
              tt <-  .Predict(object = object,
                              newdata = x,
                              params = params, ...)

              ## stores suvFunc (averaged over all trees)

              res[[ "survFunc" ]][[ i ]] <- tt$survFunc

              ## append mean survival time and survival probability (averaged over all trees) predictions using cbind
              ## this means we combine these quantities side by side for each treatment level
              ## res$mean: contains mean survival times for observations in the new data for each treatment level AGGREGATED
              ## tt$mean: represents new set of mean times for each treatment

              res[[ "mean" ]] <- cbind(res$mean, tt$mean)
              res[[ "survProb" ]] <- cbind(res$survProb, tt$survProb)

              ## increase index

              i <- i + 1L
            }

            ## calculate the optimal treatment based on the predictions
            ## defined later in this script that returns the optimal treatment, outcome, etc.

            opt <- .optimal(params = params, predicted = res, txLevels = txLevels)

            ## return a list containing the predictions for each treatment level and the optimal treatment decision
            ## also containing survFunc, mean, and? survProb

            return( list("predicted" = res, "optimal" = opt) )
          })






#-------------------------------------------------------------------------------
# Calculating optimal results
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# method returns the element containing the maximum of the critical value criterion
#-------------------------------------------------------------------------------


## initializes a new function called .optimal

.optimal <- function(params,
                     ## the results from the tree-averaged predictions
                     predicted, txLevels) {

  ## uses .CriticalValueCriterion which is initialized in class_CriticalValue.R as a generic
  ## also sued in class_CriticalValueMean.R and class_CriticalValueSurvival.R
  ## this returns "mean" for CriticalValueMean objects, and surv.mean for CriticalValueSurvival objects of type mean
  ##                    also returns surv.prob for CriticalValueSurvival objects of type prob

  crit <- .CriticalValueCriterion(params)

  ## if the criterion returned is mean,

  if (crit == "mean") {

    ## iterates over each row (observation) in predicted$mean and finds the treatment level maximizing the mean survival time for each observation

    # identify which element contains the maximum expected survival time
    ## if there's only one treatment, only use that treatment
    ## NOTE: the if else: statement is new
    if (is.vector(predicted$mean) == TRUE) {
      ## turn optimal into a column
      opt_predicted <- matrix(predicted$mean, nrow = length(predicted$mean), ncol = 1)

      optTx <- apply(X = opt_predicted,
            MARGIN = 1L,
            FUN = .whichMax,
            tieMethod = params@tieMethod)


    } else {
      optTx <- apply(X = predicted$mean,
                     MARGIN = 1L,
                     FUN = .whichMax,
                     tieMethod = params@tieMethod)
    }

    ## if the criterion returned is surv.mean,

  } else if (crit == "surv.mean") {
    # identify which element
    # first, mean survvial probability; if ties, then use mean survival time to break tie


    ## iterates over each row (observation) in predicted$mean and finds the treatment level maximizing the mean survival time for each observation

    ## specifies how to break ties when multiple treatemnts result in the same maximum mean survival time
    ## if there are ties, the tie breaking method uses the input tieMethod for "params"


    optTx <- apply(X = predicted$mean,
                   MARGIN = 1L,
                   FUN = .whichMax,
                   tieMethod = params@tieMethod)


    # index for which the survival probability is max
    ## first, look at mean survival probability

    tmp <- apply(X = predicted$survProb,
                 MARGIN = 1L,

                 ## defined at the bottom

                 FUN = .whichMax,

                 ## if there's a tie in survival probability, this will contain NA
                 tieMethod = "NA")

    # for those that are not tied in survival probability, replace mean survival time
    # with mean survival probability
    ## giving survival probability precedence in resolving ambiguities or ties in mean survival time

    ## it tmp not NA means there's no tie (there's a clear maximumO)
    isna <- is.na(x = tmp)

    ## updates the index of optTx (best mean survival times) for values where there IS a tie
    ## then, look at the mean survival times that are maximum
    optTx[!isna] <- tmp[!isna]


    ## if the criterion is surv.prob,

  } else if (crit == "surv.prob") {

    ## iterates over each row (observation) in predicted$survProb and finds the treatment level maximizing the expected survival probability for each observation


    # identify which element contains the maximum expected survival probability
    optTx <- apply(X = predicted$survProb,
                   MARGIN = 1L,

                   ## defined at the bottom

                   FUN = .whichMax,
                   tieMethod = params@tieMethod)
  }

  ## initializes a matrix "optSv" to store survival functions that correspond to the optimal treatment for each observation

  # extract the survival function at optimal tx
  ## number of columns same as number of timepoints, using .NTimes from class_IH.TimeInfo.R
  optSv <- matrix(data = 0.0, nrow = length(x = optTx), ncol = .NTimes(params))

  ## iterates through each observation and gets the survival function from the optimal treatment index
  ## each patient will have an optimal treatment, so we basically loop through each patient's optimal treatment index
  ## take the survival functions from those

  for (i in 1L:length(optTx)) {
    optSv[i,] <- predicted$survFunc[[optTx[i]]][,i]
  }

  ## returns an object of class "Optimal"
  ## which outputs the optimal treatment level for each observation, the survival fnctions at these optimal treatments, and the optimization criterion

  return( new(Class = "Optimal",
              "optimalTx" = txLevels[optTx],
              "optimalY" = optSv,
              "type" = crit) )
}



#-------------------------------------------------------------------------------
# internal function to identify the maximum value of input vector x
#  @param x a vector of values
#  @param tieMethod a character indicating the method to be used to
#    breaks {first, random, NA}
#-------------------------------------------------------------------------------
# function returns a single numeric or NA
#-------------------------------------------------------------------------------

## defines an internal function .whichMax to identify index of the maximum value in vector "x"
.whichMax <- function(x, tieMethod) {

  # Check if x has only one element. If so, just return the input
  if (length(x) == 1) return(1)

  ## identifies indices of values in x within a very small range of the maximum value
  ## treats values as equal to the maximum if they're very close

  ind <- which(x >= {max(x) - 1e-8})

  ## if there's exactly one index identified, function returns this index directly

  if (length(x = ind) == 1L) return( ind )

  ## if tieMethod == first, return the first index among the tied maximum values (prioritize earliest occurrence)

  if (tieMethod == "first") return( ind[1L] )

  ## if tieMethod == random, randomly select an index from those tied for the maximum value
  ## resample defined below

  if (tieMethod == "random") return( resample(x = ind, size = 1L) )

  ## if neither tie method are specfied or tie method isn't one of these, return an NA

  return( NA )
}


## input: a vector x and an optional set of arguments (...)
## resamples the elements of x using sample.int to generate sample indices from 1:x
## uses indices to reorder the elements in x
resample <- function(x, ...) x[sample.int(n = length(x = x), ...)]
