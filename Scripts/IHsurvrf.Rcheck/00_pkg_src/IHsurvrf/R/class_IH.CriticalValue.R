

# Virtual class to store information regarding critical value selection
#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueAsList(object, ...) {not allowed}
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {not allowed}
#


## defining a virtual class called "CriticalValueBase"
## virtual class:
##                 class that isn't intended to be directly made, serves as a base class that other classes can inherit
##                 "child" classes can inherit properties from the parent class; useful for defining
##                         common structure or behavior that various derived classes can share


setClass(Class = "CriticalValueBase",
         contains = c("VIRTUAL"))


#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method is not defined for a general CriticalValue object
#-------------------------------------------------------------------------------

## creates generic method called ".CriticalValueCriterion"
## generic method:
##                 functions that behave differently based on the class they act upon
##  this allows method dispatch based on the class of the first argument to the function
##  allows the same function name to be used for different types of objects
##  the specifics for how the function behaves is determined by the class of object passed to it


setGeneric(name = ".CriticalValueCriterion",
           def = function(object, ...) { standardGeneric(".CriticalValueCriterion") })

## sets a method for the generic function ".CriticalValueCriterion"
## when the function is called with ANY object class, it should execute the function which outputs an error for general objects
## AKA irrespective of the class of object passed, the method defined for object = ANY will throw an error
## this catches errors for any object class where a more specific method is not defined


setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })


#-------------------------------------------------------------------------------
# CRITICAL VALUE MEAN
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


## defines a new class called "CriticalValueMean" which is a subclass of CriticalvalueBase
## meaning, the new class, CriticalValueMean, inherits properties from the CriticalValueBase class (class_CriticalValue.R)


setClass(Class = "CriticalValueMean",
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------


setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueMean"),
          definition = function(object, ...) { return( "mean" ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueMean
#-------------------------------------------------------------------------------
# function returns a CriticalValueMean object
#-------------------------------------------------------------------------------


## function creates and returns a new object of class "CriticalValueMean"
## uses object-oriented programming
.criticalValueMean <- function(...) { return( new("CriticalValueMean") ) }

#-------------------------------------------------------------------------------
# CRITICAL VALUE SURVIVAL
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------



## defines a new class called "CriticalValueSurvival" which is a subclass of CriticalValueBase
## AKA, inherits properties from CriticalValueBase


setClass(Class = "CriticalValueSurvival",

         ## includes 4 slots (variables/attributes) specific to this class
         slots = c(
           ## survivalTime of any type of data
           survivalTime = "ANY",
           ## integer
           sIndex = "integer",
           ## numeric
           sFraction = "numeric",
           ## character
           type = "character"),
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------

setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueSurvival"),
          definition = function(object, ...) {
            if (object@type == "mean") return( "surv.mean" )
            if (object@type == "prob") return( "surv.prob" )
          })

setMethod(f = "initialize",
          signature = c(.Object = "CriticalValueSurvival"),
          def = function(.Object, ..., survivalTime, sIndex, sFraction, type) {

            obj <- list(...)
            tst <- sapply(X = obj,
                          FUN = function(x){
                            is(object = x,
                               class2 = "CriticalValueSurvival")
                          })

            if (any(tst)) {
              .Object <- obj[[ which(tst) ]]
            } else if (missing(x = survivalTime) && missing(x = sIndex) &&
                       missing(x = sFraction) && missing(x = type)) {
              .Object@survivalTime <- Inf
              .Object@sIndex <- -1L
              .Object@sFraction <- 0
              .Object@type <- "none"
            } else {
              if (missing(x = survivalTime) || missing(x = sIndex) ||
                  missing(x = sFraction) || missing(x = type)) {
                gn <- unlist(lapply(list(...),is))
                if( "CriticalValueSurvival" %in% gn ) return(.Object)
                stop("insufficient inputs provided")
              }
              if (type %in% c("surv.prob", "prob")) {
                .Object@type <- "prob"
              } else if (type %in% c("surv.mean", "mean")) {
                .Object@type <- "mean"
              }
              .Object@survivalTime <- survivalTime
              .Object@sIndex <- sIndex
              .Object@sFraction <- sFraction
            }
            return( .Object )
          })

#-------------------------------------------------------------------------------
# method to identify if critical value is of survival type
#-------------------------------------------------------------------------------
# method returns a logical
#-------------------------------------------------------------------------------

setGeneric(name = ".IsSurvival",
           def = function(object, ...) { standardGeneric(".IsSurvival") })


setMethod(f = ".IsSurvival",
          signature = c(object = "ANY"),
          definition = function(object, ...) { return( FALSE ) })


setMethod(f = ".IsSurvival",
          signature = c(object = "CriticalValueSurvival"),
          definition = function(object, ...) { return( TRUE ) })


#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueSurvival
#-------------------------------------------------------------------------------
# function returns a CriticalValueSurvival object
#-------------------------------------------------------------------------------


## defines a function to create & initiaize an object of class "CriticalValueSurvival"
## function inputs 3 parameters: survivalTime, timepoints, type
.criticalValueSurvival <- function(survivalTime, timePoints, type) {

  ## length of timePoints vector: total number of time points provided
  nTimes <- length(x = timePoints)

  ## sum of all time points that are less tham or equal to survival time
  ## AKA how many time points occur before or at survivalTime
  # index of last time point <= SurvivalTime
  sIndex <- sum(timePoints <= survivalTime)

  ## if the index < nTimes, means survivalTime falls between two points
  if (sIndex < nTimes) {
    # if SurvivalTime is below tau, determine the fraction
    sFraction <- {survivalTime - timePoints[sIndex]} /
      ## calculates difference between time point at indixes sIndex + 1 and sIndex
      {timePoints[sIndex + 1L] - timePoints[sIndex]}

    ## if the index is 0, survivalTime < the minimum time point,
  } else if (sIndex == 0L) {

    ## function stops with error message indicating extrapolation not possible
    # if it is below the minimum, stop -- cannot extrapolate
    stop("survival time is below minimum timepoint")
  } else {

    ## if the index is equal to nTimes, survivalTime is a or beyond the last time point
    ## set sFraction to 0, and set survivalTime to the maximum time point

    # if it is above tau, use tau -- cannot extrapolate
    sFraction <- 0.0
    survivalTime <- max(timePoints)

    ## message to inform adjustment
    message("survival time reset to tau")
  }

  ## function returns a new object of class "CriticalValueSurvival"
  ## object is initialized with calculated values for survivalTime, sIndex, sFraction, and the input type
  ## output class is the same one that was created earlier in the script
  return( new("CriticalValueSurvival",
              "survivalTime" = survivalTime,
              "sIndex" = sIndex,
              "sFraction" = sFraction,
              "type" = type) )

}


