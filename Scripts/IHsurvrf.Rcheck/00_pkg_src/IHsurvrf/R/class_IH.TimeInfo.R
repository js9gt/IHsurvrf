

# Class to store information regarding time points
#
# Class is not exported and is for internal convenience only
#
#  @slot timePoints A numeric vector; the timepoints used in the analysis
#
#  @slot timeDiff A numeric vector; the time differences between the timepoints
#
#  @slot tau A numeric object; maximum time
#
# Getters
#   .TimePoints(object, ...) {new; defined}
#   .NTimes(object, ...) {new; defined}
#   .Tau(object, ...) {new; defined}
#   .TimeDiff(object, ...) {new; defined}
#
# Methods
#   .TimeInfoAsList(object, ...) {new; defined}
#
# Functions
#   .timeInfo(timePoints, nTimes, response)
#

source("R/IH.VerifyTimePoints.R")

## define a new S4 class called "TimeInfo"
setClass(Class = "TimeInfo",

         ## A numeric vector; the timepoints used in the analysis
         slots = c(timePoints = "numeric",

                   ## timeDiff A numeric vector; the time differences between the timepoints
                   timeDiff = "numeric",

                   ## A numeric object; maximum time
                   tau = "numeric"))


#-------------------------------------------------------------------------------
# Method to retrieve the difference in timepoints
#-------------------------------------------------------------------------------
# Method returns the vector of the differences in timepoints used for the analysis
#-------------------------------------------------------------------------------

## create a new generic function called .TimeDiff

setGeneric(name = ".TimeDiff",
           def = function(object, ...) { standardGeneric(".TimeDiff") })

## establish a new method for .TimeDiff() that prevents function from being used to object without specific class

setMethod(f = ".TimeDiff",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })


## establish a new method for .TimeDiff working on objects of class Timeinfo

setMethod(f = ".TimeDiff",
          signature = c(object = "TimeInfo"),

          ## returns the value from the timeDiff slot of the "TimeInfo" object

          definition = function(object, ...) { return( object@timeDiff ) })


#-------------------------------------------------------------------------------
# Method to retrieve the maximum timepoint in class_IH.Parameters
# used in .criticalValue function
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------

## create new generic function called .Tau

setGeneric(name = ".Tau",
           def = function(object, ...) { standardGeneric(".Tau") })


## establish a new method for .Tau() that prevents function from being used to object without specific class
setMethod(f = ".Tau",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## establish a new method for .Tau working on objects of class TimeInfo
setMethod(f = ".Tau",
          signature = c(object = "TimeInfo"),

          ## returns the value from the "tau" slot of the "TimeInfo" object

          definition = function(object, ...) { return( object@tau ) })


#-------------------------------------------------------------------------------
# Method to retrieve timepoints in class_IH.Parameters
# used in .criticalValue function
#-------------------------------------------------------------------------------
# Method returns the vector of timepoints used for the analysis
#-------------------------------------------------------------------------------

## define generic function called .TimePoints()

setGeneric(name = ".TimePoints",
           def = function(object, ...) { standardGeneric(".TimePoints") })


## establish a new method for .TimePoints() that prevents function from being used to object without specific class

setMethod(f = ".TimePoints",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## implement a new method for .TimePoints()  applying to objects of TimeInfo class
setMethod(f = ".TimePoints",
          signature = c(object = "TimeInfo"),

          ## returns the timePoints slot of the object

          definition = function(object, ...) { return( object@timePoints ) })



#-------------------------------------------------------------------------------
# Method to retrieve the number of timepoints
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------

## create a new generic function called .NTimes()
## used in class_IH.DTRSurvStep.R

setGeneric(name = ".NTimes",
           def = function(object, ...) { standardGeneric(".NTimes") })

## create a method for .NTimes function that runs on ANY class to ensure .NTimes is only used with objects of defined classes
setMethod(f = ".NTimes",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })


## create a method for .NTimes that works on objects of class "TimeInfo"
setMethod(f = ".NTimes",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) {

            ## returns the length of the vector in timePoints slot
            ## AKA retrieves the number of time points
            return( length(x = object@timePoints) )
          })


#-------------------------------------------------------------------------------
# Function to verify inputs and create a TimeInfo object
#-------------------------------------------------------------------------------
# Function returns a TimeInfo object
#-------------------------------------------------------------------------------


## defines a new function called .timeInfo that creates & returns a "TimeInfo" object

.timeInfo <- function(timePoints, nTimes, tau, response) {

  ## defined in.VerifyTimePoints.R


  # ensure that timePoints and nTimes are appropriate. Methods return a vector
  # of unique time points that are sorted in ascending order.
  ## return a list containing the timePoints and tau
  timePoints <- .VerifyTimePoints(timePoints = timePoints,
                                  tau = tau,
                                  nTimes = nTimes,
                                  response = response)

  ## update tau and timePoints based on the output from .VertifyTimePoints

  tau <- timePoints$tau
  timePoints <- timePoints$timePoints

  # the total number of times points
  nTimes <- length(x = timePoints)

  # deltaT; T_{i} - T_{i-1} i = 1:nTimes with T_0 = 0
  # timeDiff <- c(timePoints[-1L] - timePoints[-nTimes], 0)
  # this is used for the truncated mean calculation in which it is
  # assumed that all t > tau equals tau to t(nTimes+1)-t(nTimes) = 0

  ## computes differences between consecutive time points to capture their intervals

  timeDiff <- timePoints[-1L] - timePoints[-nTimes]

  ## append difference after the last time point to 0
  timeDiff[nTimes] <- 0.0

  ## construct & return a new timeInfo object

  return( new(Class = "TimeInfo",

              ## assign tau to the "tau" slot
              "tau" = tau,

              ## assign timePoints to the "timePoints" slot
              "timePoints" = timePoints,

              ## assign calculated differences to "timeDiff" slot
              "timeDiff" = timeDiff) )

}

