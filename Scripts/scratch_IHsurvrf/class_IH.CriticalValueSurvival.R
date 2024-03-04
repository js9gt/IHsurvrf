
# Class extends CriticalValue to indicate that critical value is survival
#
# Class is not exported and is for internal convenience only
#
# @slot SurvivalTime A numeric object. The time at which the survival 
#   probability is to be estimated
#
# @slot sIndex An integer object. The index of the timePoint vector above
#   which the survival time lies (and it is below the sIndex + 1 element)
#
# @slot sFraction A numeric object. The fractional location of SurvivalTime in
#   t[sIndex] and t[sIndex+1]
#
# @slot type A character object. Indicates of mean is to be used to break ties.
#
# Methods
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {defined}
#  .IsSurvival(object, ...) {new; defined}
#
# Functions
# .criticalValueSurvival(SurvivalTime, timePoints)
#


## defines a new class called "CriticalValueSurvival" which is a subclass of CriticalValueBase
## AKA, inherits properties from CriticalValueBase


#' @include class_IH.CriticalValue.R
source("class_IH.CriticalValue.R")





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