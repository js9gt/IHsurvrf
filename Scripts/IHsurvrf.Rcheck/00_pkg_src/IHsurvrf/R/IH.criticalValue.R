

# Function to verify inputs and create appropriate CriticalValue object
#
#  Function is not exported and for internal convenience only
#
# Function returns an object of class CriticalValueMean or CriticalValueSurvival
#


source("R/class_IH.CriticalValue.R")
source("R/IH.VerifyCriticalValue.R")
source("R/IH.VerifySurvivalTime.R")

## define a new function .criticalValue to verify input for critical values

.criticalValue <- function(criticalValue,
                           survivalTime,
                           tau,
                           timePoints) {

  ## in VerifyCriticalValue.R script to make sure input criticalValue is of the correct class(?)

  # ensure criticalValue is one of {'mean', 'surv.prob', 'surv.mean'}.
  # Methods return the original character possibly modified to be lower case.

  criticalValue <- .VerifyCriticalValue(criticalValue = criticalValue)

  ## in VerifySurvivalTime.R script to make sure that the survival time is provided when criticalValue is "surv.prob" or "surv.mean"
  ## also ensures survivalTime is within valid bounds and >= 0 and < tau


  # ensure that survivalTime is provided if criticalValue is one of
  # {'surv.prob', 'surv.man'}. Methods return the numeric survivalTime or NULL.

  ## used in VerifySurvivalTime.R

  survivalTime <- .VerifySurvivalTime(survivalTime = survivalTime,
                                      criticalValue = criticalValue,
                                      tau = tau)

  ## if survivalTime is specified (not null)

  if (!is.null(x = survivalTime)) {
    # if survivalTime is given as input, verify it and the timePoints input and
    # create a CriticalValueSurvival object

    ## in class_CriticalValueSurvival.R which create & initiaize an object of class "CriticalValueSurvival"
    ##  object is initialized with calculated values for survivalTime, sIndex, sFraction, and the input type

    return( .criticalValueSurvival(survivalTime = survivalTime,
                                   timePoints = timePoints,
                                   type = criticalValue) )
  } else {
    # if survivalTime is not given as input, create a CriticalValueMean object

    ## otherwise if survivalTime input is NULL, create CriticalValueMean object
    ## creates and returns a new object of class "CriticalValueMean"
    ## function in class_CriticalValueMean.R

    return( .criticalValueMean() )
  }

}
