

# Function to verify inputs and create appropriate CriticalValue object
#
#  Function is not exported and for internal convenience only
#
# Function returns an object of class CriticalValueMean
#


source("R/class_IH.CriticalValue.R")
#source("~/survrf/Scripts/IHsurvrf/R/class_IH.CriticalValue.R")
source("R/IH.VerifyCriticalValue.R")
#source("~/survrf/Scripts/IHsurvrf/R/IH.VerifyCriticalValue.R")

## define a new function .criticalValue to verify input for critical values

.criticalValue <- function(criticalValue,
                           tau,
                           timePoints) {

  ## in VerifyCriticalValue.R script to make sure input criticalValue is of the correct class(?)

  # ensure criticalValue is one of {'mean', 'surv.prob', 'surv.mean'}.
  # Methods return the original character possibly modified to be lower case.

  criticalValue <- .VerifyCriticalValue(criticalValue = criticalValue)


  return( .criticalValueMean() )

}
