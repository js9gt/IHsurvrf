# Verify input 'survivalTime'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'survivalTime' is provided if criticalValue is one of 
# {'surv.prob', 'surv.mean'}. 
#
# successful methods return the numeric survivalTime or NULL.
#

## used in IH.criticalValue.R

## create a new generic function called .VerifySurvvialTime

setGeneric(name = ".VerifySurvivalTime",
           def = function(survivalTime, ...) { 
             standardGeneric(".VerifySurvivalTime") 
           })


## output an error if function isn't called on a more specific class

# the default method generates an error
setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "ANY"),
          definition = function(survivalTime, ...) { 
            stop("evalTime must be numeric or NULL", 
                 call. = FALSE)
          })

## create a new method of .VerifySurvivalTime for classes of type "numeric"

setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "numeric"),
          
          ## accounts for criticalValue and Tau as inputs
          
          definition = function(survivalTime, ..., criticalValue, tau) { 
            
            ## check if critical value is one of surv.prob or surv.mean
            
            if (!{criticalValue %in% c("surv.prob", "surv.mean")}) {
              
              ## if not, then criticalValue is "mean", so ignore and output error message
              
              message("evalTime is ignored if critical value is mean")
              return( NULL )
            }
            
            ## if more than one value for survival time is given,
            
            if (length(x = survivalTime) > 1L) {
              
              ## output an error
              
              stop("only 1 value for evalTime can be given",
                   call. = FALSE)
            }
            
            ## if survival time is < 0 or greater than tau,
            
            
            if (survivalTime <= 0.0 || survivalTime > tau) {
              
              ## output an error message and stop execution
              
              stop("evalTime must be between 0 and tau", call. = FALSE)
            }
            
            ## if survivalTime passes all checks, display message of its value
            
            message("evalTime ", survivalTime)
            
            ## also return the value
            
            return( survivalTime ) 
          })


## new method for .VerifySurvivalTime for cases where survivalTime is NULL

setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "NULL"),
          definition = function(survivalTime, ..., tau) { 
            
            ## set survival time to half of tau to provide a default value for survivalTime
            ## then call .VerifySurvivalTime to make sure this input is appropriate
            
            return( .VerifySurvivalTime(survivalTime = tau/2.0,
                                        tau = tau, ... ) )
          })
