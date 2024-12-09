

# Verify input 'criticalValue'
#
# methods are not exported and are only for internal convenience
#
# ensures that 'criticalValue' is one of {'mean', 'surv.prob', 'surv.man'}.
#
# successful methods return the original character possibly modified to be
# all lower case.


## function used in criticalValue.R


## create a new generic function called ".VerifyCriticalValue"

setGeneric(name = ".VerifyCriticalValue",
           def = function(criticalValue, ...) {
             standardGeneric(".VerifyCriticalValue")
           })

## method throws an error if it gets executed for a class that isn't appropriate

# the default method generates an error
setMethod(f = ".VerifyCriticalValue",
          signature = c(criticalValue = "ANY"),
          definition = function(criticalValue, ...) {
            stop("criticalValue must be one of ",
                 "{'mean'}",
                 call. = FALSE)
          })

## define a new method used or class "character"

setMethod(f = ".VerifyCriticalValue",
          signature = c(criticalValue = "character"),
          definition = function(criticalValue, ...) {

            ## converts the input "criticalValue" to lowercase to standardize input

            criticalValue <- tolower(x = criticalValue)

            ## checks if the input is one of these

            if (criticalValue %in% c("mean"))

              ## if so, returns the lowercase ctitical value
              return( criticalValue )

            ## otherwise, if it doesn't match one of the accepted values, returns an error

            stop("criticalValue must be one of ",
                 "{'mean'}",
                 call. = FALSE)

          })
