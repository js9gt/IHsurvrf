
# Verify input 'splitRule'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'splitRule' is one of {'logrank', 'mean'}.
#
# successful methods return the original character object with possible
#  modification to make all lower case
#


## used in class_TreeType.R

## define a new generic function called .VerifySplitRule

setGeneric(name = ".VerifySplitRule",
           def = function(splitRule, ...) {
             standardGeneric(".VerifySplitRule")
           })

## method generates an error when being called on non-specific class

# the default method generates an error
setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "ANY"),
          definition = function(splitRule, ...) {
            stop("splitRule must be one of {'mean'}", call. = FALSE)
          })


## if splitRule (input) is null, set the splitRule based on the criticalValue critieria

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "NULL"),
          definition = function(splitRule, ..., criticalValue) {

            ## if the critical value is mean,

            if (criticalValue == "mean") {

              ## set the splitRule as mean
              splitRule = "mean"
            }

            ## return the input splitrule in lowercase from the character input below

            return( .VerifySplitRule(splitRule = splitRule, ...) )
          })

## new method for input object of type "character"

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "character"),
          definition = function(splitRule, ...) {

            ## changes input slitrule to lowercase

            splitRule <- tolower(x = splitRule)

            ## checks if splitRule is one of these

            if (!(splitRule %in% c("mean"))) {

              ## if not, stop and output an error
              stop("splitRule must be one of {'mean'}")
            }

            ## otherwise return an unmodified splitRule

            return( splitRule )
          })
