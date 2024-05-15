

# Verify input 'tieMethod'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'tieMethod' is in {'random', 'first'}
#
# successful methods return the original input possibly modified to be lower
#   case
#

## used in class_TreeType.R


## create a new generic function called .VerifyTieMethod

setGeneric(name = ".VerifyTieMethod",
           def = function(tieMethod, ...) { 
             standardGeneric(".VerifyTieMethod") 
           })

## return an error if the function is called on an object without a more specific class

# the default method generates an error
setMethod(f = ".VerifyTieMethod",
          signature = c(tieMethod = "ANY"),
          definition = function(tieMethod, ...) { 
            stop("tieMethod must be one of {'random', 'first'}", 
                 call. = FALSE)
          })

## create a new method for objects of class "character"

setMethod(f = ".VerifyTieMethod",
          signature = c(tieMethod = "character"),
          definition = function(tieMethod, ...) { 
            
            ## turn the input tieMethod to lowercase
            
            tieMethod <- tolower(x = tieMethod)
            
            ## check to see if tieMethod is one of these
            
            if (!{tieMethod %in% c("random", "first")}) {
              
              ## if not, output an error
              
              stop("tieMethod must be one of {'random', 'first'}",
                   call. = FALSE)
            }
            
            ## if tieMethod meets the check, return it without modification
            
            return( tieMethod ) 
          })
