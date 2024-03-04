
# Verify input 'ERT'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'ERT' is provided as a logical or NULL object.
#
# successful methods return a logical indicating if the 
# Extremely Randomized Tree method is to be used
#

## used in script class_IH.TreeType.R


## create new generic method called .VerifyERT

setGeneric(name = ".VerifyERT",
           def = function(ERT, ...) { standardGeneric(".VerifyERT") })

# the default method generates an error
## error generation if not applied on a more specific data type
setMethod(f = ".VerifyERT",
          signature = c(ERT = "ANY"),
          definition = function(ERT, ...) { 
            stop("ERT must be logical or NULL", call. = FALSE)
          })

## define a new method for objects of class "logical"

setMethod(f = ".VerifyERT",
          signature = c(ERT = "logical"),
          definition = function(ERT, ...) { 
            
            ## if the value of ERT is NA,
            
            if (is.na(x = ERT)) {
              
              ## output an error
              
              stop("ERT must be logical or NULL", call. = FALSE)
            }
            
            ## otherwise, return the unmodified input 
            return( ERT ) 
          })

## define a new method for class "NULL" 

setMethod(f = ".VerifyERT",
          signature = c(ERT = "NULL"),
          
          ## when the input of ERT is "NULL", return a value of "TRUE"
          definition = function(ERT, ...) { return( TRUE ) })
