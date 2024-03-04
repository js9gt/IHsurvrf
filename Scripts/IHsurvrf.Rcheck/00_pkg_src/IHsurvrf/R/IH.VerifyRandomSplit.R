
# Verify input 'randomSplit'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'randomSplit' is numeric satisfying 0 < rs < 1. 
#
# successful methods return the object
#


## used in class_IH.TreeType.R script 


## create a new generic function called ".VerifyRandomSplit"

setGeneric(name = ".VerifyRandomSplit",
           def = function(randomSplit, ...) { 
             standardGeneric(".VerifyRandomSplit") 
           })

## output an error if called on object without specified class 
# the default method generates an error
setMethod(f = ".VerifyRandomSplit",
          signature = c(randomSplit = "ANY"),
          definition = function(randomSplit, ...) { 
            stop("randomSplit must obey 0 < randomSplit < 1", call. = FALSE)
          })

## define a new method to use .VerifyRandomSplit on objects of class "numeric"

setMethod(f = ".VerifyRandomSplit",
          signature = c(randomSplit = "numeric"),
          definition = function(randomSplit, ...) { 
            
            ## if the input consists of multiple values (vector), 
            
            if (length(x = randomSplit) > 1L) {
              
              ## then return an error
              
              stop("only 1 value for randomSplit can be given", call. = FALSE)
            }
            
            ## if the input is less than 0 or greater than 1, 
            
            if (randomSplit <= 0.0 || randomSplit >= 1.0) {
              
              ## then return an error
              
              stop("randomSplit must obey 0 < randomSplit < 1", call. = FALSE)
            }
            
            ## if randomSplit input passes these checks, return the original input
            
            return( randomSplit ) 
          })
