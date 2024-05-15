

# Verify input 'uniformSplit'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'uniformSplit' is provided logical or NULL. If NULL set to value
#  of ERT 
#
# successful methods return a logical object.
#

## define a new generic function called .VerifyUniformSplit
## used in class_TreeType.R


setGeneric(name = ".VerifyUniformSplit",
           def = function(uniformSplit, ...) { 
             standardGeneric(".VerifyUniformSplit") 
           })

## throw an error if used on a non-specific class

# the default method generates an error
setMethod(f = ".VerifyUniformSplit",
          signature = c(uniformSplit = "ANY"),
          definition = function(uniformSplit, ...) { 
            stop("uniformSplit must be logical or NULL", call. = FALSE)
          })

## create method applying to objects of "logical" class

setMethod(f = ".VerifyUniformSplit",
          signature = c(uniformSplit = "logical"),
          definition = function(uniformSplit, ...) { 
            
            ## if the input (uniformSplit) is NA, 
            if (is.na(x = uniformSplit)) {
              
              ## stop running and return an error message
              
              stop("uniformSplit must be logical or NULL", call. = FALSE)
            }
            
            ## otherwise if it's not NA, return the output
            return( uniformSplit ) 
          })

## create method applying to class of NULL objects, also takes into account ERT input

setMethod(f = ".VerifyUniformSplit",
          signature = c(uniformSplit = "NULL"),
          definition = function(uniformSplit, ..., ERT) { 
            
            ## if ERT is null, as well as uniformSplit, 
            if (is.null(x = ERT)) {
              
              ## stop running and output an error
              
              stop("if uniformSplit = NULL, ERT must be set", call. = FALSE)
            }
            
            ## otherwise, if ERT isn't null, call uniformSplit, while setting the value of uniformSplit to ERT (a logical)
            return( .VerifyUniformSplit(uniformSplit = ERT, ...) ) 
          })
