

# Verify input 'data'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'data' is provided as a data.frame or a matrix with named
# columns and that the object does not contain NaN values.
#
# successful methods return a data.frame object containing the data



## used in dtrSurv.R script

## create  new generic function called ".VerifyData"

setGeneric(name = ".VerifyData",
           def = function(data, ...) { standardGeneric(".VerifyData") })

## generate an error if the "data" argument doesn't match a more specific class
## AKA enforces requirement that "data" must be a specific type

# the default method generates an error
setMethod(f = ".VerifyData",
          signature = c(data = "ANY"),
          definition = function(data, ...) {
            stop("data must be a data.frame or a matrix with named columns",
                 call. = FALSE)
          })

## create a new method for .VerifyData for input of types data.frame

setMethod(f = ".VerifyData",
          signature = c(data = "data.frame"),
          definition = function(data, ...) {

            ## if there are any NaN values,

            if (any(sapply(X = data, FUN = is.nan))) {

              ## throw an error to ensure data doesn't have NaN values

              stop("data cannot include NaN values", call. = FALSE)
            }

            ## if data doesn't contain NaN values, it's returned without modification

            return( data )
          })

## create a new method for .VerifyData for input of type "matrix"

setMethod(f = ".VerifyData",
          signature = c(data = "matrix"),
          definition = function(data, ...) {

            ## if the matrix has unnamed columns,

            if (is.null(x = colnames(x = data))) {

              ## throw an error to enforce requirement for column names

              stop("if a matrix, data must include column headers",
                   call. = FALSE)
            }

            ## otherwise, if the matrix has column names, it's converted to a data.frame
            ## also, the previous .VerifyData function is called on it to check for NaN values
            ## returns an unmodified data frame

            return( .VerifyData(data = as.data.frame(x = data), ...) )
          })
