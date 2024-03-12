
# Verify input 'txName'
#
# method is not exported and is for internal convenience only
#
# ensures that 'txName' is provided as a character or character vector and
# that the provided names are present in data. This input defines the
# number of decision points for the analysis.
#
# successful methods return the original input without modification.
#


## define a new generic function called .VerifyTxName
## used in dtrSurv.R Script

setGeneric(name = ".VerifyTxName",
           def = function(txName, ...) { standardGeneric(".VerifyTxName") })

## create method so that when called on objects of class "ANY", an error is returned

# the default method generates an error
setMethod(f = ".VerifyTxName",
          signature = c(txName = "ANY"),
          definition = function(txName, ...) {
            stop("txName must be a vector of character objects",
                 call. = FALSE)
          })

## create method to operate on objects of class "character"

setMethod(f = ".VerifyTxName",
          signature = c(txName = "character"),
          definition = function(txName, ..., data) {

            ## if there is no txName provided, stop execution with an error message

            if (length(x = txName) == 0L) {
              stop("txName must be provided", call. = FALSE)
            }



            # if treatment name is in data, it is *the* treatment name
            # if it is not, it is a single item indicating the treatment
            # variable name in the common formula
            # We assume a dot between name and decision point

            test <- tryCatch(
              ## try subsetting the data using "txName"
              expr = data[,txName,drop = FALSE],

              ## if an error occurs during the execution process (txName isn't found), return a NULL
              error = function(e) { return( NULL ) })

            ## if the subset is null, and the treatment name is a single item,
            ## assumes txName indicates treatment variable name and processes accordingly

            if (is.null(x = test) && length(x = txName) == 1L) {

              ## retrieve the column names of the "data" object to store in "dataNames"

              dataNames <- colnames(x = data)

              # split column ames based on dots which results in  list of character vectors

              # split the data.frame names on dots which indicate stages
              cov <- strsplit(x = dataNames, split = ".", fixed = TRUE)

              ## applies function to each element of "cov"
              ## compares the first element of each character vector in "cov" with provided "txName"
              ## result is a logical list

              # assume that first element is the common name
              areAs <- lapply(X = cov, FUN = function(x){x[[ 1L ]] == txName})

              ## checks if at least one element in this is true (if first element is common name)

              if (sum(areAs) > 0L) {

                ## calculates number of decision points by summing the true values

                nDP <- sum(areAs)

                ## return message showing number of decision points

                message("detected ", nDP, "decision points")

                ## updates txName to be common names found

                txName <- dataNames[areAs]

                ## calls function again for verification

                return( .VerifyTxName(txName = txName, data = data) )
              }
            }

            ## try subsetting the data again. If this is unsuccessful, stop execution with an error message

            test <- tryCatch(expr = data[,txName,drop = FALSE],
                             error = function(e) {
                               stop("unable to retrieve 'txName' from data",
                                    e$message, call. = FALSE)
                             })

            ## check if txName includes any NaN values
            ## if it does, stop execution with error message

            if (any(sapply(X = test, FUN = is.nan))) {
              stop("txName cannot include NaN values", call. = FALSE)
            }

            ## check if treatment variable is a factor of integer like
            ## if not, stop execution with an error message

            # ensure tx is factor or integer-like

            ## a loop that iterates over each column of the "test" data frame
            ## generates a sequence from 1:number of columns in "test"
            for (i in 1L:ncol(x = test)) {

              ## if column "i" isn't a factor,
              if (!is.factor(x = test[,i])) {

                ## and if column "i' is a numeric,
                if (is.numeric(x = test[,i])) {

                  ## and if the original numeric values in column "i" are approximately equal oto the rounded values
                  if (!isTRUE(all.equal(target = test[,i],
                                        current = round(x = test[,i], digits = 0L)))) {

                    ## STOP if the original and rounded values are not approximately equal
                    ## this means that the column contains non-integer numeric values

                    stop("treatment variable must be integer or factor",
                         call. = FALSE)
                  }

                  ## return this error if the column "i" isn't a factor and is NOT a numeric
                } else {
                  stop("treatment variable must be integer or factor",
                       call. = FALSE)
                }
              }
            }

            ## if all checks pass, return the processed "txName"

            return( txName )
          })
