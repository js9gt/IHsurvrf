

# Verify input 'sampleSize'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'sampleSize' is numeric, numeric vector, or NULL. If NULL, set
#   to a default value based on ERT selection.
#
# successful methods return a numeric vector object
#

## used in dtrSurv.R script

## create a new generic method called ".VerifySampleSize"
setGeneric(name = ".VerifySampleSize",
           def = function(sampleSize, ...) {
             standardGeneric(".VerifySampleSize")
           })

## output an error if called on something without a more specific class

# the default method generates an error
setMethod(f = ".VerifySampleSize",
          signature = c(sampleSize = "ANY"),
          definition = function(sampleSize, ...) {
            stop("sampleSize must be a numeric or NULL", call. = FALSE)
          })

## create new method for objects of class "numeric"

setMethod(f = ".VerifySampleSize",
          signature = c(sampleSize = "numeric"),

          ## also takes into the input nDP

          definition = function(sampleSize, ..., nDP) {

            ## if any of the sampleSize is less than 1e-8 or greater than 0,

            # sampleSize must be a fraction
            if (any({sampleSize < 1e-8 || sampleSize > 1.0})) {

              ## stop running and output an error

              stop("sampleSize must be 0 < sampleSize <= 1", call. = FALSE)
            }

            ## if the sample size is only a single value,

            # sampleSize must be provided for each decision point
            if (length(x = sampleSize) == 1L) {

              ## replicate this value for each decision point
              sampleSize <- rep(x = sampleSize, times = nDP)

              ## otherwise, if the input number of sampleSizes doesn't match the numbers of decision points,

            } else if (length(x = sampleSize) != nDP) {

              ## stop running and output an error

              stop("if sampleSize provided as vector, ",
                   "must provide values for all decision points",
                   call. = FALSE)
            }

            ## if sampleSize passes all these checks, it gets returned without modification

            return( sampleSize )
          })

## create new method for .VerifySampleSize for cases where the input value of sampleSize is NULL

setMethod(f = ".VerifySampleSize",
          signature = c(sampleSize = "NULL"),

          ## take as input, ERT and nDP

          definition = function(sampleSize, ..., ERT, nDP) {

            ## if the input for ERT (logical) is null,

            if (is.null(x = ERT)) {

              ## stop running code and output error
              ## default sampleSize might depend on whether ERT is used
              stop("if sampleSize = NULL, ERT must be set", call. = FALSE)
            }

            ## if ERT isn't null,

            ## if ERT is used, defaults sampleSize to 1.0 (using all data)
            ## uf ERT is not used, uses a default fraction of 0.632 of the data used

            sampleSize <- ifelse(test = ERT,
                                 yes = 1.0,
                                 no = 0.632)

            ## then, run the .VerifySampleSize with the input value of 1 or 0.632, and nDP

            return( .VerifySampleSize(sampleSize = sampleSize, ...,
                                      nDP = nDP) )

          })
