

# Verify input 'replace'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'replace' is logical or NULL. If NULL, replace = ERT.
#
# successful methods return a logical indicating if replacement should be
# used in sampling
#


## used in class_IH.TreeType.R

## create a new generic function called .VerifyReplace
setGeneric(name = ".VerifyReplace",
           def = function(replace, ...) { standardGeneric(".VerifyReplace") })


## generate an error if not called on more specific object classes

# the default method generates an error
setMethod(f = ".VerifyReplace",
          signature = c(replace = "ANY"),
          definition = function(replace, ...) {
            stop("replace must be logical or NULL", call. = FALSE)
          })

## new method for .VerifyReplace to act on objects of class "logical"

setMethod(f = ".VerifyReplace",
          signature = c(replace = "logical"),
          definition = function(replace, ...) {

            ## if the inpput value of replace is NA,

            if (is.na(x = replace)) {

              ## output an error

              stop("replace must be logical or NULL", call. = FALSE)
            }

            ## otherwise if input isn't NA, return the original value
            return( replace )
          })

## new method for .VerifyReplace to act on objects of class "NULL"


setMethod(f = ".VerifyReplace",
          signature = c(replace = "NULL"),

          ## function also considers the ERT argument
          definition = function(replace, ..., ERT) {

            ## if replace = NULL and ERT is null as well,
            if (is.null(x = ERT)) {

              ## output an error: there neds to be a value of ERT to decide the behavior of "replace"
              stop("if replace = NULL, ERT must be set", call. = FALSE)
            }

            ## otherwise, if replace = NULL and ERT is not null,
            ## run the .VerifyReplace function using the opposite of ERT as the "replace" value
            ## if ERT are used, the default would be to not use replacement sampling

            return( .VerifyReplace(replace = !ERT, ...) )
          })
