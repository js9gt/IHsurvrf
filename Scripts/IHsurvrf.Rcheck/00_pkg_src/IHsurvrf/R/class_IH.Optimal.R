
# Class for storing estimated optimal treatment and value
#
# Class is not exported and is for internal convenience only
#
#  @slot optimalTx A vector object. The index of the estimated optimal tx
#
#  @slot optimalY A vector. The estimated value
#
#  @slot type A character. One of "mean" or "prob"
#
# Getters
#  .OptimalY(object, ...) {new; defined}
#  .OptimalAsList(object, ...) {new; defined}
#

## defines a new S4 class called "Optimal"
setClass(Class = "Optimal",

         ## 3 slots
         slots = c(

           ## vector that shows index of the estimated optimal treatment
           "optimalTx" = "vector",

           ## vector that shows estimated value
           "optimalY" = "matrix",

           ## nature of optimal values-- "mean" or "prob"
           "type" = "character"))

#-------------------------------------------------------------------------------
# Method to retrieve the estimated optimal value
#-------------------------------------------------------------------------------
# Method returns a vector object
#-------------------------------------------------------------------------------

## define a generic function ".OptimalY" used to retrieve estimated optimal values in "optimalY" slot
## generic function can operate on different classes
setGeneric(name = ".OptimalY",
           def = function(object, ...) { standardGeneric(".OptimalY") })


## define a method for function .OptimalY that throws error if called on objects of ANY class
## probably designed to have error if called on objects NOT of class "Optimal"

setMethod(f = ".OptimalY",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## method for .OptimalY specifically for "Optimal" class
## returns the optimalY slot of an "Optimal" object

setMethod(f = ".OptimalY",
          signature = c(object = "Optimal"),
          definition = function(object, ...) { return( object@optimalY ) })
