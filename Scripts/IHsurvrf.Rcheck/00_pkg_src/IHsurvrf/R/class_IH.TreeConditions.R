

# Class to store parameters that regulate tree and specify analysis preferences
#
# Class is not exported and is for internal convenience only
#
#  @slot nTree An integer object. The number of trees to be generated in forest.
#
#  @slot nodeSize An integer object. The minimum number of cases allowed in a
#    node
#
#  @slot minEvent An integer object. The minimum number of events allowed in a
#    node
#
#  @slot pooled A logical object. TRUE indicates that treatment groups are to
#    be considered in combination.
#
#  @slot stratifiedSplit A number object. The coefficient phi for stratified
#    random spliting
#
# Getters
#   .NTree(object, ...) {new; defined}
#   .NodeSize(object, ...) {new; defined}
#   .MinEvent(object, ...) {new; defined}
#   .Pooled(object, ...) {new; defined}
#   .StratifiedSplit(object, ...) {new; defined}
#
# Methods
#   .TreeConditionsAsList(object, ...) {new; defined}
#
# Functions
# .treeConditions(..., nTree, nodeSize, minEvent,
#                 pooled, stratifiedSplit)
#

## define a new S4 class called "TreeConditions"

setClass(Class = "TreeConditions",

         ## An integer object. The number of trees to be generated in forest.

         slots = c(nTree = "integer",

                   ## nodeSize An integer object. The minimum number of cases allowed in a node
                   nodeSize = "integer",

                   ## An integer object. The minimum number of events allowed in a node
                   minEvent = "integer",

                   ## TRUE indicates that treatment groups are to be considered in combination.
                   pooled = "logical",

                   ## A number object. The coefficient phi for stratified random spliting
                   stratifiedSplit = "numeric"))

## Getters

#-------------------------------------------------------------------------------
# Method to retrieve number of trees in the forest
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------

## create new generic function called .NTree
setGeneric(name = ".NTree",
           def = function(object, ...) {
             standardGeneric(".NTree")
           })

## establish a new method for .NTree() that prevents function from being used to object without specific class

setMethod(f = ".NTree",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## establish a new method for .Tau working on objects of class TreeConditions

setMethod(f = ".NTree",
          signature = c(object = "TreeConditions"),

          ## return the nTree slot (integer) : number of trees to be generated in a forest

          definition = function(object, ...) { return( object@nTree ) })

## write a function called .treeConditions which outputs a "TreeConditions" object
## this basically acts as a parameter check that all the inputs are within the valid range and in the proper format

.treeConditions <- function(...,

                            ## each of these params corresponds to slot in "TreeConditions" class
                            nTree,
                            nodeSize,
                            minEvent,
                            pooled,
                            stratifiedSplit) {

  ## checks if minimum number of events is numeric,

  # minimum number of events must be integer and > 0
  ## if not, stop running the function and output an error message
  if (!is.numeric(x = minEvent)) stop("minEvent must be integer", call. = FALSE)

  ## convert minimum number of events to an integer
  minEvent <- as.integer(x = minEvent)

  ## checks if minimum number of events is < 1 (it must be), otherwise stop running and output an error message
  if (minEvent < 1L) stop("minEvent must be non-zero positive", call. = FALSE)

  ## checks if the node size is numeric

  # minimum number of cases in each node, must be integer and > 0
  if (!is.numeric(x = nodeSize)) stop("nodeSize must be integer", call. = FALSE)

  ## converts it to an integer
  nodeSize <- as.integer(x = nodeSize)

  ## checks if minimum number of events is < 1 (it must be), otherwise stop running and output an error message
  if (nodeSize < 1L) stop("nodeSize must be non-zero positive", call. = FALSE)

  ## perform the same check on the number of trees

  # number of trees to grow in forest, must be integer and > 0
  if (!is.numeric(x = nTree)) stop("nTree must be integer", call. = FALSE)
  nTree <- as.integer(x = nTree)
  if (nTree < 1L) stop("nTree must be non-zero positive", call. = FALSE)

  ## checks if the "pooled" input is a logical and not NA

  if (!is.logical(x = pooled) || is.na(x = pooled)) {

    ## otherwise, output an error message
    stop("pooled must be logical", call. = FALSE)
  }

  ## check if stratifiedSplit is null or too small

  if (is.null(x = stratifiedSplit) || stratifiedSplit <= 1e-8) {

    ## if this is the case, adjust the value to 0
    stratifiedSplit <- 0.0
  } else {

    ## otherwise, if the value is greater than 1
    if (stratifiedSplit > 1.0) {

      ## outut an error message

      stop("stratifiedSplit must be [0,1]", call. = FALSE)
    }
  }

  ## return an object of class TreeConditions

  return( new(Class = "TreeConditions",
              "nTree" = nTree,
              "nodeSize" = nodeSize,
              "minEvent" = minEvent,
              "pooled" = pooled,
              "stratifiedSplit" = stratifiedSplit) )

}

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of cases allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------

## used in IH.dtrSurv.R in setUpBasics

## create new generic function called .NodeSize
setGeneric(name = ".NodeSize",
           def = function(object, ...) { standardGeneric(".NodeSize") })


## establish a new method for .NodeSize() that prevents function from being used to object without specific class

setMethod(f = ".NodeSize",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## establish a new method for .NodeSize working on objects of class TreeConditions

setMethod(f = ".NodeSize",
          signature = c(object = "TreeConditions"),

          ## return the nodeSize slot (integer) : The minimum number of cases allowed in a node
          definition = function(object, ...) { return( object@nodeSize ) })


#-------------------------------------------------------------------------------
# Method to retrieve minimum number of events allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------

## used in IH.dtrSurv.R in setUpBasics
## create new generic function called .MinEvent

setGeneric(name = ".MinEvent",
           def = function(object, ...) { standardGeneric(".MinEvent") })


## establish a new method for .MinEvent() that prevents function from being used to object without specific class

setMethod(f = ".MinEvent",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## establish a new method for .MinEvent working on objects of class "TreeConditions"

setMethod(f = ".MinEvent",
          signature = c(object = "TreeConditions"),

          ## return the minEvent slot (integer): minimum number of events allowed in a node

          definition = function(object, ...) { return( object@minEvent ) })


#-------------------------------------------------------------------------------
# Method to retrieve flag for pooled analysis
#-------------------------------------------------------------------------------
# Method returns a logical
#-------------------------------------------------------------------------------

## create new generic function called .Pooled()

setGeneric(name = ".Pooled",
           def = function(object, ...) { standardGeneric(".Pooled") })

## defined new method for .Pooled that prevents function from being used  object without specific class

setMethod(f = ".Pooled",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

## establish a new method for .Pooled working on objects of class "TreeConditions"

setMethod(f = ".Pooled",
          signature = c(object = "TreeConditions"),

          ## returns the pooled slot (logical): TRUE indicates that treatment groups are to be considered in combination.
          definition = function(object, ...) { return( object@pooled ) })
