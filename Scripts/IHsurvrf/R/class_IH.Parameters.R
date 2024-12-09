
# Class to store parameters that regulate tree and specify analysis preferences
#
# Class is not exported and is for internal convenience only
#
# Methods
#   .ParametersAsList(object, ...) {new; defined}
#
# Functions
# .parameters(timePoints, nTimes, response, nTree, ERT, uniformSplit,
#                      randomSplit, splitRule, replace, nodeSize,
#                      minEvent, tieMethod, criticalValue,
#                       nSamples, stratifiedSplit)



source("R/class_IH.TimeInfo.R")
#source("~/survrf/Scripts/IHsurvrf/R/class_IH.TimeInfo.R")
source("R/class_IH.CriticalValue.R")
#source("~/survrf/Scripts/IHsurvrf/R/class_IH.CriticalValue.R")
source("R/IH.criticalValue.R")
#source("~/survrf/Scripts/IH.criticalValue.R")
source("R/class_IH.TreeType.R")
#source("~/survrf/Scripts/IHsurvrf/R/class_IH.TreeType.R")
source("R/class_IH.TreeConditions.R")
#source("~/survrf/Scripts/IHsurvrf/R/"class_IH.TreeConditions.R)

#
#-------------------------------------------------------------------------------
# Class creation for "Parameters_mean" and "Parameters_Survival" (for the output)
# then create a class union called "Parameters"
#-------------------------------------------------------------------------------
# contain slots for time information, mean critical values, tree type info, tree splitting
#-------------------------------------------------------------------------------


## create a new S4 class called "Parameters_mean"
## note: creating a new class is different than creating a generic function

setClass(Class = "Parameters_Mean",

         ## contains slots for time information, mean critical values, tree type information, tree splitting conditions

         contains = c("TimeInfo", "CriticalValueMean", "TreeType", "TreeConditions"))



## create a class union called "Parameters" that can refer to either "Parameters_Mean" or "Parameters_Survival"
## functions can accept either class type with this

setClassUnion(name = "Parameters",
              members = c("Parameters_Mean"))

setMethod(f = "initialize",
          signature = c(.Object = "Parameters_Mean"),
          def = function(.Object, ...) {

            obj <- list(...)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("unrecognized object sent to Parameters object")
              }
            }
            return( .Object )

          })



#-------------------------------------------------------------------------------
# Function to verify inputs and create a Parameters object
#-------------------------------------------------------------------------------
# Function returns a Parameters object
#-------------------------------------------------------------------------------

## defines a function to verify inputs, and create a "Parameters" object (class union)
.parameters <- function(
    # A character object or a numeric vector object. If a character
    #'   object, must be one of \{"quad", "uni", "exp"\} indicating the distribution
    #'   from which the time points are to be calculated. For character input,
    #'   input 'nTimes' must also be provided. If a numeric vector, the
    #'   time points to be used. If 0 is not the first value, it will be
    #'   concatenated by the software.
    #'   ############## not sure if we need this
    timePoints,

    # study legnth
    tau,

    # integer: total number of timepoints to be used
    ## this is picked to be the number between 0 and tau
    nTimes,

    ## extracts the "response variable": matrix of survival response variables from the model that will be input
    response,

    ## An integer object. The number of trees to grow.
    nTree,

    ## Logical object. TRUE = use ERT to set the candidate variable
    ERT,

    ## Logical object. if TRUE, and ERT = TRUE, sample random cutoff from uniform distribution
    ## otherwise if FALSE, and ERT = TRUE, select random case & select cutoff point from mean cutoff between it and next largest covar
    ## ignore this input if ERT = FALSE
    ## if this is NULL, ERT must be set
    ## if this is NULL and there is ERT value, use the same as ERT input
    uniformSplit,

    # [0,1], probability that a random split will occur
    randomSplit,

    # character object of "logrank" or "mean"
    ## if value is NULL and criticalValue == "mean"; this will use "mean"
    ## if value is NULL and criticalValue == "surv.prob" or "surv.mean"; this will use logrank

    splitRule,

    ## logical object. TRUE = sample with replacement (indiv can be present more than once in the sample)
    ## if null, 'replace' = !'ERT'
    replace,

    ## integer object. Minimum number of individuals present in a node
    nodeSize,

    ## integer object. Minimum number of events present in a node
    minEvent,

    ## character object of "first" or "random"
    tieMethod,

    ## character object of "mean"
    criticalValue,

    ## gotten from code: number of individuals in the dataset
    nSamples,

    ## numeric. Stratified random split coefficient.
    stratifiedSplit) {

  ## calls a function ".timeInfo" in class_TimeInfo.R script that returns info about specified time points, tau, nTime


  # initialize TimeInfo
  # function returns a TimeInfo object
  timeInfo <- .timeInfo(timePoints = timePoints,
                        tau = tau,
                        nTimes = nTimes,
                        response = response)


  ## converts criticalValue argument to lowercase

  cv <- tolower(criticalValue)

  ## .criticalValue function in criticalValue.R
  ## initializes a critical value object of either "CriticalValueMean" or "CriticalValueSurvival" based on input "survivalTime" parameter

  # initialize CriticalValue which verifies inputs and
  # function returns an object of class CriticalValueMean or
  #   CriticalValueSurvival depending on input survivalTime
  ## CriticalValueSurvival would output "sIndex" and "sFraction" to be used in tSurvTree to calculate survival probability
  ## these values are used in setUpBasics
  ## class CriticalValueMean or CriticalValueSurvival

  criticalValue <- .criticalValue(criticalValue = criticalValue)



  ## .treeType function in class_TreeType.R script
  ## stores information about tree structure, splitting critiera,

  # initialize tree type info
  # function returns a TreeType object
  treeType <- .treeType(ERT = ERT,
                        nSamples = nSamples,
                        uniformSplit = uniformSplit,
                        replace = replace,
                        randomSplit = randomSplit,
                        splitRule = splitRule,
                        tieMethod = tieMethod,
                        criticalValue = cv)



  ## .treeConditions function in class_TreeConditions.R script
  # Method returns a list containing 6 elements
  #   "nTree" an integer, the total number of trees in forest
  #   "nodeSize" an integer, the minimum number of cases in a node
  #   "minEvent" an integer, the minimum number of events in a node
  #   "stratifiedSplit" the coefficient phi for stratified


  # initialize tree conditions info
  # function returns a TreeConditions object
  treeConditions <- .treeConditions(

    ## these are all input parameters
    nTree = nTree,
    nodeSize = nodeSize,
    minEvent = minEvent,
    stratifiedSplit = stratifiedSplit)


    ## otherwise, return "Parameters_Mean" class
    ## CriticalValueMean class initiated in class_CriticalValueMean.R Script

    return( new(Class = "Parameters_Mean",
                timeInfo,
                criticalValue,
                treeType,
                treeConditions) )
  }


