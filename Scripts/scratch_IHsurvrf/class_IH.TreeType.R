
# Class to store information regarding tree type ERT or Breiman
#
# Class is not exported and is for internal convenience only
#
#  @slot replace A logical object; TRUE indicates that samples can include
#    duplicate records
#
#  @slot randomSplit A numeric object; The probability of a random split
#
#  @slot ERT A logical object; TRUE indicates that extremely randomized trees
#    methods are to be used
#
#  @slot uniformSplit A logical object; TRUE indicates that cutoffs are to 
#    be selected based on a uniformed random number
#
#  @slot splitRule A character object; must be one of {'logrank', 'mean'}
#
#  @slot tieMethod A character object; must be one of {'first', 'random', 'NA'}
#
# Getters
#
# Methods
#
# Functions
#
#' @include IH.VerifyERT.R IH.VerifyRandomSplit.R 
#' @include IH.VerifyUniformSplit.R IH.VerifyReplace.R IH.VerifySplitRule.R
#' @include IH.VerifyTieMethod.R

source("IH.VerifyERT.R")
source("IH.VerifyRandomSplit.R")
source("IH.VerifyUniformSplit.R")
source("IH.VerifyReplace.R")
source("IH.VerifySplitRule.R")
source("IH.VerifyTieMethod.R")



## create a new S4 class called "TreeType"

setClass(Class = "TreeType",
         
         ## A logical object; TRUE indicates that samples can include duplicate records
         slots = c("replace" = "logical",
                   
                   ## A numeric object; The probability of a random split 
                   
                   "randomSplit" = "numeric",
                   
                   ## A logical object; TRUE indicates that extremely randomized trees methods are to be used
                   "ERT" = "logical",
                   
                   ## A logical object; TRUE indicates that cutoffs are to be selected based on a uniformed random number
                   "uniformSplit" = "logical",
                   
                   ## A character object; must be one of {'logrank', 'mean'}
                   "splitRule" = "character",
                   
                   ## A character object; must be one of {'first', 'random', 'NA'}
                   "tieMethod" = "character"))

## Getters

## create a new function called .treeType that returns a object of class "TreeType" with validated params

# initializer
.treeType <- function(
    ## Logical object. TRUE = use ERT to set the candidate variable
    ERT, 
                      
    ## not a slot in above function
    ## gotten from code: number of individuals in the dataset 
    nSamples,  
    
    ## Logical object. if TRUE, and ERT = TRUE, sample random cutoff from uniform distribution
    ## otherwise if FALSE, and ERT = TRUE, select random case & select cutoff point from mean cutoff between it and next largest covar
    ## ignore this input if ERT = FALSE
    uniformSplit,  
    
    ## logical object. TRUE = sample with replacement (indiv can be present more than once in the sample)
    ## if null, 'replace' = !'ERT'
    replace,  
    
    # character object of "logrank" or "mean"
    ## if value is NULL and criticalValue == "mean"; this will use "mean"
    ## if value is NULL and criticalValue == "surv.prob" or "surv.mean"; this will use logrank
    splitRule,  
    
    ## character object of "first" or "random"
    tieMethod,
    
    # [0,1], probability that a random split will occur
    randomSplit,
    
    ## not a slot in above function
    ## character object of "mean", "surv.prob" or "surv.mean"
    criticalValue) {
  
  # ensure that ERT is logical or NULL. Methods return a logical.
  
  ## .VerifyERT function used in VerifyERT.R script 
  
  ERT <- .VerifyERT(ERT = ERT)
  
  ## .VerifyRandomSplit.R script
  
  # ensure that randomSplit is 0 <= rs < 1. Methods return a numeric.
  randomSplit <- .VerifyRandomSplit(randomSplit = randomSplit)
  
  
  ## .VerifyUniformSplit.R script
  
  # ensure that uniformSplit is logical or NULL. Methods return a logical.
  uniformSplit <- .VerifyUniformSplit(uniformSplit = uniformSplit, ERT = ERT)
  
  ## .VerifyReplace.R Script
  
  # ensure that replace is logical or NULL. Methods return a logical.
  replace <- .VerifyReplace(replace = replace, ERT = ERT)
  
  ## ensures the input split rule is either "logrank" or "mean", convert to lowercase
  ## VerifySplitRule.R script
  
  # verify splitRule. methods return the original character object with possible
  # modification to make all lower case
  ## if split rule is NULL, go off the input for critical value
  
  splitRule <- .VerifySplitRule(splitRule = splitRule, 
                                criticalValue = criticalValue)
  
  ## check tieMethod is one of the allowed values "first", "random", "NA
  
  # successful methods return the original character input possibly modified to
  # lower case
  
  ## used in VerifyTieMethod.R Script (must be "random" or "first")
  
  tieMethod <- .VerifyTieMethod(tieMethod = tieMethod)
  
  ## create an object of class "TreeType" with the verified input parameters
  
  return( new(Class = "TreeType", 
              "replace" = replace, 
              "randomSplit" = randomSplit,
              "ERT" = ERT,
              "uniformSplit" = uniformSplit,
              "splitRule" = splitRule,
              "tieMethod" = tieMethod) )
  
}
