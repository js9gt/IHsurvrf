

# Class extends CriticalValue to indicate that critical value is non-survival
#  mean

#### AKA mean survival time 

#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueCriterion(object, ...) {defined}
#  .CreateValueObject(object, ...) {defined}
#
# Functions
#  .criticalValueMean(...)
#

## this part includes any documentation from the class_CriticalValue.R function, keeps related documentation together
## for the sake of package development, using the Roxygen2 comment generates the R documentation files
#' @include class_IH.CriticalValue.R

source("class_IH.CriticalValue.R")


## defines a new class called "CriticalValueMean" which is a subclass of CriticalvalueBase
## meaning, the new class, CriticalValueMean, inherits properties from the CriticalValueBase class (class_CriticalValue.R)
setClass(Class = "CriticalValueMean",
         contains = c("CriticalValueBase"))


#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueMean
#-------------------------------------------------------------------------------
# function returns a CriticalValueMean object
#-------------------------------------------------------------------------------

## function creates and returns a new object of class "CriticalValueMean"
## uses object-oriented programming 
.criticalValueMean <- function(...) { return( new("CriticalValueMean") ) }