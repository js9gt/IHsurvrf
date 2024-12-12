

# Virtual class to store information regarding critical value selection
#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueAsList(object, ...) {not allowed}
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {not allowed}
#


## defining a virtual class called "CriticalValueBase"
## virtual class:
##                 class that isn't intended to be directly made, serves as a base class that other classes can inherit
##                 "child" classes can inherit properties from the parent class; useful for defining
##                         common structure or behavior that various derived classes can share


setClass(Class = "CriticalValueBase",
         contains = c("VIRTUAL"))


#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method is not defined for a general CriticalValue object
#-------------------------------------------------------------------------------

## creates generic method called ".CriticalValueCriterion"
## generic method:
##                 functions that behave differently based on the class they act upon
##  this allows method dispatch based on the class of the first argument to the function
##  allows the same function name to be used for different types of objects
##  the specifics for how the function behaves is determined by the class of object passed to it


setGeneric(name = ".CriticalValueCriterion",
           def = function(object, ...) { standardGeneric(".CriticalValueCriterion") })

## sets a method for the generic function ".CriticalValueCriterion"
## when the function is called with ANY object class, it should execute the function which outputs an error for general objects
## AKA irrespective of the class of object passed, the method defined for object = ANY will throw an error
## this catches errors for any object class where a more specific method is not defined


setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })


#-------------------------------------------------------------------------------
# CRITICAL VALUE MEAN
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


## defines a new class called "CriticalValueMean" which is a subclass of CriticalvalueBase
## meaning, the new class, CriticalValueMean, inherits properties from the CriticalValueBase class (class_CriticalValue.R)


setClass(Class = "CriticalValueMean",
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------


setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueMean"),
          definition = function(object, ...) { return( "mean" ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueMean
#-------------------------------------------------------------------------------
# function returns a CriticalValueMean object
#-------------------------------------------------------------------------------


## function creates and returns a new object of class "CriticalValueMean"
## uses object-oriented programming
.criticalValueMean <- function(...) { return( new("CriticalValueMean") ) }



