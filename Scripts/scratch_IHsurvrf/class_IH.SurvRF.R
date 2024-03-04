
# Virtual class to denote objects are arising from survRF step
#
# Methods
#   .Predict(object, newdata, ...) {new; not allowed}


## create a new, virtual class called "SurvRFObject"
## since it's virtual, you cannot create objects of this class, used as a base class from which other classes inherit
setClass(Class = "SurvRFObject",
         contains = c("VIRTUAL"))


## defines a new S4 class called "SurvRF" 

setClass(Class = "SurvRF",
         slots = c(
           
           ## A list object. The results of the tree building algorithm for
           ##    each tree in the forest
           "trees" = "list",
           
           ## A list object. The values averaged across all trees in the
           ##    forest
           "forest" = "list",
           
           ## A character vector. The variables considered in the
           ##    analysis
           "variables" = "character",
           
           ## An integer. The maximum number of covariates considered for 
           #    splitting
           "mTry" = "integer",
           
           ## An integer vector. The number of categories for each covariate
           #    considered. >=2 unordered factor, 1 ordered factor, 0 continuous
           "nCat" = "integer",
           
           ## A list object. The categories in each covariate considered.
           "xLevels" = "list"),
         
         ## SurvRF will inherit from SurvRFObject (virtual class)
         
         contains = c("SurvRFObject"))