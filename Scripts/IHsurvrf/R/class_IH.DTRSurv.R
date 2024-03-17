

source("R/class_IH.DTRSurvStep.R")


## define a new S4 class called "DTRSurv"
setClass(Class = "DTRSurv",

         ## slots (contain attributes/properties) that objects of this class contains
         slots = c(
           ## stores function call that created the output (tracks how an object was generated)
           ## slot named call & has type call
           "call" = "call",

           ## slot named stageResults that stores a list
           "stageResults" = "list",

           ## slot named value that stores ANY type
           "value" = "ANY",

           ## slot named "params" that stores another class called "Parameters" created in another file class_Parameters.R
           "params" = "Parameters"))



#' Prediction Method: used in data generation process after estimation, when we want evaluation
#'
#' Method to estimate the value for new data or to retrieve estimated value for
#'  training data
#'
#' param object A DTRSurv object. The object returned by a call to dtrSurv().
#'
#' param ... Ignored. Used to require named inputs.
#'
#' param newdata NULL or a data.frame object. If NULL, this method retrieves
#'   the estimated value for the training data. If a data.frame, the
#'   value is estimated based on the data provided.
#'
#' param stage An integer object. The stage for which predictions are desired.
#'
#' param findOptimal A logical object. If TRUE, the value is estimated for
#'   all treatment options and that leading to the maximum value for each
#'   individual is used to estimate the value.
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' result <- dtrSurv(data = dt,
#'                   txName = c("A.1", "A.2"),
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1,
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' tt <- predict(object = result)
#' tt <- predict(object = result, stage = 1)
#' tt <- predict(object = result, findOptimal = FALSE)
#' tt <- predict(object = result, newdata = dt)
#' tt <- predict(object = result, newdata = dt, stage = 1)
#' tt <- predict(object = result, newdaata = dt, findOptimal = FALSE)

### .Predict() takenfrom class_IH.DTRSurvStep.R which then uses the function from class_IH.SurvRF.R

setMethod(f = "predict",
          signature = c(object = "DTRSurv"),
          definition = function(object,
                                ...,
                                newdata,
                                stage = 1,
                                findOptimal = TRUE) {

            if (stage > length(x = object@stageResults)) {
              stop("requested stage not present in analysis", call. = FALSE)
            }

            ## if there's no new data (it's null), jusr retrieve the fitted values
            ## method returns the contents of the "forest" slot from the SurvRF object
            ## this object contains list object containing survFunc, mean, and? survProb

            if (missing(x = newdata)) {
              ## the input is of class DTRSurvStep for the individual stage
              ### this function is defined in class_DTRSurvStep.R

              return( .Predict(object = object@stageResults[[ stage ]],
                               newdata = NULL,
                               params = object@params,
                               findOptimal = findOptimal) )
            } else {

              ## otherwise, if findOptimal = TRUE, return the optimal predictions using .PredictAll():
              ## for all possible treatment options-- go through each treatment ant see what the patient's estimated survival would have been
              ## ## return a list containing the predictions for each treatment level and the optimal treatment decision

              ## if findOptimal = FALSE, use .Predict() calculate the estimatated survival values that the patient receives using the actions they already got

              return( .Predict(object = object@stageResults[[ stage ]],
                               newdata = newdata,
                               params = object@params,
                               findOptimal = findOptimal) )
            }

          })
