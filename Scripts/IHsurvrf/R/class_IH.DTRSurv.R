


## Roxygen2 documentation tag. Documentation for this class (DTRSurv) will include contents from class_DTRSurvStep.R file
#' @include class_IH.DTRSurvStep.R


## define a new S4 class called "DTRSurv"
setClass(Class = "DTRSurv",

         ## slots (contain attributes/properties) that objects of this class contains
         slots = c(
           ## stores function call that created the output (tracks how an object was generated)
           ## slot named call & has type call
           "call" = "call",

           ## the results of the finalized, overall pooled forest from each of the stages
           "FinalForest" = "ANY",

           ## slot named "params" that stores another class called "Parameters" created in another file class_Parameters.R
           "params" = "Parameters",

           ## this stores the number of iterations until the forest reaches convergence
           "n_it" = "ANY",

           ## intput and output long_data where the actions get updated
           "long_data" = "ANY",

           ## input and output matrix of optimal survival probs from previous iterations
           ## this is only used for convergence

           "prev_probs"= "matrix"))



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
#' dt <- data.frame("Y_1" = sample(1:100,100,TRUE), "Y_2" = sample(101:200,100,TRUE),
#'                  "D_1" = rbinom(100, 1, 0.9), "D_2" = rbinom(100,1,0.9),
#'                  "A_1" = rbinom(100, 1, 0.5), "A_2" = rbinom(100,1,0.5),
#'                  "X_1" = rnorm(100), "X_2" = rnorm(100), subj.id = 1:100)
#'
#' result <- IHsurvrf(data = dt,
#'                   txName = c("A_1", "A_2"),
#'                   models = Surv(Y,D)~X+A, stageLabel = "_",)
#'
#' tt <- predict(object = result@Forest1)
#' tt <- predict(object = result@Forest1, findOptimal = FALSE)

### .Predict() takenfrom class_IH.DTRSurvStep.R which then uses the function from class_IH.SurvRF.R

setMethod(f = "predict",
          signature = c(object = "DTRSurv"),
          definition = function(object,
                                ...,
                                newdata,
                                findOptimal = TRUE) {


            ## if there's no new data (it's null), jusr retrieve the fitted values
            ## method returns the contents of the "forest" slot from the SurvRF object from the final forest
            ## this object contains list object containing survFunc, mean

            if (missing(x = newdata)) {
              ## the input is of class DTRSurvStep for the individual stage
              ### this function is defined in class_DTRSurvStep.R

              ### use the contents from the FINAL FOREST

              return( .Predict(object = object@FinalForest,
                               newdata = NULL,
                               params = object@params,
                               findOptimal = findOptimal) )
            } else {

              ## otherwise, if findOptimal = TRUE, return the optimal predictions using .PredictAll():
              ## for all possible treatment options-- go through each treatment ant see what the patient's estimated survival would have been
              ### NOTE: we use the Final Forest
              ## ## return a list containing the predictions for each treatment level and the optimal treatment decision

              ## if findOptimal = FALSE, use .Predict() calculate the estimatated survival values that the patient receives using the actions they already got



              return( .Predict(object = object@FinalForest,
                               newdata = newdata,
                               params = object@params,
                               findOptimal = findOptimal) )
            }

          })
