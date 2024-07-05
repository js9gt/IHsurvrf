

## defines a new S4 class called pool1 for storing a single stage of IH random forest Q-learning survival analysis

setClass(
  Class = "IHsurvrf",

  ## slots (contain attributes/properties) that objects of this class contains

  slots = c(
    ## dataframe object which stores the results of the analysis in terms of actions and stages and actions at pooling steps
    "data" = "data.frame",

    ## stores the optimal estimated survival probabilities from the prior stage's analysis
    ## this includes the optimal ones as well as the appended ones from prior stage in a pattern
    "PMVoptimal" = "matrix",


    ## stores SurvRF object-- contains primary results of tree building algorithm

    ## pool1_results: the results of the pooled data for ALL prior data points (excluding current iteration of loop)
    "survRF1" = "ANY",

    ## append1_results: results including newly appended data (current iteration of loop)-- but this doesn't include fitting whole forest w newly appended data

    "survRF2" = "ANY"

  )
)
