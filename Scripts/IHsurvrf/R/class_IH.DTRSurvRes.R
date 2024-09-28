

## define a new S4 class called "DTRSurv"
setClass(Class = "DTRSurvRes",

         ## slots (contain attributes/properties) that objects of this class contains
         slots = c(
           ## stores function call that created the output (tracks how an object was generated)
           ## slot named call & has type call
           "call" = "call",

           ## the results of the finalized, overall pooled forest from each of the stages
           "Forest1" = "ANY",

           ## slot named value that stores ANY type
           "Forest2" = "ANY",

           ## slot named value that stores ANY type
           "Forest3" = "ANY",

           ## slot named value that stores ANY type
           "Forest4" = "ANY",

           ## slot named value that stores ANY type
           "Forest5" = "ANY",

           ## slot named "params" that stores another class called "Parameters" created in another file class_Parameters.R
           "params" = "Parameters",

           ## intput and output long_data where the actions get updated
           "long_data" = "ANY",

           ## input and output matrix of optimal survival probs from previous iterations
           ## this is only used for convergence

           "prev_probs"= "matrix",

           ## this stores the number of iterations until the forest reaches convergence
           "n_it" = "ANY",

           ## this stores the matrix of the average differences between the survival curves for patients between iterations
           "avgKM_diff" = "matrix",

           ## we will need to call different slots in the matrix to be used as cutoffs

           "cutoff" = "numeric"
         ))
