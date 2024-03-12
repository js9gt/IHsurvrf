

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
