

# Verify inputs 'timePoints', 'tau', and 'nTimes'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'timePoints' is numeric or character, generates time points if
# appropriate, ensures that 'nTimes' is appropriate. 

# A character object or a numeric vector object. If a character
#'   object, must be one of \{"quad", "uni", "exp"\} indicating the distribution
#'   from which the time points are to be calculated. For character input,
#'   input 'nTimes' must also be provided. If a numeric vector, the
#'   time points to be used. If 0 is not the first value, it will be
#'   concatenated by the software.
#
# successful methods return a vector of unique time points that are sorted in 
#   ascending order and the maximum time point, tau
#

## used in class_TimeInfo.R script

## define a generic function in R called .VerifyTimePoints
setGeneric(name = ".VerifyTimePoints",
           def = function(timePoints, nTimes, ...) { 
             standardGeneric(".VerifyTimePoints") 
           })

## function throw error if not used on more specific class 
## we expect the function to specify a distribution to generate time points

# the default method generates an error
setMethod(f = ".VerifyTimePoints",
          signature = c(timePoints = "ANY",
                        nTimes = "ANY"),
          definition = function(timePoints, nTimes, ...) { 
            stop("timePoints input must be one of {'quad', 'uni', 'exp'} ",
                 "or a numeric vector", call. = FALSE)
          })

## define a new method for .VerifyTimePoints acting on objects of class "character"

setMethod(f = ".VerifyTimePoints",
          signature = c(timePoints = "character",
                        nTimes = "numeric"),
          definition = function(timePoints, nTimes, ..., tau, response) { 
            
            # if timePoints is provided as a character, the character indicates
            # which distribution should be used to obtain the time points from
            # the response data. Current options are {quad, uni, exp} for
            # quadratic, uniform, and exponential, respectively.
            
            ## convert input nTimes to an integer
            
            nTimes <- as.integer(x = nTimes)
            
            ## check if it's positive. If not, stop running and throw error
            
            if (nTimes <= 0L) stop("nTimes must be positive", call. = FALSE)
            
            # if timePoints are generated internally, must use one
            # of "quad", "uni" or "exp"
            
            ## convert timePoints to lowercase
            
            timePoints <- tolower(x = timePoints)
            
            ## check if timePoints matches one of the distribution times
            if (!timePoints %in% c("quad", "uni", "exp")) {
              
              ## if not, calls .VerifyTimePoints with null values (which I think throws error)
              .VerifyTimePoints(timePoints = NULL, nTimes = NULL)
            }
            
            ## check if "response" (input) is a matrix
            
            # identify the 75% percentile
            if (is.matrix(x = response)) {
              
              ## if so, calculate the row sums to determine the times
              
              times <- rowSums(x = response, na.rm = TRUE)
              
              ## if it's not a matrix, 
            } else {
              
              ## use the "response" as observed times
              times <- response
            }
            
            ## check if "tau" (input) is provided (not null)
            
            if (is.null(x = tau)) {
              
              ## if it's null, return this message 
              message("the study length (tau) was not provided ",
                      "and was set as the third quartile of the observed times")
              
              ## also set the value of the third quartile to the observed times
              
              timeThirdQuartile <- stats::quantile(x = times, probs = 0.75)
            } else {
              
              ## otherwise if tau is provided, set that as the third quartile of the observed times
              timeThirdQuartile <- tau
            }
            
            ## if the time point is a "quad" distribution, 
            
            if (timePoints == "quad") {
              
              # quadratically spaced times points between 0 and ceiling of square root of third quartile
              ## create a sequence and square each element to achieve quadratic spacing
              timePoints <- seq(from = 0.0, 
                                to = ceiling(x = sqrt(x = timeThirdQuartile)), 
                                length.out = nTimes)^2
              
              ## if the time point is a uniform distribution, generate points based on third quartile value
              
            } else if (timePoints == "uni") {
              
              ## and third quartile > 1
              
              if (timeThirdQuartile > 1.0) {
                # evenly spaced time points between 0 and ceiling third quartile, rounded up to nearest integer
                timePoints <- seq(from = 0.0, 
                                  to = ceiling(x = timeThirdQuartile), 
                                  length.out = nTimes)
              } else {
                
                ## if third quartile < 1
                # evenly spaced time points between 0 and rounded third quartile
                timePoints <- seq(from = 0.0, 
                                  to = round(x = timeThirdQuartile, digits = 4L), 
                                  length.out = nTimes)
              }
              
              ## if the time point is an exponential distribution, generate times based on third quartile value
              
            } else if (timePoints == "exp") {
              
              ## and third quartile > 1
              
              if (timeThirdQuartile > 1.0) {
                # exponentially spaced time points between 0 and ceiling of the log of third quartile + 1
                ## then subtract 1 from each element of the exponential sequence
                timePoints <- exp(x = seq(from = 0, 
                                          to = log(x = ceiling(x = timeThirdQuartile) + 1L), 
                                          length.out = nTimes)) - 1.0
                
              } else {
                
                ## here, don't need to subtract 1 since the value of timeThirdQuartile < 1
                # evenly spaced time points between 0 and rounded third quartile
                timePoints <- exp(x = seq(from = 0, 
                                          to = log(x = round(x = timeThirdQuartile, digits = 4L) + 1L), 
                                          length.out = nTimes))
              }
            }
            
            
            ## if the minimum time point is significantly greater than 0, 
            
            if (min(timePoints) > 1e-8) {
              
              ## make sure it starts with 0 by adding a 0
              timePoints <- c(0.0, timePoints)
            }
            
            ## if "tau" (input) is null,
            
            if (is.null(x = tau)) {
              
              ## set the value of "tau" to be the maximum number of time points
              tau <- max(timePoints)
              
              ## if the value of "tau" is specified and is greater than the max of timePoints
            } else if (tau > max(timePoints)) {
              
              ## include tau with the timePoints
              timePoints <- c(timePoints, tau)
              
              ## if tau is smaller than the max, but greater than the min,
              
            } else if (tau < max(timePoints) && tau > min(timePoints)) {
              
              ## removes all elements in timePoints greater than tau
              timePoints <- timePoints[timePoints < tau]
              
              ## appends tau to the timePoints
              timePoints <- c(timePoints, tau)
              
              ## if tau is smaller than minimum value, stops execution with error
            } else if (tau < min(timePoints)) {
              stop("tau cannot be < all time points", call. = FALSE)
            }
            
            ## return a list containing the timePoints and tau
            
            return( list("timePoints" = timePoints, "tau" = tau) )
            
          })



## define a new method called .VerifyTimePoints operating when timePoints is a numeric and nTimes is of any time 
setMethod(f = ".VerifyTimePoints",
          signature = c(timePoints = "numeric",
                        nTimes = "ANY"),
          definition = function(timePoints, nTimes, ..., tau, response) { 
            
            ## if timePoints has 0 length, stops execution with an error
            
            if (length(x = timePoints) == 0L) {
              stop("timePoints is of zero length", call. = FALSE)
            }
            
            ## removes any duplicate values and sorts
            
            # if timePoints provided, sort them and ensure uniqueness of values
            timePoints <- sort(x = unique(x = timePoints))
            
            ## if the minimum of timePoints is greater than a very small positive value,
            
            if (min(timePoints) > 1e-8) {
              
              ## add 0 to the beginning to ensure the sequence starts with 0
              timePoints <- c(0.0, timePoints)
            }
            
            ## if "tau" (input) is null,
            
            if (is.null(x = tau)) {
              
              ## set the value of "tau" to be the maximum number of time points
              
              tau <- max(timePoints)
              
              ## if the value of "tau" is specified and is greater than the max of timePoints
              
            } else if (tau > max(timePoints)) {
              
              ## include tau with the timePoints
              timePoints <- c(timePoints, tau)
              
              ## if tau is smaller than the max, but greater than the min,
              
            } else if (tau < max(timePoints) && tau > min(timePoints)) {
              
              ## removes all elements in timePoints greater than tau
              timePoints <- timePoints[timePoints < tau]
              
              ## appends tau to the timePoints
              timePoints <- c(timePoints, tau)
              
              ## if tau is smaller than minimum value, stops execution with error
            } else if (tau < min(timePoints)) {
              stop("tau cannot be < all time points", call. = FALSE)
            }
            
            
            ## return a list containing the timePoints and tau
            return( list("timePoints" = timePoints, "tau" = tau) )
          })
