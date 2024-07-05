

# @param timePoints A numeric vector. The timepoints used in the analysis
#
# @param survVector A numeric vector. The survival function for a single
#   individual
#
# @param by A numeric. The amount of time by which the survival function is
#   to be shifted
#

## shifts a given survival function by a specified amount of time and returns adjusted survival probabilities
## ex) if survival probability at time 5 was 0.8, time 10 was 0.6,and we shift by -5
## then the shifted probabilities are now 0 and 5
## at 0, this is before the first original time point, so survival prob set to 1
## at 5, this is the same as the first original time, so survival prob set to 0.8 which was the original one at that time
## if the shifted time is in the middle of two, uses linear interpolation to return the resulting survival time

# Generate a sample dataset with uniform time points
#set.seed(123)
#timePoints <- seq(1, 15)  # 15 uniform time points
#survVector <- pmin(1, exp(-0.2 * timePoints))  # Exponential survival function
### we shift by observed survival time for each individual from response[elig]
#by <- 3
### we would do these for a single shift AKA use only that patient's "by"


# Function returns the shifted survival function
.shift <- function(

  ## numeric vector. time points used in the analysis
  timePoints,

  ## numeric vector. Survival function for a single individual
  survVector,

  ## numeric. Amount of time by which survival function is to be shifted
  by) {

  # the number of timepoints in the analysis
  nTimes <- length(x = timePoints)

  # shift the timepoints by the specified amount
  timePrime <- timePoints - by

  ## initialize numeric vector "survPrime" to store shifted probabilities with the same length as the original number of points

  # initialize the shifted survival function
  survPrime <- numeric(length = nTimes)

  ## loop through each time point starting with the first

  for (i in 1L:nTimes) {

    ## calculate the shift index-- AKA determine index of the last original time point <= shifted time point[i]
    ## sum these up

    # for each shifted time value, count the number of timepoints <= timeprime_i
    sIndex <- sum(timePoints <= timePrime[i])

    ## if the index of the last original time point isn't the last (nTimes) and is not before the first time point

    if (sIndex < nTimes && sIndex > 0L) {
      # if the number of timepoints <= timeprime_i is not all or none
      # determine the fraction of dt[i] to which the shift corresponds

      ## calculate the fraction of the interval between the two time points
      ## AKA quantifies where between the two time points the shifted time falls as a fraction of their interval

      sFraction <- {timePrime[i] - timePoints[sIndex]} /
        {timePoints[sIndex + 1L] - timePoints[sIndex]}

      # interpolate the survival function to the shifted time value

      ## use linear interpolation so that the survival probability at shifted time point is the weighted average of the probabilities
      ## weights are determined by how close the shifted time is to each of the original time points
      survPrime[i] <- survVector[sIndex] * (1.0-sFraction) +
        survVector[sIndex+1L] * sFraction

      ## if there are no original time points before the shifted point (shifted time is before first time point),

    } else if (sIndex == 0L) {

      ## this means that at the first shifted time point there is full survival

      # if no timepoints are <= to the shifted time, survival set to 1
      survPrime[i] <- 1.0

      ## if all time points are <= the shifted times, this is not possible so throw an error

    } else if (sIndex == nTimes) {
      # it should never happen that all timepoints are <= shifted times
      stop("this condition should not happen -- contact maintainer")
    }

    ## ensures survival function is non-increasing. Survival probabilities should not increase as time progresses, they should stay the same or decrease
    ## if i > 1 (only enforce this check after the first time point bc no previous time to compare against)
    ## and if the survival probability at the current shifted time point is greater than at the previous time point, this is a violation
    ## if there is an increase, we correct the increase by setting the current time point equal to the previous time point

    if (i > 1L && survPrime[i] > survPrime[i-1L]) survPrime[i] <- survPrime[i-1L]

  }

  return( survPrime )
}


## function to shift survival functions for multiple individuals simultaneously and convert these into probability mass vectors
.shiftMat <- function(timePoints,

                      ## matrix where each column represents survival function for an individual
                      survMatrix,

                      ## vector specifying how much to shift survival function for each individual
                      shiftVector,

                      ## boolean: whether to convert shifted survival functions into probability mass vectors
                      surv2prob) {

  # number of timepoints
  nTimes <- length(x = timePoints)

  # number of individuals in subset
  nSamples <- length(x = shiftVector)

  ## initializes matrix to hold shifted survival functions for all individuals filled with 0s at first

  # default survival function to zero
  survShifted <- matrix(data = 0.0, nrow = nTimes, ncol = nSamples)

  ## loop through each of the individuals

  for (i in 1L:nSamples) {
    # for each individual, shift their survival function down in time by
    # the specified amount given in shiftVector
    survShifted[,i] <- .shift(timePoints = timePoints,
                              survVector = survMatrix[,i],
                              by = shiftVector[i])
  }

  ## if the shifted survival functions should be converted into probability mass vectors (surv2prob = TRUE)
  ## turns the continuous survival function (probability of surviving up to each time point) into discrete distribution
  ## probability of event occurring within each time interval.

  if (surv2prob) {
    # if survival function is to be converted to probability mass vector
    # return the change in the survival function at each time step

    ## calculate the change in survival probability at each consecutive pait of time points
    ## subtracting each row of survShifted from the row above it
    ## append a row of 0s at the end to align with matrix dimensions
    ## probability mass vector representing the change in survival probabilities at each time point
    ## NOTE: we coerce survShifted[-1L] to be a matrix in case we have only 1 observation to make stub/double stub
    survShifted <- survShifted - rbind(as.matrix(survShifted[-1L,]), 0.0)
  }

  return( survShifted )
}




