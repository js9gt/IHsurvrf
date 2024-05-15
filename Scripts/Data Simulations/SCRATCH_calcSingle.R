


#!! create a new subroutine called "kaplan" for the KM estimator
#
#! Calculate the Kaplan Meier estimator
#! ns integer, the number of time points
#! nj real(:), at risk
#! oj real(:), events
#! z real(:), estimator

kaplan <- function(ns, nj, oj) {
  # Initialize array to store KM estimator
  z <- numeric(ns)

  # Calculate difference between number at risk and number of events
  num <- nj - oj

  # Calculate the KM estimator for the first time point and store it
  z[1] <- num[1] / nj[1]

  # Loop through each time point
  for (i in 2:ns) {
    # If the number at risk is greater than 0
    if (nj[i] > 1e-8) {
      # Calculate the KM estimator based on difference between number at risk and number of events
      z[i] <- (num[i] / nj[i]) * z[i - 1]
    } else {
      # Otherwise, set the KM estimator to be the estimate at the previous value
      z[i] <- z[i - 1]
    }
  }

  # Return the KM estimator
  return(z)
}


#! estimate the survival function and mean survival time
#!   nCases: integer, the number of elements in casesIn
#!   casesIn: integer, the indices of the subset for which the value is
#!     calculated

casesIn <- c(13, 14, 2, 5, 10)
nCases <- length(casesIn)

calcValueSingle <- function(nCases, casesIn, survFunc, mean) {
  # Initialize variables
  ## surv func is rep (0, nTimes)
  survFunc <- rep(0, nTimes)
  mean <- 0

  # used to calculate number of at risk cases at each time point
  Rb <- rowSums(pr[,casesIn])

  # Initialize Nj(1) with the number of cases at the first time point
  Nj <- numeric(nTimes)
  Nj[1] <- nCases

  # Calculate number of cases at risk at each subsequent time point

  ### 1:26 are fine, then 27 is negative
  for (i in 2:nTimes) {
    Nj[i] <- Nj[i - 1] - Rb[i - 1]
  }

  # Number of events at each time point
  Oj <- numeric(nTimes)
  for (i in 1:nTimes) {
    # Calculate the number of events
    Oj[i] <- sum(pr[i, casesIn] * delta[casesIn])
  }

  # Kaplan-Meier estimate survival function
  survFunc <- kaplan(ns = nTimes, nj = Nj, oj = Oj)

  # Mean survival time
  mean <- sum(survFunc * dt)

  # Return values
  return(list(survFunc = survFunc, mean = mean))
}
