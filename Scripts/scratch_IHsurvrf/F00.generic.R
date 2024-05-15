

library(dplyr)
search2 <- function(s, xvec, include = TRUE) {
  # Finding the index of the largest value inxvec <= input value of"s"
  # if include = TRUE, consider values = s
  # if include = FALSE only consider strict inequality
  # Fn is the cdf table (x (should be ordered) and y)
  # index = which (s >= c(-Inf, xvec) & s <= c(xvec, Inf)) - 1
  index = if (include) {which (s >= c(-Inf, xvec)) - 1} else {which (s > c(-Inf, xvec)) - 1}
  
  # returns index for largest value in xvec <= to "s"
  max(index)
}

### test : search2(49.25, xvec =  48:52)

# vectorizes the function above to return the index for the largest value in xvec <= "s"
## in this case, s is turned into a vector
## search2.vec can be be applied to each element of vector "s"
## returns vector of indices for largest values in xvec <= each element of "s", a vector
search2.vec <- Vectorize(search2, vectorize.args = "s")

#################
##             ##
#################

## linear interpolation: curve fitting using linear polynomials to construct new datapoints
##### from discrete set of known data points
##### this function is designed for interpolation between 2 points with x-coordinates
##### can also return interpolated y-values or exponential tail for y-values

# s: value for which interpolation is sought
# xvec: vector containing x coordinates
# yvec: OPTIONAL, vector containing y-coordinates. If provided, interpolation will be done on y-vals too
# return.y.only: if TRUE, function only returns interpolated y-value
# exponential.tail: if TRUE, considers exponential tails in y-vector

lin.interpolate <- function(s, xvec, yvec = NULL, return.y.only = FALSE, exponential.tail = TRUE) {
  
  # find length of x-vector
  x.len = length(xvec)
  
  # recall: search2 returns index for the largest value in xvec <= to "s" which is a scalar
  low.index <- search2(s, xvec)
  
  # get value of xvec corresponding to largest value <= s
  lower <- xvec[low.index]
  
  # check if this index is at the end of xvec (which is ordered)
  if (low.index == x.len) {
    
    # if it's at the end, we rewrite this as the upper value & have proportion 0
    upper <- lower
    proportion <- 0
  } else {
    
    # otherwise, the upper x-value is the next higher value
    upper <- xvec[low.index + 1]
    
    # proportion between the lower and upper values
    proportion <- (s - lower) / (upper - lower)
  }
  
  # create a result with the index of the largest xvec <=s, and the proportion
  result <- c(low.index = low.index, proportion = proportion)
  
  # if the y vector is provided as an input, perform interpolation for the y-values
  if (!is.null(yvec)) {
    lower.y <- yvec[low.index]
    upper.y <- yvec[min(x.len, low.index + 1)]
    
    
    if (is.infinite(upper.y) && exponential.tail) {
      # method of calculating proportion for infinite upper value of y and if the tail is exponential 
      y.interpolate <- lower.y * log (1 - s) / log (1 - lower)
    } else {
      y.interpolate <- lower.y + proportion * (upper.y - lower.y)
    }
    
    #if return.y.only = TRUE, only return interpolated y-value
    if (return.y.only) return(y.interpolate)
    
    # add interpolated y to the result vector
    result["y.interpolate"] <- y.interpolate
  }
  # print the result vector
  return(result)
}

# vectorizes the linear interpolation function which can accept a vector of values s
# performs interpolation on each element of s simultaneously
# returns different interpolation results for dif values of s
lin.interpolate.vec <- Vectorize(lin.interpolate, vectorize.args = "s")

#################
##             ##
#################

# time calculator
# s:takes integer values to denote dif operations
# reset: TRUE means reset time tracking
# units: species units for time calculation
## tt(1): records start time
## tt(2): records elapsed time
## tt(2, reset = TRUE): resets timer for later measurements

tt <- function(s, reset = FALSE, units = "auto"){
  # if s = 1, records current system time & assigns to variable time.tmp in GLOBAL ENVIRONMENT
  if (s==1) {time.tmp <<- Sys.time() # record time
  
  # if s = 2, calculates elapsed time since time when s = 1
  } else if (s==2) { # calculate time
    
    # begin: time when s was 1
    result <- data.frame(begin = time.tmp, 
                         # end: current time
                         end = Sys.time(), 
                         # elapsed time, formatted according to the specified units
                         elapsed = difftime(Sys.time(), time.tmp, units = units))
    
    # if reset = TRUE, when s = 2, resets time.tmp to current time for subsequent tracking
    if (reset) time.tmp <<- Sys.time()
    return(result)
  }
}

#################
##             ##
#################

## checks whether input value (x) is NA or 0
# either 0 or na gives TRUE logical output
# used in G01.Goldberg.R::auxiliary()
# is.na0(c(1:3,NA, 0, 0, 3))
is.na0 <- function(x) {
  is.na(x) | (x == 0)
}

#################
##             ##
#################


# complete the delta vector so that the final status is inherited to the next NA values
## AKA, ensures any NA values in "vec" are replaced with last known non-NA

# vec:intput vector with NA values
# nstages: OPTIONAL param, number of stages or by default, length of the vector

#complete.delta.vec <- function(vec, n.stages = length(vec)) {
#  
#  # identifies index of first NA values in the vector, or if all values are non-NA sets to n.stages + 1
#  # if all elements are non-NA, ensures "ref" exceeds the vector length by 1
#  # -1 applied to ensure index refers to last non-NA value or last position of vector
#  
#  ref <- min(which(is.na(vec)), n.stages + 1) - 1
#  
#  # stores last non-NA value 
#  delta.final <- vec[ref]
#  
#  # if reference index less than total number of stages, there are subsequent stages after first NA value
#  if (ref < n.stages)
#    # then,replace all values from reference index to end of vector with the last non-NA value
#  vec[ref:n.stages] <- delta.final
#  
#  # returns the updated vector with NA values filled in
#  vec
#}

## Example: complete.delta.vec(c(1,0,NA, 0))

#################
##             ##
#################

# operates on a matrix with missing values along its rows
# uses the complete.delta.vec from before to fill missing values in within each row of the matrix

# mat: input matrix with missing values in ows
# n.stages: OPTIONAL param determining number of stages (columns). By default, set fo # of columns in matrix

#complete.delta.mat <- function(mat, n.stages = dim(mat)[2]) {
#  
#  # apply the complete.delta.vec to each row of the matrix
#  # to fill in missing values for each row
#  # then transpose the matrix to return the matrix where missing vals in each row are replaced
#  apply(mat, 1, complete.delta.vec) %>% t
#}
## tmp[1:10,"delta",]
# complete.delta.mat(tmp[1:10,"delta",])


#################
##             ##
#################

# calculates survival probabilities at specific time point based on time, t & S
# t: time point at which survival prob needs to be calculated
# S: survival object containing info about survival probabilities 
# rightcts: TRUE = includes the right censored time points

S.t <- function(t, S, rightcts = TRUE) {
  # uses search2.vec to find indices for largest values in S$time vector <= each element of "t", a vector
  # AKA find largest survival time up to time t
  index = search2.vec(t, S$time, include = rightcts)
  
  # retrieves survival probability from S$surv vector for the closest time point up to time t
  # pmax() used to prevent invalid indices bc if index contains value < 1
    # then return vector where each element is maximum of the corresponding value from index & 1
  S$surv[pmax(index, 1)]
}

#################
##             ##
#################

# generates a summary dataframe from an "outputData" object

flowchart <- function(outputData) {
  n.stages = dim(outputData)[3]
  n.sample = dim(outputData)[1]
  
  # initializes df called "result" with columns for stage numbers, total count, percentage, 
          # censored count, deceased count, died count, count for next stage
  result <- data.frame(stage = c(1:n.stages, "Total", "Percent"),
                       total = rep(NA, n.stages + 2), percentage = NA,
                       censored = NA, died = NA, nextStage = NA)
  
  # loop through stages. For each stage,
    # calculate total count for that stage based on number of rows in "outputData"
    # compute percentage of samples in that stage compared to total SS
    # identifies censored cases (delta = 0)
    # identifies failure/death (gamma = 1) (Tk < Uk)
    # updates "outputData" to exclude censored & failure cases
    # updates "result" DF with counts of censore, failures, count for next stage
    # gets total aggregated counts for total and percent
  for (stage in 1:n.stages) {
    
    # calculate total count for each stage
    result[stage, "total"]  = dim(outputData)[1]
    
    # calculate percentage for each stage
    result[stage, "percentage"]  = round(result[stage, "total"] / n.sample,2)
    
    # Identify censored & deceased cases for the current stage
    cens <- outputData[, "delta", stage] == 0
    # drop those censored from the outputData
    outputData <- outputData[!cens,,, drop = FALSE]
    died <- outputData[, "gamma", stage] == 1
    # drop those died from the outputData
    outputData <- outputData[!died,,, drop = FALSE]
    
    # update results DF with the counts of censored, deceased, and next stage
    result[stage, "censored"]  = sum(cens, na.rm = TRUE) 
    result[stage, "died"]      = sum(died, na.rm = TRUE)
    result[stage, "nextStage"] = sum(!died, na.rm = TRUE)
  }
  
  # calculate total censored & deceased counts across all stages
  result[n.stages + 1, c("censored", "died")] <- apply(result[, c("censored", "died")], 2, sum, na.rm = TRUE)
  
  # calculate percentage of censored & deceased cases across all stages
  result[n.stages + 2, c("censored", "died")] <- round(result[n.stages + 1, c("censored", "died")]/n.sample, 2)
  # result[c("Total", "Percent"), 1:2]
  
  # result of all summarized information in each stage & overall totals
  result
}

############################# Example of using flowchart #################
# Simulated data generation
#set.seed(123)
#n_stages <- 5  # Number of stages
#n_patients <- 100  # Number of patients

# Creating empty output data structure
#outputData <- array(NA, dim = c(n_patients, 2, n_stages),
                    #dimnames = list(paste0("Patient_", 1:n_patients),
                                     #c("delta", "gamma"),
                                    #paste0("Stage_", 1:n_stages)))


# Simulating random data for delta and gamma for each stage
#for (stage in 1:n_stages) {
  #outputData[, 1, stage] <- sample(c(0, 1), n_patients, replace = TRUE)  # Simulating 'delta'
  #outputData[, 2, stage] <- sample(c(0, 1), n_patients, replace = TRUE)  # Simulating 'gamma'
#}

# Applying the flowchart function on the simulated data
#result_summary <- flowchart(outputData)
#result_summary


###########################################################################


#################
##             ##
#################

## Convert multi-stage dataset, "output" into observable format for analysis

output2observable <- function(output, stage = NULL, cumulative = FALSE) {
  
  # gets number of stages from the 3rd dimension of the "output" array
  
  n.stages = dim(output)[3]
  
  # if no stage provided, considers all stages
  if (is.null(stage)) {
    stage = 1:n.stages
    
    # if cumulative = TRUE, sets stage to all stages up to the maximum
  } else if (cumulative) {
    stage = 1:max(stage)
  }
  
  # number of rows in output dataset (n)
  n = dim(output)[1]
  
  # defining variable names related to covariates & stage variables
  nm = dimnames(output)[[2]]
  nm.covar = c("lB", grep("Z[0-9]+", nm, value = TRUE))
  nm.stage =  c("event.time", "delta", "action", nm.covar)
  
  # initializes empty DF with subject ID column for the first stage, taken from "output"
  df <- data.frame(subject.id = output[, "subj.id", 1])
  
  # iterates through each stage to make
    # df.i: DF containing columns related to Time, delta, action, covariates for stage i
    # mergecolumns of this into main dataframe (DF)
  for (i in stage) {
    df.i <- data.frame(output[, nm.stage, i])
    names(df.i) <- paste(c("T", "delta", "A", nm.covar), i, sep = "_")
    df <- cbind(df, df.i)
  }
  
#  # calculates new delta value in DF
#  df$delta = 
#    # selects columns starting with "delta_"
#    df %>% dplyr::select(starts_with("delta_")) %>% 
#    # converts selected cols into matrix
#    as.matrix %>% 
#    # applies function to each row to see if any values == 0 (Censoring)
#      # this returns a 1 if there's no 0, returns a 0 if there is a 0
#    apply(1, function(s) 1 - any(s == 0, na.rm = TRUE))
  df
}

# output2observable(obs.data$output)



#################
##             ##
#################

# find index of maximum value within a vector, x
# handles scenarios where there might be ties for the maximum value
# random: selects one of indices randomly for the set of tied indices

whichMax <- function (x, tie.method = "random") {
  
  # indentifies indices where value in x is equal to the maximum
  ind = which(x == max(x))
  
  # if there's only one index, return that one
  if (length(ind) == 1L) return(ind)
  
  # if there are multiple indices (ties)
  if (tie.method == "random") {
    sample(ind, 1)
    
    # if the tie method = "first" then return index of the first occurence of the max value
  } else if (tie.method == "first") {
    ind[1]
    
    # if tie method isn't "random" or "first", return an NA
  } else {
    NA
  }
}