

source("R/class_IH.pool1.R")



### initialize a list of results which hold:
## the final data
## the estimated survival probabilities
## the results of each of the two forests for replication:
##### pool1_results: the results of the pooled data for ALL prior data points (excluding current iteration of loop)
##### append1_results: results including newly appended data (current iteration of loop)-- but this doesn't include fitting whole forest w newly appended data


## define an internal function .IHsurvstep
.IHsurvstep <- function(...,
                         ## input arguments
                        data,
                         models,
                         params,
                        ## IHsurvresults[[q + 1L]]@PMVoptimal
                         priorIH,
                        ## k is the stage
                        k,
                        ## uses the current iteration counter
                        iter_counter,
                        nDP){

## uses the input k as the stage
k <- k

iter_counter <- iter_counter

nDP <- nDP
## initialize an empty lis to store results of the analysis to hold output from each decision point
IHsurvresults <- list()

pr_poolappend2 <- priorIH@PMVoptimal


## calls the .dtrSurvStep() function for the pooled data, only on stages ndp - 2, ndp -1, ndp
## defined in class_DTRSurvStep.R
## use model from last time point
## returns object of class"DTRSurvStep" with the predicted survival curves for each treatment & also the one that's the optimal (based on whatever criteria of interest)
## input the HC optimal survival probabilities

pool1_results <- .dtrSurvStep(
  ## use the model for pooled data
  model = models,
  ## only include the data where the stage >= nDP - 2, or for the next iteration in stage nDP - 4, use nDP - 3
  data = data[data$stage >= (k + 1), ],
  priorStep = NULL,
  params = params,
  txName = "A",
  ## set to the same as the input ((NOTE THIS GETS CHANGED EARLIER IN THE CODE, so we temporarily set it to null))
  mTry = NULL,
  ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
  sampleSize = 1,
  ## we are in the first step of pooling patient data after running HC code
  pool1 = TRUE,
  appendstep1 = TRUE,

  ## using the full survival probabilities from the last cycle
  ## ex) includes nDP - 3, nDP - 2, nDP - 1; where nDP - 4 is the new data
  inputpr = pr_poolappend2
)

## now, we want to retrieve the optimal actions for each patient for stages nDP - 2, nDP - 1, nDP from the first pooling
# pool1_results@optimal@optimalTx

# Step 1: Create a new column "A.pool1" in long_data
# Generate column name dynamically based on iteration counter
column_name1 <- paste("A.pool", iter_counter, sep = "")

# Set the column name in long_data
data[[column_name1]] <- NA

# Step 2: Identify indices for all the stages we just used (which don't include the current stage bc we haven't used that for appending yet)
## meaning, if nDP - 4 is the new data and current loop, we want nDP - 3, nDP - 2, nDP - 1, nDP
stage_indices <- which(data$stage %in% (k + 1):nDP)

# Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
## also update the "A" column
data[[column_name1]][stage_indices] <- pool1_results@optimal@optimalTx
data$A[stage_indices] <- pool1_results@optimal@optimalTx

## retrieve matrix of optimal only for stage nDP - 3:
## recall: we can start from seq 1 since we haven't yet incorporated the new data
## current is nDP - 4 loop...
# Calculate the indices of the columns you want to retrieve-- note this matrix needs to be transposed to get to our dimensions
## when retrieved for some reason it has patients for rows, instead of columns
pool1pr_indices <- seq(1, ncol(t(pool1_results@optimal@optimalY)), by = (nDP - k))


# Retrieve the columns of optimal survival probabilities using the calculated indices
## this is equivalent to class_IH.DTRSurvStep.R: t(.OptimalY(object = priorStep@optimal))

shiftedprob1 <- t(pool1_results@optimal@optimalY)[, pool1pr_indices]


#########################################################
################# append 1 alg ############################
#########################################################

### to do: only do this for eligible patients for the new stage nDP - 3, and then perform shifting

# identify individuals with complete data in new stage nDP - 3

## prepares model frame using the formula (input) and the data (input)-- NOTE we use the original model input
## models = Surv(Y,D)~X + A
## missing values are not specifically handled (not omitted)
## only include the variables that are used in the model formula mod
## If there are variables in your data that are not used in the formula, they won't be included in the resulting x.

x_append1 <- stats::model.frame(formula = models,
                                ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                data = data %>% filter(stage == k) %>% select(-matches("^A\\.")),
                                na.action = na.pass)


## we want to exclude the A.opt.HC column and A.pool1 column
elig_append1 <- stats::complete.cases(x_append1)

# extract response and delta from model frame

## extract survival response
response_append1 <- stats::model.response(data = x_append1)

## extract censoring indicator (delta) from the second column of the "response" data
## "L" is used to indicate that 2 is an integer
delta_append1 <- response_append1[, 2L]

## updates the "response" variable to only include the first column of the original "response" data which represents survival times
response_append1 <- response_append1[, 1L]

# remove response from x

## if first column of the model frame (x) is the response variable, remove this column
## probablhy to construct predicte response from the predictors, since the response has nothing to do with the prediction itself
if (attr(x = terms(x = models), which = "response") == 1L) {
  x_append1 <- x_append1[,-1L, drop = FALSE]
}

## marks zeroed survival times and updates eligibility

# responses that are zero (effectively) indicate censored at a previous stage

## 1e-8 is the tolerance, if these responses are smaller than a very small number, this is marked as TRUE
zeroed <- abs(x = response_append1) < 1e-8

## update eligibiity vector so that cases that are eligible can't have been marked as zeroed

elig_append1 <- elig_append1 & !zeroed

## if there are no eligible cases, output an error message

if (sum(elig_append1) == 0L)
  stop("no cases have complete data", call. = FALSE)

## displays message indicating the number of cases that are still eligible for analysis at this stage

message("cases in stage: ", sum(elig_append1))


## retrieve the number of timepoints from the params object
## defined in class_TimeInfo.R

nTimes <- .NTimes(object = params)

# create an empty matrix for all previously eligible cases

## initializes a survival matrix called "survMatrix" with 0's

survMatrix <- matrix(data = 0.0,
                     nrow = nTimes,
                     ncol = nrow(x = x_append1))

## sets the first row to 1.0

survMatrix[1L,] <- 1.0

## updates survMatrix with survival functions estimated from the previous step
## priorStep: A DTRSurvStep object. The analysis from a previous step

# retrieve estimated OPTIMAL survival function from previous step
## then accesses the "eligibility" slot (logical) to select only the columns where the patient was still eligible from the previous time
## replace these with the estimated optimal value from the"optimal" slot of the "priorStep" object
## transpose the output of .OptimalY to align with the structure
## .OptimalY defined in class_Optimal.R
## .OptimalY acts on optimal slot of priorStep (class DTRSurvStep)-- this comes from the .PredictAll()
## optimal slot is  object of class optimal
## .OptimalY acts to return the optimalY slot

## Currently, we are inloop nDP - 4 is the new data-- we have not incorporated this yet
## we select the eligible patients in the stage prior: stage nDP + 1 (nDp - 3)
priorStep_elig <- pool1_results@eligibility[seq(1, length(pool1_results@eligibility), by = nDP - (k))]

## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
survMatrix[, priorStep_elig] <- shiftedprob1


# shift the survival function down in time (T_i - Tq) and
# transform to a probability mass vector for only those
# eligible for this stage
# .shiftMat is an internal function defined in shiftMat.R

## transforms the survival functions in survMatrix to probability mass format for the current stage
## shifts the survival function based on the observed survival times
## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage

append1_pr <- .shiftMat(
  timePoints = .TimePoints(object = params),

  ## extracts columns from survMatrix corresponding to cases that are eligible
  ## this is a matrix matrix where each column represents survival function for an individual
  survMatrix = survMatrix[, elig_append1, drop = FALSE],

  ## extracts survival times corresponding to eligible cases
  ## this is how much to shift survival function for each individual
  shiftVector = response_append1[elig_append1],

  ## probably transforming survival times into probabilities?
  surv2prob = TRUE
)

## sets very small values in pr to 0

append1_pr[abs(append1_pr) < 1e-8] <- 0.0

##### NOTE: this step is equivalent to shifting by 0, essentially we just transform the optimal survival probabilities into probability mass vectors
## this step is essential to match with the probabilities that we've just computed
####################
####################
####################
####################
prev.op1 <- t(pool1_results@optimal@optimalY)
## calculate the change in survival probability at each consecutive pait of time points
## subtracting each row of survShifted from the row above it
## append a row of 0s at the end to align with matrix dimensions
## probability mass vector representing the change in survival probabilities at each time point
prev.op1 <- prev.op1 - rbind(prev.op1[-1L,], 0.0)
####################
####################
####################
####################


## now, we need to add append1_pr into the pool1_results@optimal as an input-- adding the shifted optimals basically into the pooled optimals

# Initialize a matrix to store the combined results, with a matrix of NAs
pr_poolappend <- matrix(NA, nrow = nTimes, ncol = ncol(prev.op1) + ncol(append1_pr))

# Define the number of columns to insert from append1_pr
num_cols_insert <- ncol(append1_pr)

# Define the index to insert columns from append1_pr
insert_index <- seq(from = 1, to = ncol(pr_poolappend), by = nDP - (k) + 1)

# Initialize a counter for columns from append1_pr
col_counter <- 1

# Loop over each index to insert columns from append1_pr into combined_results
for (i in seq_along(insert_index)) {
  # Insert columns from append1_pr into combined_results
  pr_poolappend[, (insert_index[i]):(insert_index[i] + (nDP - k))] <-
    cbind(append1_pr[, col_counter], prev.op1[, ((i - 1) * (nDP - k) + 1):(i * (nDP - k))])

  # Increment the counter for columns from append1_pr
  col_counter <- col_counter + 1
}


## after appending, we fit another pooled random forest, including the data from nDP - 3

append1_results <- .dtrSurvStep(
  ## use the model for pooled data
  model = models,
  ## only include the data where the stage >= the current stage in the loop which includes the new data
  ## so, our new data is nDP - 4, which is the current stage that the loop is in
  data = data[data$stage >= k, ],
  priorStep = NULL,
  params = params,
  txName = "A",
  ## set to the same as the input
  mTry = NULL,
  ## NOTE: sampleSize is a vector of inputs from 0-1, but we use the whole sample size so we just specify 1
  sampleSize = 1,
  ## this is used as TRUE for processing the treatment levels
  pool1 = TRUE,
  appendstep1 = TRUE,
  inputpr = pr_poolappend
)

## retrieve the optimal treatments: these are for all 400 "patients" with 4 stages-- we add this to a column called "A_append1"
## we want to only retrieve the optimal for stages nDP - 2, nDP - 1, nDP as well as the optimal survival curves
## therefore, starting from index 2, 6, 10, etc


column_name2 <- paste("A.append", iter_counter, sep = "")

# Set the column name in data
data[[column_name2]] <- NA

# Step 2: Identify indices for stages 3, 4, and 5
## NOTE: currently we are at stage nDP - 4
## this means we only want to select the treatments we've already used for pooling AKA not the current stage
stage_indices <- which(data$stage %in% (k + 1):nDP)

# Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5-- we want to NOT select column 1, 5, etc.
### AKA we select all indices except for these
## also update the "A" column


data[[column_name2]][stage_indices] <- append1_results@optimal@optimalTx[!(1:length(append1_results@optimal@optimalTx) %% (nDP - (k) + 1) %in% 1)]
data$A[stage_indices] <- append1_results@optimal@optimalTx[!(1:length(append1_results@optimal@optimalTx) %% (nDP - (k) + 1) %in% 1)]

## retrieve the optimal probabilities again at the last stage (nDP - 2)
# Retrieve the columns of optimal survival probabilities using the calculated indices

shiftedprob2 <- t(append1_results@optimal@optimalY)[, seq(2, nrow(append1_results@optimal@optimalY), by = (nDP - (k) + 1))]


#########################################################
################# append 2 alg ############################
#########################################################

## we are using the same data to append, so we need to re-initialize the survival matrix

## initializes a survival matrix called "survMatrix" with 0's

survMatrix <- matrix(data = 0.0,
                     nrow = nTimes,
                     ncol = nrow(x = x_append1))

## sets the first row to 1.0

survMatrix[1L,] <- 1.0

## we also re-update the eligibility. We need to select eligibility from stage nDP - 2, which start from column 2, and then goes every nDP - 3
## if we are in stage nDP - 4 as the current iteration, we want to select the ones corresponding to nDP - 3
priorStep_elig <- append1_results@eligibility[seq(2, length(append1_results@eligibility), by = (nDP - (k) + 1))]

## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
survMatrix[, priorStep_elig] <- shiftedprob2


append2_pr <- .shiftMat(
  timePoints = .TimePoints(object = params),

  ## extracts columns from survMatrix corresponding to cases that are eligible
  ## this is a matrix matrix where each column represents survival function for an individual
  survMatrix = survMatrix[, elig_append1, drop = FALSE],

  ## extracts survival times corresponding to eligible cases
  ## this is how much to shift survival function for each individual
  shiftVector = response_append1[elig_append1],

  ## probably transforming survival times into probabilities?
  surv2prob = TRUE
)


## sets very small values in pr to 0
append2_pr[abs(append2_pr) < 1e-8] <- 0.0

### now, we need to insert these into the optimal probabilities from the append1_results random forest
## we replace the estimated optimal probabilities for stage nDP - 3 in the forest with ALL patients:
## 1st column from appended pr, REPLACE 5th column 9th column, etc-- every 4th column should be REPLACED (not inserted) from appended pr

## now, we need to add append1_pr into the pool1_results@optimal as an input-- adding the shifted optimals basically into the pooled optimals
### NOTE: we need to essentially shift the pooled optimal survival function by 0-- convert these into probability mass vectors
## this is essential to match with the shifted probabilities we just calculated

prev.op2 <- t(append1_results@optimal@optimalY)
## calculate the change in survival probability at each consecutive pait of time points
## subtracting each row of survShifted from the row above it
## append a row of 0s at the end to align with matrix dimensions
## probability mass vector representing the change in survival probabilities at each time point
prev.op2 <- prev.op2 - rbind(prev.op2[-1L,], 0.0)
####################
####################
####################
####################



# Initialize a matrix to store the combined results, with a matrix of NAs
pr_poolappend2 <- matrix(NA, nrow = nTimes, ncol = ncol(prev.op2))



# Define the index to insert columns from append1_pr
insert_index <- seq(from = 1, to = ncol(pr_poolappend2), by = (nDP - (k) + 1))

# Exclude columns from pr_poolappend2 that we will be "overwriting"
prev.op_future <- prev.op2[, -insert_index]

# Initialize a counter for columns from append1_pr
col_counter <- 1

# Loop over each index to insert columns from append1_pr into combined_results
## what we will do, is insert the column from append2_pr for the current stage, then the past optimal ones

for (i in seq_along(insert_index)) {
  # Insert columns from append1_pr into combined_results
  pr_poolappend2[, (insert_index[i]):(insert_index[i] + (nDP - k))] <-
    cbind(append2_pr[, col_counter], prev.op_future[, ((i - 1) * (nDP - k) + 1):(i * (nDP - k))])

  # Increment the counter for columns from append1_pr
  col_counter <- col_counter + 1
}


iter_counter <- iter_counter + 1

## create a new object of class "IHsurvrf"

IHsurvresults <- new(
  Class = "IHsurvrf",
  ## data uses input data
  "data" = data,
  ## oprimal probability mass vector from the newly appended data
  "PMVoptimal" = pr_poolappend2,
  "survRF1" = pool1_results,
  "survRF2" = append1_results
)

## return newly created "IHsurvrf" object

return(IHsurvresults)

}


