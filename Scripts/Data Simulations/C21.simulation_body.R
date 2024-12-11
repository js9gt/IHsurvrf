
library(parallel)

library(IHsurvrf)

setwd("~/survrf/Scripts/Data Simulations")
## generic is sourced in multistage_sim.R
source("F02.multistage_sim.R")

## checks value of the variable "criterion" to define value function
if (criterion == "mean") {
  val.fn <- function(x) {mean(x, na.rm = TRUE)}
} 


## initialize matrix called stat.stage  with NA values with a row for each of the simulations run, and a column for each stage
stat.stage <- matrix(NA, nrow = n.sim, ncol = n.stages,
                     dimnames = list(1:n.sim, 1:n.stages))

## initialize a list called arg.obs with various parameters and settings used for analysis
## assigns the same value of the initialization list to all of these different variables:
## used as arguments for observed data, IHsurvrf, observed data with no censoring, and ZOM
## these are to be input into multistage data generation


## we want the no censor part, IHsurvRV, ZOM to have no censoring
arg.obs <- arg.IHsurvrf <- arg.obs.no.censor <- arg.trt1 <- arg.trt0 <- arg.trueopt <-
  list(
    n.sample = n, max_stages = n.stages, tau = tau,
    ## initially all patients are at risk
    at.risk = 1,
    ## Life so far lived at the beginning of the stage
    cumulative.time = 0,
    ## inital length of previous visits
    prior.visit.length = 0,
    rho = 0.5,
    # action-specific effects
    D0 = 0,
    D1 = 1,
    g = 0.25,
    a1 = a1, b1 = b1, z1 = z1, p1 = p1, g1 = g1, h1 = h1, r1 = r1,
    a2 = a2, b2 = b2, z2 = z2, p2 = p2, g2 = g2, h2 = h2, r2 = r2,
    a3 = a3, b3 = b3, z3 = z3, p3 = p3, g3 = g3, h3 = h3, r3 = r3,

    ## a logical for if we want to include censoring (besides administrative censoring)
    ### for no censoring argument, we will change this part to FALSE: for everything except "observed"
      censoringyesno = TRUE,

    ## for inputting dif policies used to generate data
    ## this is a DTRSurv object
    policy = NULL,
    summary = TRUE

  )

## for the observed data, we allow for censoring, everything else doesn't have censoring
arg.obs$censoringyesno <- TRUE


## setting n.sample attribute of ll these objects (created above) to "n.eval"
#### n.eval is a default setting read into C21.simulation_run
#### this is the number of samples per simulation
### Q: why do do we overwrite this from n to n.eval?
## observed data still has the regular number of samples-- I think this is "training data" so we have more samples

## only arg.obs still has the original number of samples
arg.obs.no.censor$n.sample <- arg.IHsurvrf$n.sample <- arg.trt1$n.sample <- arg.trt0$n.sample <- arg.trueopt$n.sample <-
  n.eval


## setting nodesize to 5
nodesize = 3
## setting minimum number of deaths by rounding sqrt of "nodesize" to nearest integer (no decimals)
mindeath = round(sqrt(c(nodesize)), 0)


# flow: obs (1) -> optimal (n.mc), obs.no.censor (n.mc)
print(Sys.time())

# Define the function that runs a single simulation
run_simulation <- function(sim){
  result <- data.frame(no = sim, observed = NA, IHsurvrf = NA, trt1 = NA, trt0 = NA, trueopt = NA, 
                      time.obs = NA,
                       time.IHsurvrf = NA, time.trt1 = NA, time.trt0 = NA, proportion.censor = NA, num.convergence.it = NA)  # time for each method (both policy est and eval)
                       
  
  attr(result, "criterion") <- list(criterion = criterion)
  
  for (i in 1:n.stages) result[[paste0("n_", i)]] = NA
  
  
  cat("###########################  simulation ", sim, "########################### \n")
  cat("########################### (criterion ", criterion,  ")################## \n")

  cat ("1. Data \n")


## tt(1): records start time
## in F00.generic.R
tt(1)

# generate multistage obs.data (allows for censoring) using arg.obs defined earlier
set.seed(sim*10000)
## do.call: calls a function and a list of arguments to be passed to it
## observed data uses "n" number of patients which is the sample size
obs.data <- do.call(simulate_patients, arg.obs)


# observed policy value
##### simulate a set of data again, when there is no censoring (AKA) using a censoring value of -10 for all surv.prev
## also use n.eval as number of samples: 10000
set.seed(sim*10000 + 10)
obs.data.rep <- do.call(simulate_patients, arg.obs.no.censor)

### observed policy estimation:
## val.fn depends on the criterion employed
## observed value for mean = take mean of the total event time (across all patients)
## in this case, the critical value will be null
## this is the observed policy from the data we generated
## if criterion is "surv.prob", "surv.mean", count the proportion of elements >= crit.value (ignoring NAs)
## crit.value is put at 5 (5 year survival rate) in C21.simulation_run.R

## this observed policy is based on n.eval number of reps to calculate a single value of observed policy
## we do this for the number of stages set to ss-- this is already accounted for in the cumulative.event.time
result["observed"] <- val.fn(obs.data.rep$summary$cumulative.event.time)

## note: 0 = censored, 1 = not censored, so we calculate the proportion of 1's from the generated observed data, then subtract from 1
####### USES OBSERVED DATA, not policy data
result["proportion.censor"] <- 1 - mean(obs.data$summary$censor.status, na.rm = TRUE)


# select columns named "n_1", "n_2"... "n_nstages" which we created earlier
result[paste0("n_", 1:n.stages)] <-
  ## for each stage, calculate the proportion of observations in each stage (s = whatever iteration of the n.stages we are in)
  ### AKA this is the proportion of patients who have their terminal stage at stage s
  sapply(1:n.stages, function(s) mean(obs.data$summary$terminal.stage == s))



result["time.obs"] <- tt(2, reset = TRUE)["elapsed"]

print(flowchart(obs.data$output))

# transforming data from an array format to a data.frame format
data.df = output2observable(obs.data$output)

## removes the data generated for the policy (with n.eval iterations) from the R environment to free up memory space,also runs garbage collection to free up memory
rm(obs.data.rep); gc()


#########  IHsurvrf policy estimation ####################
cat ("2. IHsurvrf \n")

## skip.IHsurvrf is set in C21.simulation_run.R (default is TRUE)
## if !skip.IHsurvrf == TRUE, execute the code inside the block (meaning, skip.IHsurvrf = FALSE)
if (!skip.IHsurvrf) {
  cat ("  2. IHsurvrf - Policy estimation \n")

  # Generate the formulas for each stage-- NOTE that this formula is the same across ALL stages
  models = Surv(`T`, delta) ~ A + baseline1 + baseline2 + state1 + state2 + prior.visit.length + cumulative.time + nstages + action.1.count + action.0.count


  ##### CHECK THESE ARGUMENTS

  ## create new list of arguments to use IHsurvrf
  ## used the observed data generated earlier (once it's in observable form)
  arg.IHsurvrf2 = list(data = data.df,
                  ## will give A_1, A_2, ... A_n.stages
                  txName = paste("A", 1:n.stages, sep = "_"),

                  ## list of models generated above for each stage
                  models = models,
                  ## use the observed lengths of the previous stages in the current stage's model as covariates
                  ## we don't need this as our data simulations currently track the cumulative time
                  #usePrevTime = FALSE,
                  ## maximum study length, set as 10 in C21.simulation_run.R
                  tau = tau,
                  ## distribution to draw timepoints from
                  timePoints = "uni",
                  ## number of timepoints from 0 to tau
                  nTimes = tau,

                  ## mean,surv.prob, or surv.mean
                  criticalValue = criterion,

                  ## if the criterion is mean, use truncated mean, otherwise, use logrank test
                  splitRule = "mean",

                  ## use ERT
                  ERT = TRUE,

                  ## If 'ERT' and 'uniformSplit' are TRUE,
                  #'   the random cutoff is sampled from a uniform distribution over the range
                  #'   of available covariate values for tree splitting
                  uniformSplit = TRUE,
                  
                  tieMethod = "random",
                  sampleSize = 1,

                  ## individuals can't be present in sample more than once (sample without replacement)
                  replace = FALSE,

                  ## probability that randomSplit will occur-- fed into setUpBasics as "rs" a global parameter
                  randomSplit = 0.2,

                  ## number of trees to grow
                  ## for normal sims set to 300
                  nTree = 300,

                  ## maximum number of covariates to consider-- a vector of covars for each decision point
                  mTry = rep(sqrt(10), n.stages),


                  ## Covariates for which the number of splits (s_i) is less
                  #'    than s*stratifiedSplit/d are explored preferentially
                  #     (s is the total number of splits, d is the
                  #'    total number of covariates under consideration).
                  stratifiedSplit = 0.1,
                  stageLabel = "_",
                  
                  #########
                  ######### NOTE: we need to change this nstrata if we want to test different things
                  nstrata = 2,
                  windowsize = 10)

  ## setting a different seed
  set.seed(sim*10000 + 5)

  ## feed arguments to dtrSurv() function: loops through each of the stages using Q-learning, using the predicted optimal from the previous stage and carrying it back
  ## the result should be the predicted optimal
  optimal.IHsurvrf <- do.call(IHsurvrf, c(arg.IHsurvrf2, list(nodeSize = nodesize, minEvent = mindeath )))

  ## if the first element of the class == try-error, this value will be TRUE
  IHsurvrf.error <- class(optimal.IHsurvrf)[1] == "try-error"

  ##if csk.error is FALSE (code works properly), then assign the results of the dtrSurv() which are optimal, to the "policy" slot of arg.csk
  ## this is of class DTRSurv
  arg.IHsurvrf$policy <- if (!IHsurvrf.error) optimal.IHsurvrf



  ## remove the results from dtrSurv()
  rm(optimal.IHsurvrf); gc()

  ############### evaluating the results of the estimated policy: we used the observed data with 100 patients (or whatever n is) to train the RF


  cat ("  2. IHsurvrf - Evaluation \n")
        ## setting a different seed
        set.seed(sim*10000 + 10)

        ## generate new data if IHsurvrf.error is FALSE (aka code runs properly)
        ## use arguments provided in original arg.csk to generate new multistage data-- with POLICY added
        ### this is the same set of arguments we used to generate multistage data for the observed data
        ### basically, we predict the patient's optimal action after feeding them through the trained RF
        ### the number of patients generated here is of n.eval which is 10000

        if (!IHsurvrf.error) IHsurvrf.data.rep <- do.call(simulate_patients, arg.IHsurvrf)


        ## val.fn depends on the criterion employed
        ## observed value for mean = take mean of the total event time (across all patients)
        ## in this case, the critical value will be null
        ## this is the observed policy from the data we generated
        ## if criterion is "surv.prob", "surv.mean", count the proportion of elements >= crit.value (ignoring NAs)
        ## crit.value is put at 5 (5 year survival rate) in C21.simulation_run.R

        ## calculate value function of the new generated data-- generated based on predicted optimal policy
        ## this should acount for ss as one of the inputs in arg.IHsurvrf which truncates the summary data to only look at 1:ss stages
        if (!IHsurvrf.error) result["IHsurvrf"] <- val.fn(IHsurvrf.data.rep$summary$cumulative.event.time)

        ## tt(2, reset = TRUE): resets timer for later measurements, and retreives the elapsed time in terms of mins for the time it takes to estimate and evaluate the policy
        result["time.IHsurvrf"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]

        ####### gets the number of iterations from the training data for convergence
        result["num.convergence.it"] <- unique(IHsurvrf.data.rep$summary$convergence.iterations)

        ## reset the policy to NULL and clean
        arg.IHsurvrf$policy <- NULL; gc()

        ## re,pve the data generated
        rm(IHsurvrf.data.rep); gc()
}

###############################################################
####################### ZOM SECTION ###########################
###############################################################


cat ("3. Estimation - Treatment 1 Only \n")
## all pts treated identically regardless of their characterstics
## chosen regime yields most favorable outcomes across pt population based on model predictions
## it can be seen as equivalent to the standard of care (aka every pt receives standard of care)
if (!skip.trt1) {

  ## create an object of claass "trt1" to input as the policy
  trt1 <- structure(1, class = "trt1")

  ## add the input policy which is "trt1"
  arg.trt1$policy <- trt1


  cat ("  3. Treatment 1 only (zero-order model) - Evaluation \n")
  set.seed(sim*10000 + 10)
  ## simulate multistage data using the estimated policy (trained tree) from the ZOM above
  trt1.data.rep <- do.call(simulate_patients, arg.trt1)

  result["trt1"] <- val.fn(trt1.data.rep$summary$cumulative.event.time)
  result["time.trt1"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]

  arg.trt1$policy <- NULL; gc()
  rm(trt1.data.rep); gc()
}

cat ("4. Estimation - Treatment 0 Only \n")
## all pts treated identically regardless of their characterstics
## chosen regime yields most favorable outcomes across pt population based on model predictions
## it can be seen as equivalent to the standard of care (aka every pt receives standard of care)
if (!skip.trt0) {

  ## create an object of claass "trt1" to input as the policy
  trt0 <- structure(0, class = "trt0")

  ## add the input policy which is "trt1"
  arg.trt0$policy <- trt0


  cat ("  3. Treatment 0 only (zero-order model) - Evaluation \n")
  set.seed(sim*10000 + 10)
  ## simulate multistage data using the estimated policy (trained tree) from the ZOM above
  trt0.data.rep <- do.call(simulate_patients, arg.trt0)

  result["trt0"] <- val.fn(trt0.data.rep$summary$cumulative.event.time)
  result["time.trt0"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]

  arg.trt0$policy <- NULL; gc()
  rm(trt0.data.rep); gc()
}



### saving and cleaning
print(result)

## mean of each column
#print(apply(result, 2, mean, na.rm = TRUE))
#saveRDS(result, filename.tmp) # saving the temporary results
gc()
# }

return(result)
}


# Number of cores to use for parallelization
num_cores <- 7

# Run the simulations in parallel
results_list <- mclapply(1:200, run_simulation, mc.cores = num_cores)

# Combine the results into a single dataframe
final_results <- do.call(rbind, results_list)

write.csv(final_results, "/nas/longleaf/home/js9gt/survrf/Outputs/testing_simplecode_300pt_hicens", row.names=FALSE)


## ------------------------------------- Function to look at the characteristics of our simulated data ------------------------------------ ##

simulate_and_analyze <- function(sim_values, arg.obs) {
  # Initialize placeholders for results
  censoring_percentages <- numeric(length(sim_values))
  final_stage_counts <- list()
  
  for (sim in sim_values) {
    # Set seed for reproducibility
    set.seed(sim * 10000)
    
    # Generate observed data
    obs.data <- do.call(simulate_patients, arg.obs)
    
    # Save the flowchart for reference
    flowchart_output <- flowchart(obs.data$output)
    #flowchart_filename <- paste0("flowchart_sim_", sim, ".png")
    #png(flowchart_filename)
    #print(flowchart_output)
    #dev.off()
    
    # Extract information from the last row
    censoring_percentages[sim] <- flowchart_output[nrow(flowchart_output), 4]
    
    # Extract information for final stage counts
    stages <- flowchart_output[1:n.stages, 1] # First column has stage indices (10 rows)
    censored <- flowchart_output[1:n.stages, 4] # Column 4 has number of censored patients
    died <- flowchart_output[1:n.stages, 5] # Column 5 has number of patients who died
    
    # Sum censored and died for each stage
    final_stage_counts[[sim]] <- censored + died
    
    
  }
  
  # Identify indices of `NULL` elements in final_stage_counts
  valid_indices <- !sapply(final_stage_counts, is.null)
  
  # Filter out `NULL` elements from both final_stage_counts and censoring_percentages
  final_stage_counts <- final_stage_counts[valid_indices]
  censoring_percentages <- censoring_percentages[valid_indices]
  
  
  
  # Combine results into a list for easy access
  results <- list(
    mean_censoring_percentage = mean(censoring_percentages, na.rm = TRUE),
    n_mean_total = length(censoring_percentages),
    final_stage_counts = final_stage_counts
  )
  
  return(results)
}


