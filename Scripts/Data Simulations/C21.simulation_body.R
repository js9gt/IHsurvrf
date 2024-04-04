
setwd("~/survrf/Scripts/IHsurvrf")
source("R/IH.dtrSurv.R")

setwd("~/survrf/Scripts/Data Simulations")
## generic is sourced in multistage_sim.R
source("F02.multistage_sim.R")

## checks value of the variable "criterion" to define value function
if (criterion == "mean") {
  val.fn <- function(x) {mean(x, na.rm = TRUE)}
} else if (criterion %in% c("surv.prob", "surv.mean")) {
  ## calculates the proportion of elements in x that are greater than or equal to crit.value
  val.fn <- function(x) {mean(x >= crit.value, na.rm = TRUE)}
}

## modifies a temporary filename by replacing ".rds" with "_tmp.rds"
#filename.tmp <- gsub("\\.rds", "_tmp.rds", filename)

## initialize matrix called stat.stage  with NA values with a row for each of the simulations run, and a column for each stage
stat.stage <- matrix(NA, nrow = n.sim, ncol = n.stages,
                     dimnames = list(1:n.sim, 1:n.stages))

########## simulation ##################
## initialize a data.frame called "result"  with columns for various statistics
## gets a result for each iteration of the simulation (AKA one row = values for one sim)
result <- data.frame(no = 1:n.sim, observed = NA, IHsurvrf = NA, zom = NA, time.obs = NA,
                     time.IHsurvrf = NA, time.zom = NA, percent.censor = NA)  # time for each method (both policy est and eval)

## assign an attribute to "result" containing a list of 2 elements
## stores value of criterion used
## stores value of critical value used
attr(result, "criterion") <- list(criterion = criterion, crit.value = crit.value)

## loop from 1 to n.stages to add columns to result, including columns named n_1, n_2, ...n.stages initialized with NA
for (i in 1:n.stages) result[[paste0("n_", i)]] = NA

## initialize a list called arg.obs with various parameters and settings used for analysis
## assigns the same value of the initialization list to all of these different variables:
## used as arguments for observed data, IHsurvrf, observed data with no censoring, and ZOM
## these are to be input into multistage data generation


## we want the no censor part, IHsurvRV, ZOM to have no censoring
arg.obs <- arg.IHsurvrf <- arg.obs.no.censor <- arg.zom <-
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
      censoringyesno = FALSE,

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
arg.obs.no.censor$n.sample <- arg.IHsurvrf$n.sample <- arg.zom$n.sample <-
  n.eval

## setting nodesize to 5
nodesize = 5
## setting minimum number of deaths by rounding sqrt of "nodesize" to nearest integer (no decimals)
mindeath = round(sqrt(c(nodesize)), 0)


# flow: obs (1) -> optimal (n.mc), obs.no.censor (n.mc)
print(Sys.time())
for (sim in 1:n.sim) {
  cat("###########################  simulation ", sim, "########################### \n")
  cat("########################### (criterion ", criterion, crit.value, ")################## \n")

  cat ("1. Data \n")


## tt(1): records start time
## in F00.generic.R
tt(1)

# generate multistage obs.data (allows for censoring) using arg.obs defined earlier
set.seed(sim*10000)
## do.call: calls a function and a list of arguments to be passed to it
obs.data <- do.call(simulate_patients, arg.obs)

# observed policy value
##### simulate a set of data again, when there is no censoring (AKA) using a censoring value of -10 for all surv.prev
## also use n.eval as number of samples
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
result[sim, "observed"] <- val.fn(obs.data.rep$summary$cumulative.event.time)

## note: 0 = censored, 1 = not censored, so we calculate the proportion of 1's from the generated observed data, then subtract from 1
####### USES OBSERVED DATA, not policy data
result[sim, "percent.censor"] <- 1 - mean(obs.data$summary$censor.status, na.rm = TRUE)

# select columns named "n_1", "n_2"... "n_nstages" which we created earlier
result[sim, paste0("n_", 1:n.stages)] <-
  ## for each stage, calculate the proportion of observations in each stage (s = whatever iteration of the n.stages we are in)
  ### AKA this is the proportion of patients who have their terminal stage at stage s
  sapply(1:n.stages, function(s) mean(obs.data$summary$terminal.stage == s))


result[sim, "time.obs"] <- tt(2, reset = TRUE)["elapsed"]

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

  # Generate the formulas for each stage-- NOTE that this does not include action (since we are not doing pooled analysis)
  models = lapply(1:n.stages, function(x) {
    # Create the formula with all terms for this stage
    paste0("Surv(T_", x, ", delta_", x, ") ~ baseline1_", x, " + baseline2_", x, " + state_", x, " + prior.visit.length_", x,
           " + cumulative.time_", x, " + nstages_", x, " + action.1.count_", x, " + action.0.count_", x,
           ## comment this line out if we don't want action
           " + A_", x
           ) %>%
      as.formula
  })

  ## create new list of arguments to use IHsurvrf
  ## used the observed data generated earlier (once it's in observable form)
  arg.IHsurvrf2 = list(data = data.df,
                  ## will give A_1, A_2, ... A_n.stages
                  txName = paste("A", 1:n.stages, sep = "_"),

                  ## list of models generated above for each stage
                  models = models,
                  ## use the observed lengths of the previous stages in the current stage's model as covariates
                  usePrevTime = TRUE,
                  ## maximum study length, set as 10 in C21.simulation_run.R
                  tau = tau,
                  ## distribution to draw timepoints from
                  timePoints = "uni",
                  ## number of timepoints from 0 to tau
                  nTimes = 200,

                  ## mean,surv.prob, or surv.mean
                  criticalValue = criterion,
                  ## the time at which to compare the survival probabilities, here, looking at 5 year survival (from the set value of crit.value for the criterion-- or NULL)
                  evalTime = crit.value,

                  ## if the criterion is mean, use truncated mean, otherwise, use logrank test
                  splitRule = ifelse(criterion == "mean", "mean", "logrank"),

                  ## use ERT
                  ERT = TRUE,

                  ## If 'ERT' and 'uniformSplit' are TRUE,
                  #'   the random cutoff is sampled from a uniform distribution over the range
                  #'   of available covariate values for tree splitting
                  uniformSplit = TRUE,

                  ## individuals can't be present in sample more than once (sample without replacement)
                  replace = FALSE,

                  ## probability that randomSplit will occur-- fed into setUpBasics as "rs" a global parameter
                  randomSplit = 0.2,

                  ## number of trees to grow
                  nTree = 300,

                  ## maximum number of covariates to consider-- a vector of covars for each decision point
                  mTry = rep(sqrt(p), n.stages),

                  ## false = using stratified analysis
                  pooled = TRUE,

                  ## Covariates for which the number of splits (s_i) is less
                  #'    than s*stratifiedSplit/d are explored preferentially
                  #     (s is the total number of splits, d is the
                  #'    total number of covariates under consideration).
                  stratifiedSplit = 0.1)

  ## setting a different seed
  set.seed(sim*10000 + 1)

  ## feed arguments to dtrSurv() function: loops through each of the stages using Q-learning, using the predicted optimal from the previous stage and carrying it back
  ## the result should be the predicted optimal
  optimal.IHsurvrf <- do.call(dtrSurv, c(arg.IHsurvrf2, list(nodeSize = nodesize, minEvent = mindeath )))

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
        ### same number of patients too
        if (!IHsurvrf.error) IHsurvrf.data.rep <- do.call(simulate_patients, arg.IHsurvrf)


        ## val.fn depends on the criterion employed
        ## observed value for mean = take mean of the total event time (across all patients)
        ## in this case, the critical value will be null
        ## this is the observed policy from the data we generated
        ## if criterion is "surv.prob", "surv.mean", count the proportion of elements >= crit.value (ignoring NAs)
        ## crit.value is put at 5 (5 year survival rate) in C21.simulation_run.R

        ## calculate value function of the new generated data-- generated based on predicted optimal policy
        if (!IHsurvrf.error) result[sim, "IHsurvrf"] <- val.fn(IHsurvrf.data.rep$summary$cumulative.event.time)

        ## tt(2, reset = TRUE): resets timer for later measurements, and retreives the elapsed time in terms of mins for the time it takes to estimate and evaluate the policy
        result[sim, "time.IHsurvrf"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]

        ## reset the policy to NULL and clean
        arg.IHsurvrf$policy <- NULL; gc()

        ## re,pve the data generated
        rm(IHsurvrf.data.rep); gc()
}

###############################################################
####################### ZOM SECTION ###########################
###############################################################


cat ("6. Estimation - zero-order model\n")
## all pts treated identically regardless of their characterstics
## chosen regime yields most favorable outcomes across pt population based on model predictions
## it can be seen as equivalent to the standard of care (aka every pt receives standard of care)
if (!skip.zom) {
  cat ("  6. zero-order model - Policy estimation \n")
  set.seed(sim*10000 + 5)
    ## If an error or warning occurs during the execution, it will be captured, and the execution won't halt
    ## Instead, the program will continue to execute
  optimal.zom <- try(dtrSurv(data = data.df,
                                       txName = paste("A", 1:n.stages, sep = "_"),
                                       models = models,
                                       usePrevTime = TRUE,
                                       tau = tau,
                                       timePoints = "uni",
                                       nTimes = 200,
                                       criticalValue = criterion, evalTime = crit.value,
                                       splitRule = ifelse(criterion == "mean", "mean", "logrank"),
                                       ERT = TRUE, uniformSplit = TRUE, replace = FALSE,
                                       randomSplit = 0.2,

                                       ##### these are different -- these were:
                                       ## minEvent: round(sqrt(c(nodesize)), 0)
                                       # this means that a node can be split even if there's only 1 event in a node
                                       ## nodeSize = 5
                                       # large node size: 10,000 to indicate that all pts remain in root node (unless > 10,000)
                                       minEvent = 1, nodeSize = 1e+5, # zero-order model
                                       #########

                                       nTree = 300,

                                       #### this is different--
                                       ## was rep(sqrt(p), n.stages)
                                       # we only consider 1 covariate for splitting
                                       mTry = rep(1, n.stages),
                                       pooled = FALSE,

                                       #### this is different--
                                       ## was 0.1
                                       # means there's no preferential splitting for using covariates that are underused
                                       stratifiedSplit = 0))

  zom.error <- class(optimal.zom)[1] == "try-error"
  ## add the estimated policy from ZOM to the argumentto simulate multistage data
  arg.zom$policy <- if (!zom.error) optimal.zom
  rm(optimal.zom); gc()


  cat ("  6. zero-order model - Evaluation \n")
  set.seed(sim*10000 + 10)
  ## simulate multistage data using the estimated policy (trained tree) from the ZOM above
  if (!zom.error) zom.data.rep <- do.call(simulate_patients, arg.zom)
  if (!zom.error) result[sim, "zom"] <- val.fn(zom.data.rep$summary$cumulative.event.time)
  result[sim, "time.zom"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
  arg.zom$policy <- NULL; gc()
  rm(zom.data.rep); gc()
}


### saving and cleaning
print(result[sim, ])
print(apply(result, 2, mean, na.rm = TRUE))
#saveRDS(result, filename.tmp) # saving the temporary results
gc()
# }
}
result


write.csv(result, "/nas/longleaf/home/js9gt/survrf/Outputs/pt75.csv", row.names=FALSE)




