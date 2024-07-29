

#### C21.simulation_body.R is to be run by source("C21.simulation_run.R")


# n             # training sample size
# n.eval        # evaluation sample size
# n.sim         # number of simulation replicates
# p = 5         # number of covariates (should agree with beta coefficients' dimension)
# tau = 1000    # total study length
# n.stages = 50 # number of stages

### 0. Get the setting number
## This section retrieves the settings from the command line arguments.
## It expects four arguments: 1 = beta, 2 = propensity, 3 = size, and 4 = crit.
## these are used to set the different beta, propensity, size, crit settings
## If these arguments are not provided, it sets default values.
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) < 4) {
  warning("commandArgs was not provided. Being set as c(1,1,1,1).")
  message(" 300 patients, max stages = 10, tau = 2,000, 200 simulation replicates, 10,000 eval, 
          40% censoring, obs setting (not RCT),
          with full convergence no forest 2 appending")
   arg = c(1, 1, 1, 1) # by default
  print(arg)
}

## renames the arguments for clarity
## also converts arguments to numeric format
names(arg)[1:4] = c("beta", "propensity", "size", "crit")
arg1 <- as.numeric(arg[1]) # 1..4
arg2 <- as.numeric(arg[2]) # 1..2
arg3 <- as.numeric(arg[3]) # 1..2
arg4 <- as.numeric(arg[4]) # 1..3
## extracts a date if provided, otherwise, set as the current date
arg.date <- if (is.na(arg[5]) | arg[5] == "") Sys.Date() else as.character(arg[5])


### 1. setup
# n.stages = 3 # number of stages
default <- list(
  ## evaluation sample size: 10000
  n.eval = 10000,
  ## number of simulation replicates: 200
  n.sim = 200,
  ## tau (days): total study length
  tau = 2000,
  ## maximum number of stages
  n.stages = 10,
  ## the stage we start at since we don't want issues with too small sample size
  ss = NULL)

# arg4 crit: a list containing different 2 different criterion with associated values
crit <- list(crit1 = list(criterion = "mean", crit.value = NULL, value = "truncated mean E[T]"),
             crit2 = list(criterion = "surv.mean",
                          ## 2 year survival in terms of days
                          crit.value = 730, value = "S(730)"))

# arg3 size: a list containing 2 training sample sizes
size <- list(
  ## small ss is 75
    small.sample.size = list(n = 300),
  ## large ss is 10,000
             large.sample.size = list(n = 10000))

# arg2 propensity
propensity <-   # (int), state
  list(obs = list(beta.propensity = function(p) c(0, -rep(0.5, p))),  # no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(0, rep(0, p))))    # RCT

### NOTE: *I think* that we would change different settings for alpha (similar to arg1 beta) to get different censoring rates
## maybe also dif # of baseline variables, and different number of random state variables?
# arg1 beta: a list containing beta values for 4 different settings
# Create a list for the settings
coefs <- list(
  setting1 = list(
    coef_failure = list(

      ### 10% censoring

      ##  a1 + b1 * state + c1*state2 + z1 * nstages + p1 * cumulative.length + g1 * action + h1*prior.visit.length +
      ###    r1 * action * state * nstages * cumulative.length*prior.visit.length

      a = -6, b = -0.2, c = -0.5, z = -0.025, p = -0.02, g = 0.1, h = -0.08, r = 0.05
    ),
    coef_nextvisit = list(
      a = -1.5, b = -0.2, c = -0.5, z = -0.025, p = -0.02, g = 0.1, h = -0.08, r = 0.05
    ),
    coef_censoring = list(
      ## a = -6 for 40% censoring with tau = 2000
      a = -17, b = -0.2, c = -0.5, z = -0.025, p = -0.02, g = 0.1, h = -0.08, r = 0.05
    )
  ),

  ## higher censoring
  setting2 = list(
    coef_failure = list(

      ##  a1 + b1 * state + z1 * nstages + p1 * cumulative.length + g1 * action + h1*prior.visit.length +
      ###    r1 * action * state * nstages * cumulative.length*prior.visit.length

      a = -7, b = -0.3, c = -0.5, z = -0.025, p = 0.02, g = -0.2, h = 0.008, r = 0.01
    ),
    coef_nextvisit = list(
      a = -5, b = -0.3, c = -0.05, z = -0.015, p = 0.025, g = -0.2, h = 0.008, r = 0.01
    ),
    coef_censoring = list(
      a = -9.5, b = -0.3, c = -0.6, z = -0.04, p = 0.02, g = -0.2, h = 0.008, r = 0.01
    )
  )
)

# Set names
names(coefs$setting1$coef_failure) <- paste0(names(coefs$setting1$coef_failure), "1")
names(coefs$setting1$coef_nextvisit) <- paste0(names(coefs$setting1$coef_nextvisit), "2")
names(coefs$setting1$coef_censoring) <- paste0(names(coefs$setting1$coef_censoring), "3")
names(coefs$setting2$coef_failure) <- paste0(names(coefs$setting2$coef_failure), "1")
names(coefs$setting2$coef_nextvisit) <- paste0(names(coefs$setting2$coef_nextvisit), "2")
names(coefs$setting2$coef_censoring) <- paste0(names(coefs$setting2$coef_censoring), "3")

# p.list: List of the number of covariates for the coefficient settings
## takes off 2 from intercept, nstages, cumulative.length, action, prior.visit.length, interaction between all of them
p.list <- (length(coefs[[arg1]]$coef_failure) - 6)

# setting1 n.boot 50 / n 300
setting = c(arg = list(arg), default, coefs[[arg1]]$coef_failure, coefs[[arg1]]$coef_nextvisit,
            coefs[[arg1]]$coef_censoring, p = p.list,
            propensity[[arg2]], size[[arg3]], crit[[arg4]])

### 2. put a selected setting into the global environment
cat("setting (beta, propensity, size, crit) ", arg, "\n")
list2env(setting, envir = globalenv())
beta.propensity <- beta.propensity(p)

## creates a "simResult" document with the different arguments, this is inside an "output" folder
#filename = paste0("output/simResult_", arg.date, "_beta", arg1, "_prop", arg2,
#                  "_n", arg3, "_crit", arg4, ".rds")

predictedPropensity <- function(state1, state2) {
  # Compute the predicted probability using the logistic function
  plogis(cbind(1, state1,state2) %*% beta.propensity)
}

### we don't need the predicted hazard, censoring, etc.

####### initializes 4 of the same empty lists with a length equal to n.stages
predictedPropensityvec <- vector("list", n.stages)

## iterates over each stage to basically write the same function in each slot of the vector
## basically use this so each stage uses the same functions
for (i in 1:n.stages) {
  predictedPropensityvec[[i]] <- predictedPropensity
}

### 3. Run the simulation
skip.IHsurvrf <- FALSE
skip.trt1 <- FALSE
skip.trt0 <- FALSE
skip.opt <- TRUE
### commented out since we don't need it
## skip.gk <- skip.dw <- TRUE
cv.nodesize = FALSE

source("~/survrf/Scripts/Data Simulations/C21.simulation_body.R")

