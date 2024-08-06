
#######################################################
### multi-stage data preprocessing to append stage names
#######################################################

## data to test on generated in F02.multistage_sim.R
## a <- output2observable(pts)

## Convert multi-stage dataset, "output" into observable format for analysis

output2observable <- function(output, stage = NULL) {

  # gets number of stages from the 3rd dimension of the "output" array

  n.stages = dim(output)[3]

  # if no stage provided, considers all stages
  if (is.null(stage)) {
    stage = 1:n.stages

  }

  # number of rows in output dataset (n)
  n = dim(output)[1]

  # defining variable names related to covariates & stage variables
  nm = dimnames(output)[[2]]

  ## names of the covariates
  nm.covar = c("baseline1", "baseline2", "state1", "state2", "prior.visit.length", "cumulative.time", "nstages", "action.1.count", "action.0.count")
  nm.stage =  c("event.time", "delta", "gamma", "action", nm.covar)

  # initializes empty DF with subject ID column for the first stage, taken from "output"
  df <- data.frame(subj.id = output[, "subj.id", 1])

  # iterates through each stage to make
  # df.i: DF containing columns related to Time, delta, action, covariates for stage i
  # mergecolumns of this into main dataframe (DF)
  for (i in stage) {
    df.i <- data.frame(matrix(output[, nm.stage, i], ncol = length(nm.stage), byrow = FALSE))
    names(df.i) <- paste(c("T", "delta", "gamma", "A", nm.covar), i, sep = "_")
    df <- cbind(df, df.i)
  }

  # calculates new delta value in DF
#  df$delta =
#    # selects columns starting with "delta_"
#    df %>% dplyr::select(starts_with("delta_")) %>%
#    # converts selected cols into matrix
#    as.matrix %>%
#    # applies function to each row to see if any values == 0 (Censoring)
#    # this returns a 1 if there's no 0, returns a 0 if there is a 0
#    apply(1, function(s) 1 - any(s == 0, na.rm = TRUE))
  df
}


# time calculator
tt <- function(s, reset = FALSE, units = "auto"){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    result <- data.frame(begin = time.tmp, end = Sys.time(), elapsed = difftime(Sys.time(), time.tmp, units = units))
    if (reset) time.tmp <<- Sys.time()
    return(result)
  }
}


flowchart <- function(outputData) {
  n.stages = dim(outputData)[3]
  n.sample = dim(outputData)[1]
  result <- data.frame(stage = c(1:n.stages, "Total", "Percent"),
                       total = rep(NA, n.stages + 2), percentage = NA,
                       censored = NA, died = NA, nextStage = NA, mean.trtdiff = NA)
  for (stage in 1:n.stages) {
    ## mean treatment difference is the mean of everyone in that stage
    result[stage, "mean.trtdiff"] = mean(outputData[,,stage][,"trt.diff"], na.rm = T)
    
    result[stage, "total"]  = dim(outputData)[1]
    result[stage, "percentage"]  = round(result[stage, "total"] / n.sample,2)
    cens <- outputData[, "delta", stage] == 0
    # drop those censored
    outputData <- outputData[!cens,,, drop = FALSE]
    died <- outputData[, "gamma", stage] == 1
    # drop those died
    outputData <- outputData[!died,,, drop = FALSE]
    result[stage, "censored"]  = sum(cens, na.rm = TRUE)
    result[stage, "died"]      = sum(died, na.rm = TRUE)
    result[stage, "nextStage"] = sum(!died, na.rm = TRUE)
  
    
  }
  result[n.stages + 1, c("censored", "died")] <- apply(result[, c("censored", "died")], 2, sum, na.rm = TRUE)
  result[n.stages + 2, c("censored", "died")] <- round(result[n.stages + 1, c("censored", "died")]/n.sample, 2)
  # result[c("Total", "Percent"), 1:2]
  result
}
