
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
  nm.covar = c("baseline1", "baseline2", "state", "prior.visit.length", "cumulative.time", "nstages", "action_1_count", "action_0_count")
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

  # calculates new delta value in DF
  df$delta =
    # selects columns starting with "delta_"
    df %>% dplyr::select(starts_with("delta_")) %>%
    # converts selected cols into matrix
    as.matrix %>%
    # applies function to each row to see if any values == 0 (Censoring)
    # this returns a 1 if there's no 0, returns a 0 if there is a 0
    apply(1, function(s) 1 - any(s == 0, na.rm = TRUE))
  df
}
