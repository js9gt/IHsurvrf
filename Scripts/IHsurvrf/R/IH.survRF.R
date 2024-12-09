



source("R/class_IH.Parameters.R")
#source("~/survrf/Scripts/IHsurvrf/R/IH.VerifyERT.R")
source("R/class_IH.CriticalValue.R")
#source("~/survrf/Scripts/IHsurvrf/R/IH.VerifyERT.R")
source("R/class_IH.SurvRF.R")
#source("~/survrf/Scripts/IHsurvrf/R/IH.VerifyERT.R")



## function generates a survival random forest-- collection of decision trees used to estimate survival function for individuals based on their covariates

.survRF <- function(...,

                      ## A data.frame object {n x p}. The model.frame for covariates to be
                      #   considered in splitting
                      x,

                      ## An integer vector object {n}. Vector of indicators for censoring
                      #   (1 = not censored; 0 = censored)

                      delta,

                      ## A matrix object {nt x n}. Probability mass vector of survival
                      #   function
                      pr,

                      ## An integer object. The number of trees to grow.
                      params,


                      ## An integer object. The maximum number of covariates to use
                      #   for splitting algorithm.
                      mTry,

                      ## An integer object. The number of samples to draw for each
                      #    tree

                      sampleSize
                    ) {
  ## use fortran functions "setUpInners" and "survTree"
  ## prepares and initializes data structures and variables for later steps where we actually grow the forest


  # if x_i is an unordered factor, nCat_i is the number of levels
  # if x_i is not a factor, nCat is 0
  # if x_i is an ordered factor, nCat is 1

  ## iterates over each covariate in the dataset & extracts its levels
  ## returns NULL for numerical variables
  ## store the actual levels for each factor

  xLevels <- lapply(X = x, FUN = levels)

  ## applies the "nlevels" function across all covariates in x-- counting the number of levels
  ## returns 0 for non-factor covariates
  ## stores the number of levels each category has

  nCat <- sapply(X = x, FUN = nlevels)

  ## further checks if a covariate is an ordered factor. If this is true, set to 1 for that covar, otherwise retain original value

  nCat <-
    ifelse(test = sapply(X = x, FUN = is.ordered),
           yes = 1L,
           no = nCat)

  ## calculates total number of samples (individuals) in the dataset

  # number of individuals in training data
  nSamples <- nrow(x = x)

  ## calculates total number of unique time points from the probability mass vector (designed for .shiftMatrix output)

  # number of time points: this doesn't change from the replicate pr so it doesn't matter which one we use
  nTimes <- nrow(x = pr)

  # total number of trees to be grown in the forest
  # .NTree() is a getter method defined for Parameters objects, returning the "nTree" slot
  ## from class_TreeConditions.R

  nTree <- .NTree(object = params)

  ## adjusts the "sampleSize" parameter by multiplying by the total numerb of samples & rouding up to integer
  ## sampleSize is an input value * number of individuals in training data
  ## number of samples drawn for building each tree

  # determine the number of samples to include in each tree
  ## need nSamples as input
  sampleSize <- ceiling(x = sampleSize*nSamples)

  ## Print sampleSize

  message("sampleSize: ", sampleSize)

  ## maximum number of nodes based on sample size to ensure there's enough room for binary tree structure
  ## add one tree to account for root node

  # maximum number of nodes in a tree
  maxNodes <- 2L * sampleSize + 1L

  ## converts data frame "x" into a matrix format & convert factor variables into integers for Fortran processing

  # convert factors to integers
  x = data.matrix(frame = x)

  ## after conversion to a matrix, check again (just like nSamples) the number of observations

  nr = nrow(x = x)

  ## x: data frame w/ for covars to be considered in splitting
  # pr: matrix of probability mass vector of survival function
  # delta: 1 = not censored, 0 = censored
  # mTry: integer; max number of covars to use for splitting
  # sampleSize: initeger; number of samples to draw for each tree
  # from params: input vals into dtrSurv --> input into .dtrSurvStep --> input into .survRF
  ##### nCat: calculated through nlevels () in base R
  # nTree: input into dtrSurv
  ##### maxNodes: calculated based on input sample size





  res = .Fortran(
    "setUpInners",
    # number of cases under consideration
    t_n = as.integer(x = nrow(x = x)),
    # number of covariates
    t_np = as.integer(x = ncol(x = x)),
    # the covariates
    t_x = as.double(x = x),
    # probability mass vector of the survival functions
    t_pr = as.double(x = t(x = pr)),

    # indicator of censoring
    t_delta = as.integer(x = delta),
    # maximum number of covariates to try for splitting
    t_mTry = as.integer(x = mTry),
    # number of categories in each covariate
    t_nCat = as.integer(x = nCat),
    # number of cases to sample for each tree
    t_sampleSize = as.integer(x = sampleSize),
    # number of trees in the forest
    t_nTree = as.integer(x = params@nTree),
    # maximum number of nodes
    ## need this as input
    t_nrNodes = as.integer(x = maxNodes),
    PACKAGE = "IHsurvrf"
  )


  ## grows the survival forest over "nTrees" which is a global variable which will need to be declared
  ## the inputs just set the dimensions
  ## also use input "replace" from dtrSurv.R as a global variable
  ## also use input "stratifiedSplit" from dtrSurv.R as a global variable

  survTree <- .Fortran(
    "survTree",
    ## output: survival function averaged over forest
    ## number of time points (gotten earlier) * number of samples
    forestSurvFunc = as.double(numeric(nTimes*nr)),
    ## output: mean survival averaged over forest
    ## number of samples
    forestMean = as.double(numeric(nr)),
    ## output: survival probability averaged over forest
    ## number of samples
    forestSurvProb = as.double(numeric(nr)),
    PACKAGE = "IHsurvrf"
  )

  # retrieve trees from Fortran

  ## create an empty list to store the details of each tree generated by Fortran
  trees <- list()

  ## iterate over specified number of trees to fetch its dimensions and node details in fortran
  for (i in 1L:nTree) {
    ## initializing a new list for each tree

    trees[[i]] <- list()


    ## retrieves tree dimensions
    treeDims <- .Fortran(
      "treeDim",
      ## input the i-th tree
      iTree = as.integer(x = i),

      ## output: number of rows to retrieve in node matrix.
      ## 1L input is just a placeholder to provide a starting value before its updated by the actual output
      nr = as.integer(x = 1L),

      ## output: number of columns to retrieve in node matrix
      ## 1L input is just a placeholder to provide a starting value before its updated by the actual output
      nc = as.integer(x = 1L),
      PACKAGE = "IHsurvrf"
    )


    ## retrieves the details of the i-th tree


    temp <- .Fortran(
      "getTree",
      ## input the i-th tree
      iTree = as.integer(x = i),
      ## number of rows from the treeDims object
      nr = as.integer(x = treeDims$nr),
      ## number of columns from the treeDims object
      nc = as.integer(x = treeDims$nc),
      ## number of nodes = number of rows* number of columns
      nodes = as.double(x = numeric(length = treeDims$nc*treeDims$nr)),
      ## number of rows * timepoints
      survFunc = as.double(x = numeric(length = treeDims$nr*nTimes)),
      ## number of rows (each observation)
      mean = as.double(x = numeric(length = treeDims$nr)),
      ## number of rows
      survProb = as.double(x = numeric(length = treeDims$nr)),
      PACKAGE = "IHsurvrf"
    )

    ## store the retrieved details & create separate compoenents for nodes, survFunc, mean in the "trees" list

    ## convert the "node" vector from "temp" output into a matrix representing the trees structure
    ## matrix has a row for each node and columns corresponding to node attributes
    ## the vector itself already has node information of the entire tree-- we just reshape this vector into a matrix
    trees[[i]]$nodes <- matrix(data = temp$nodes,
                               nrow = treeDims$nr,
                               ncol = treeDims$nc)

    ## turn survival function data into matrix with one row per time point, one column per node
    trees[[i]]$survFunc <- matrix(data = temp$survFunc,
                                  nrow = nTimes,
                                  ncol = treeDims$nr)
    trees[[i]]$mean <- temp$mean

    trees[[i]]$survProb <- temp$survProb
  }

  ## initializes empty list to store aggregated results from the forest

  forest <- list()

  ## create a new aspect of the list called "survFunc" that holds survival function data from ALL TREES into a matrix
  ## each row of matrix corresponds to a time point
  ## each column corresponds to an individual in dataset

  forest[["survFunc"]] <- matrix(data = survTree$forestSurvFunc,
                                 nrow = nTimes,
                                 ncol = nr)

  ## create a new aspect called "mean" that stores the mean times calculated across all trees for individuals in the dataset
  ## presents summary statistic representing expected survvial time for each individual as predicted by the forest
  ## mean survival time across all trees within the forest to provide a single mean survival time prediction for each individual

  forest[["mean"]] <- survTree$forestMean

  ## if the critical value matches one of mean of probability, include the survival probabilities of the overall forest too
  ## RECALL: critical value either mean survival time, or survival probability at that stage
  ## we want to pick the action that maximizes the mean survival time or the survival probability
  ## from class_IH.CriticalValueSurvival.R if it's a survival based object
  ## from class_IH.CriticalValueMean.R if it's a mean based object
  ## returns an error if the class of the object is ANY

  crit <- .CriticalValueCriterion(params)


  ## if this type is "surv.mean" or "surv.prob" we retrieve the survival probability
  #     - surv.prob = mean survival probability at time “evalTime”
  # surv.mean = hybrid; first use mean survival probability, then if there are ties use the mean survival time to identify optimality
  if (crit %in% c("surv.mean", "surv.prob")) {
    forest[["survProb"]] <- survTree$forestSurvProb
  }

  ## outputs object of class SurvRF

  return(
    new(
      Class = "SurvRF",

      ## list of individual trees that comprise the forest containing detailed node and survival function data

      "trees" = trees,

      ## the aggregated results from the overall forest

      "forest" = forest,
      "variables" = colnames(x = x),
      "mTry" = mTry,
      "nCat" = nCat,

      ## levels for each categorical variables
      "xLevels" = xLevels
    )
  )
}
