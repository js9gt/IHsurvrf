

# Internal function to predict value of a tree
#
# Function is not exported
#
# @param x A data.frame object {n x p}. The model.frame for covariates to be
#   considered in splitting
#
# @param params A Parameters object. All information that regulates tree and
#   specifies analysis preferences.
#
# @param nCat An integer object. Number of categories in each covariate. If
#   covariate i is continuous, nCat[i] = 0; if covariate i is an ordered factor,
#   nCat[i] = 1; if covariate i is an unordered factor, nCat[i] = levels()
#
# @param nodes A list object. The nodes of a tree.
#

#' @include class_IH.Parameters.R

## defines a new function .predictSurvTree used for predicting outcomes from a survival tree
.predictSurvTree <- function(...,

                             ## A data.frame object {n x p}. The model.frame for covariates to be
                             #   considered in splitting (new data)
                             x,

                             ## A Parameters object. All information that regulates tree and
                             #   specifies analysis preferences.
                             params,

                             ## An integer object. Number of categories in each covariate. If
                             #   covariate i is continuous, nCat[i] = 0; if covariate i is an ordered factor,
                             #   nCat[i] = 1; if covariate i is an unordered factor, nCat[i] = levels()
                             nCat,

                             ## A list object. The nodes of a tree.
                             nodes) {

  ## returns NULL if the input data is missing (function can only proceed w/ valid data)

  if (is.null(x = x)) return( NULL )

  ## converts the input data frame into a matrix format for processing

  x <- data.matrix(frame = x)

  ## idenfities which nodes are used (not NA) and counts them.
  ## extracts the first column of the matrix for non-na values

  usedNodes <- !is.na(x = nodes$nodes[,1L])

  ## calculates the number of nodes used (not NA values)
  nUsed <- sum(usedNodes)

  ## number of observations in the data frame

  n <- nrow(x = x)

  ## number of predictors in the data frame
  np <- ncol(x = x)

  ## retrieves number of time points from the "params" object
  nTimes <- .NTimes(params)

  ## set NA values in the "nodes" information to 0

  nodes$nodes[is.na(nodes$nodes)] <- 0.0

  ## invoke a Fortran subroutine called "predictSurvTree" to compare survival predictors for each obs based on the tree structure & covars
  ## pass some parameters to fortran\

  ## dtrSurv.f90 fortran code

  res <- .Fortran("predictSurvTree",

                  ## number of observations in df
                  n = as.integer(n),

                  ## number of predictors in the df
                  np = as.integer(np),

                  ## numeric matrix representation of x for fortran processing
                  xt = as.double(x = x),

                  ## vector indicating type of predictor (continuous, ordered, unordered)
                  nCat = as.integer(x = nCat),

                  ## total number of survival functions
                  nt = as.integer(x = nrow(x = nodes$survFunc)),

                  ## number of nodes used in tree
                  nNodes = as.integer(x = nUsed),

                  ## extracted survival funcs,mean survival times,survival probs in the "usedNodes" (logical)
                  ## meaning, for a single tree, extract the survival functions, means, survival probability
                  tsurvFunc = as.double(x = nodes$survFunc[,usedNodes]),
                  mean = as.double(nodes$mean[usedNodes]),
                  survProb = as.double(nodes$survProb[usedNodes]),
                  nCols = as.integer(x = ncol(x = nodes$nodes)),
                  tnodes = as.double(x = nodes$nodes[usedNodes,]),

                  ## ore-allocated vectors to be filled by fortran
                  predSurvFunc = as.double(numeric(length = n*nTimes)),
                  predMean = as.double(numeric(length = n)),
                  predSurvProb = as.double(numeric(length = n)),

                  ## specifies the R package associated with the fortran code
                  PACKAGE = "IHsurvrf")

  ## checks if the "params" critical value is a mean or a probability that matches either of surv.mean or surv.prob
  ## ouputs a logical vector

  isSurv <- .CriticalValueCriterion(params) %in% c("surv.mean", "surv.prob")

  ## construct a list to store matrices and results
  valueObj <- list()

  ## create a column called "survFunc" with the numeric vector of "res$predSurvFunc" converted into a matrix
  ## matrix where computed probabilities at each time point for each individual is contained
  ## each column is a single individual's survvial probabilities across all "nTimes"

  valueObj[[ "survFunc" ]] <- matrix(data = res$predSurvFunc, nrow = nTimes, ncol = n)

  ## storing the res$predMean
  valueObj[[ "mean" ]] <- res$predMean

  ## checks if analysis involves survival probabilities using a logical defined earlier
  ## if isSurv is TRUE, store the predicted survival probabilities for each observation
  if (isSurv) valueObj[[ "survProb" ]] <- res$predSurvProb

  return( valueObj )

}
