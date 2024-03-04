


dtrSurv <- function(data,
                    txName,
                    models, 
                    ...,
                    usePrevTime = TRUE,
                    timePoints = "quad",
                    nTimes = 100L,
                    tau = NULL,
                    criticalValue = "mean",
                    evalTime = NULL,
                    splitRule = NULL,
                    ERT = TRUE,
                    uniformSplit = NULL,
                    sampleSize = NULL,
                    replace = NULL,
                    
                    ## fed into fortran via setUpBasics as "rs" (global param)
                    randomSplit = 0.2,
                    tieMethod = "random",
                    minEvent = 3L,
                    nodeSize = 6L, 
                    nTree = 10L,
                    mTry = NULL,
                    pooled = FALSE,
                    stratifiedSplit = NULL,
                    stageLabel = ".") {
  
  
  
  # total number of individuals in dataset
  nSamples <- nrow(x = data)


params <- .parameters(timePoints = timePoints,
                      tau = tau,
                      nTimes = nTimes,
                      
                      ## response is created (not input)
                      response = response,
                      nTree = nTree, 
                      ERT = ERT, 
                      uniformSplit = uniformSplit, 
                      randomSplit = randomSplit, 
                      splitRule = splitRule,
                      replace = replace, 
                      nodeSize = nodeSize, 
                      minEvent = minEvent, 
                      tieMethod = tieMethod,
                      criticalValue = criticalValue, 
                      survivalTime = evalTime,
                      
                      ## nSamples is created (not input)
                      nSamples = nSamples,
                      pooled = pooled,
                      stratifiedSplit = stratifiedSplit)

} 
