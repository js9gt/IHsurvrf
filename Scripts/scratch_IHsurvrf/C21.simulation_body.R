

#### C21.simulation_body.R is to be run by source("C21.simulation_run.R")

#### library
# install.packages("survival")
# install.packages("randomForestSRC")
# install.packages("DTRreg")
# install.packages("dtrSurv")
setwd("~/survrf/Scripts/scratch_IHsurvrf")
source("R/IH.dtrSurv.R")

setwd("~/survrf/Scripts/scratch_IHsurvrf")
source("F00.generic.R")

setwd("~/survrf/Scripts/scratch_IHsurvrf")
source("F02.multiStage.R")

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
    result <- data.frame(no = 1:n.sim, observed = NA, csk = NA, gkLM = NA, gkRF = NA, 
                         dw = NA, zom = NA,
                         time.obs = NA, time.csk = NA, time.gkLM = NA, time.gkRF = NA,
                         time.dw = NA, time.zom = NA, percent.censor = NA)  # time for each method (both policy est and eval)
    
    ## assign an attribute to "result" containing a list of 2 elements
      ## stores value of criterion used
      ## stores value of critical value used
    attr(result, "criterion") <- list(criterion = criterion, crit.value = crit.value)
    
    ## loop from 1 to n.stages to add columns to result, including columsn named n.1, n.2, ...n.stages initialized with NA
    for (i in 1:n.stages) result[[paste0("n.", i)]] = NA
    
    ## initialize a list called arg.obs with various parameters and settings used for analysis
    ## assigns the same value of the initialization list to all of these different variables:
    ## arg.obs = arg.csk = arg.obs.no.censor = arg.gk.lm = arg.gk.rf = arg.dw = arg.zom
    arg.obs <- arg.csk <- arg.obs.no.censor <- arg.gk.lm <- arg.gk.rf <-  arg.dw <-  arg.zom <- 
      list(
        n.sample = n, n.stages = n.stages, tau = tau, tick = 0.01,  # structural parameters
        at.risk = 1,      # initial state vector
        p = p, corr = -0.5,                              # cor of two error processes
        
        ### the following are lists with n.stages elements, each element is the same function
        ## ifelse(action == 1, cbind(1, log(surv.previous + 1), covariate) %*% beta.hazard1, cbind(1, log(surv.previous + 1), covariate) %*% beta.hazard0)
        predHazardFn = predHazardFn, 
        
        ### plogis(cbind(1, log(surv.previous + 1), covariate) %*% beta.propensity)
        predPropensityFn = predPropensityFn, 
        
        
        ## ifelse(action == 1, cbind(1, log(surv.previous + 1), covariate) %*% beta.censor1, cbind(1, log(surv.previous + 1), covariate) %*% beta.censor0)
        predCensorFn = predCensorFn,             # list of predictor functions
        hidden_data = TRUE, printFlag = FALSE
      )
    
    ## we set the predCensorFn of all these objects to the value of "noCensorFn"
    ### creates matrix with 1 column and rows = number of surv.previous values
    ### matrix(-10, nrow = length(surv.previous), ncol = 1)
    ## surv.previous is the cumulative time the patient survived-- 
    ## since tau is 10, this is saying all patients survived the whole study 
    arg.obs.no.censor$predCensorFn <- arg.csk$predCensorFn <- arg.gk.lm$predCensorFn <- 
      arg.gk.rf$predCensorFn <- arg.dw$predCensorFn <- arg.zom$predCensorFn <-
      noCensorFn
    
    ## setting n.sample attribute of ll these objects (created above) to "n.eval"
    #### n.eval is a default setting read into C21.simulation_run
    #### this is the number of samples per simulation
    arg.obs.no.censor$n.sample <- arg.csk$n.sample <- arg.gk.lm$n.sample <- 
      arg.gk.rf$n.sample <- arg.dw$n.sample <- arg.zom$n.sample <-
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
        # obs.data using 
        set.seed(sim*10000)
        obs.data <- do.call(multiStageDynamics, arg.obs)
        
        # observed policy value
        ##### simulate a set of data again, when there is no censoring (AKA) using a censoring value of -10 for all surv.prev 
        ## also use n.eval as number of samples 
        set.seed(sim*10000)
        obs.data.rep <- do.call(multiStageDynamics, arg.obs.no.censor)
        
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
        
        # select columns named "n.1", "n.2"... "n.nstages" which we created earlier
        result[sim, paste0("n.", 1:n.stages)] <- 
          ## for each stage, calculate the proportion of observations in each stage (s = whatever iteration of the n.stages we are in)
          ### AKA this is the proportion of patients who have their terminal stage at stage s
          sapply(1:n.stages, function(s) mean(obs.data$summary$terminal.stage == s))
        
        ## tt(2, reset = TRUE): resets timer for later measurements, and retreives the elapsed time in terms of seconds
        
        result[sim, "time.obs"] <- tt(2, reset = TRUE)["elapsed"]
        print(flowchart(obs.data$output))
        
        # transforming data from an array format to a data.frame format
        data.df = output2observable(obs.data$output)
        
        ## removes the data generated for the policy (with n.eval iterations) from the R environment to free up memory space,also runs garbage collection to free up memory
        rm(obs.data.rep); gc()
      
      #########  estimation
      cat ("2. csk \n")
      
      ## skip.csk is set in C21.simulation_run.R (default is TRUE)
      ## if !skip.csk == TRUE, execute the code inside the block (meaning, skip.csk = FALSE)
      if (!skip.csk) {
        cat ("  2. csk - Policy estimation \n")
        
        
        # new package dtrSurv
        models = as.formula(Surv(`T`, delta) ~ A + lB + Z1 + Z2 + Z3 + Z4 + Z5)
          
        
        ## create new list of arguments to use dtrSurv
        ## used the observed data generated earlier (once it's in observable form)
        arg.csk2 = list(data = data.df, 
                        ## will give A_1, A_2, ... A_n.stages
                       txName = paste("A", 1:n.stages, sep = "_"),
                       
                       stageLabel = "_", 
                       
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
                       
                       ## using stratified analysis
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
        optimal.csk <- do.call(IHdtrSurv, c(arg.csk2, list(nodeSize = nodesize, minEvent = mindeath )))
        
        ## if the first element of the class == try-error, this value will be TRUE
        csk.error <- class(optimal.csk)[1] == "try-error"
        
        ##if csk.error is FALSE (code works properly), then assign the results of the dtrSurv() which are optimal, to the "policy" slot of arg.csk
        ## this is of class DTRSurv
        arg.csk$policy <- if (!csk.error) optimal.csk
        
        ## remove the results from dtrSurv()
        rm(optimal.csk); gc()
        
        cat ("  2. csk - Evaluation \n")
        ## setting a different seed
        set.seed(sim*10000 + 10)
        
        ## generate new data if csk.error is FALSE (aka code runs properly)
        ## use arguments provided in original arg.csk to generate new multistage data-- with POLICY added
        ### this is the same set of arguments we used to generate multistage data for the observed data
        ### basically, we predict the patient's optimal action after feeding them through the trained RF 
        if (!csk.error) csk.data.rep <- do.call(multiStageDynamics, arg.csk)
        
        
        ## val.fn depends on the criterion employed 
        ## observed value for mean = take mean of the total event time (across all patients)
        ## in this case, the critical value will be null
        ## this is the observed policy from the data we generated
        ## if criterion is "surv.prob", "surv.mean", count the proportion of elements >= crit.value (ignoring NAs)
        ## crit.value is put at 5 (5 year survival rate) in C21.simulation_run.R
        
        ## calculate value function of the new generated data-- generated based on predicted optimal policy
        if (!csk.error) result[sim, "csk"] <- val.fn(csk.data.rep$summary$cumulative.event.time)
        
        ## tt(2, reset = TRUE): resets timer for later measurements, and retreives the elapsed time in terms of mins for the time it takes to estimate and evaluate the policy
        result[sim, "time.csk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
        
        ## reset the policy to NULL and clean
        arg.csk$policy <- NULL; gc()
        
        ## re,pve the data generated
        rm(csk.data.rep); gc()
      }
      
                                      ########################################
                                      ## commented out since we don't use them 
                                      ########################################
        
    # cat ("3. Goldberg & Kosorok - lm \n")
    # if (!skip.gk) {
    #   cat ("  3. Goldberg & Kosorok - lm - Policy estimation \n")
    #   set.seed(sim*10000 + 2)
    #   optimal.gk.lm <-
    #     try(gk.separate(common.formula = formula.lm(p),
    #                     common.Tx.label = "A", stage.label = 1:n.stages, tau = tau,
    #                     data = data.df, stage.sep = "_", method = "lm", regress.prev.time = F))
    #   # optimal.gk.lm <- try(gk.Q(data = obs.data$output, est.fn = gk.lm, formula = formula.lm(p), tau = tau))
    #   gklm.error <- class(optimal.gk.lm)[1] == "try-error"
    #   arg.gk.lm$policy <- if (!gklm.error) optimal.gk.lm$survRF
    #   attr(arg.gk.lm$policy, "class") = "GKLM"
    #   rm(optimal.gk.lm); gc()
    #   
    #   cat ("  3. Goldberg & Kosorok - lm - Evaluation \n")
    #   set.seed(sim*10000 + 10)
    #   if (!gklm.error) gklm.data.rep <- do.call(multiStageDynamics, arg.gk.lm)
    #   if (!gklm.error) result[sim, "gkLM"]    <- val.fn(gklm.data.rep$summary$cumulative.event.time)
    #   result[sim, "time.gkLM"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
    #   arg.gk.lm$policy <- NULL; gc()
    #   rm(gklm.data.rep); gc()
    #   
    # cat ("4. Estimation - Goldberg & Kosorok  - rf \n")
    # 
    #   cat ("  4. Estimation - Goldberg & Kosorok  - rf - Policy estimation \n")
    #   set.seed(sim*10000 + 3)
    #   arg.gk.rf2 = list(common.formula = formula.rf(p),
    #                    common.Tx.label = "A", stage.label = 1:n.stages, tau = tau,
    #                    data = data.df, stage.sep = "_", regress.prev.time = F)
    #   optimal.gk.rf <-
    #     try(do.call(gk.separate, c(arg.gk.rf2, list(nodesize = nodesize, method = "rf"))))
    #   # optimal.gk.lm <- try(gk.Q(data = obs.data$output, est.fn = gk.lm, formula = formula.lm(p), tau = tau))
    #   gkrf.error <- class(optimal.gk.rf)[1] == "try-error"
    #   arg.gk.rf$policy <- if (!gkrf.error) optimal.gk.rf$survRF
    #   if (!gkrf.error) attr(arg.gk.rf$policy, "class") = "GKRF"
    #   rm(optimal.gk.rf); gc()
    #   
    #   cat ("  4. Estimation - Goldberg & Kosorok  - rf - Evaluation \n")
    #   set.seed(sim*10000 + 10)
    #   if (!gkrf.error) gkrf.data.rep <- do.call(multiStageDynamics, arg.gk.rf)
    #   if (!gkrf.error) result[sim, "gkRF"]    <- val.fn(gkrf.data.rep$summary$cumulative.event.time)
    #   result[sim, "time.gkRF"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
    #   arg.gk.rf$policy <- NULL; gc()
    #   rm(gkrf.data.rep); gc()
    # }
    # 
    # cat ("5. Estimation - Simoneau et al. \n")
    # if (!skip.dw) {
    #   cat ("  5. Estimation - Simoneau et al. - Policy estimation \n")
    #   set.seed(sim*10000 + 4)
    #   optimal.dw <- try(dwSurv(data = obs.data$output, tau = tau))
    #   dw.error <- class(optimal.dw)[1] == "try-error"
    #   arg.dw$policy <- if (!dw.error) optimal.dw$dw.est
    #   rm(optimal.dw); gc()
    #   
    #   cat ("  5. Estimation - Simoneau et al. - Evaluation \n")
    #   set.seed(sim*10000 + 10)
    #   if (!dw.error) dw.data.rep <- do.call(multiStageDynamics, arg.dw)
    #   if (!dw.error) result[sim, "dw"]  <- val.fn(dw.data.rep$summary$cumulative.event.time)
    #   result[sim, "time.dw"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
    #   arg.dw$policy <- NULL; gc()
    #   rm(dw.data.rep); gc()
    # }
      
      cat ("6. Estimation - zero-order model\n")
      ## all pts treated identically regardless of their characterstics
      ## chosen regime yields most favorable outcomes across pt population based on model predictions
      ## it can be seen as equivalent to the standard of care (aka every pt receives standard of care)
      if (!skip.zom) {
        cat ("  6. zero-order model - Policy estimation \n")
        set.seed(sim*10000 + 5)
        optimal.zom <- try(dtrSurv:::dtrSurv(data = data.df, 
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
                                             minEvent = 1, nodeSize = 1e+4, # zero-order model
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
        arg.zom$policy <- if (!zom.error) optimal.zom
        rm(optimal.zom); gc()
        
        cat ("  6. zero-order model - Evaluation \n")
        set.seed(sim*10000 + 10)
        if (!zom.error) zom.data.rep <- do.call(multiStageDynamics, arg.zom)
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
    print(Sys.time())
    #saveRDS(list(statistics = result, settings = setting), filename)
    