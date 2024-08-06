

######## this is a function to assess convergence using results from IH.dtrSurv.R which fits the initial forest


IHdtrConv <- function(data,

                      ## a numeric to indicate the strata we are assessing for convergence
                      strata,

                      ## the input forest results from the previous iteration
                      prev.iteration = NULL, nDP, params, nTimes, models, mTry, ..., long_data,
                      prev_probs){


  ######
  ###### NOTE: also input long_data & output long data

  ###########
  ########### have a stage-wise eligibility matrix so that each patient has same # stages
  ########### for patients who are not eligible, their matrices should just be 0


  ## using the long data that was input
  data <- data

  prev.iteration <- prev.iteration

  nDP <- nDP

  params <- params

  nTimes <- nTimes

  models <- models

  mTry <- mTry

  long_data <- long_data




  ########## for the last stage (stage 10) we don't need to append so we do this outside the loop
  ## start at stage nDP:
  message("Convergence Strata ", strata, " starting prediction from stage ", nDP)

  ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
  ### if the current iteration is nDP - 2, we want to predict for nDP - 1
  ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

  ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
  response_var <- as.character(formula(prev.iteration@FinalForest@model)[[2]][-1])
  response_with_stage <- paste0(response_var, "_", (nDP))
  terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")
  terms_with_stage <- paste0(terms, "_", (nDP))
  updated_formula <- paste("Surv(", paste(response_with_stage, collapse = " , "), ") ~ ", paste(terms_with_stage, collapse = " + "))
  x = get_all_vars(updated_formula, data %>% filter(!!sym(paste0("strata", strata,"_",nDP)) == 1))
  ## remove stage suffixto use in prediction
  new_col_names <- gsub(paste0("_", (nDP), "$"), "", colnames(x))
  colnames(x) <- new_col_names

    # Initialize a matrix to store the combined results, with a matrix of NAs
  ## putting 0's in all locations-- ones that belong to certain stages will be overwritten
  pr_pooled <- matrix(0, nrow = nTimes, ncol = nrow(data)*nDP)

  message("cases in stage: ", dim(x)[1])


  ##### we only want to run this prediction if there's available data, otherwise skip this step

  if (dim(x)[1] != 0){
  last.stage.pred <- .Predict(object = prev.iteration@FinalForest,
                              newdata = x,
                              params = params,
                              findOptimal = T)



  ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
  # Iterate over each stage of interest

  # Extract optimal treatments for the current stage
  optimal_treatments <- last.stage.pred$optimal@optimalTx

  # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
  ### patients who are eligible are ones who have a complete case for x
  eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

  # also update "A" column
  ### we select the eligible values which are in the current stage, and strata
  ### then we subset t only those eligible (have complete cases)
  long_data$A[long_data$stage == nDP & long_data[[paste0("strata", strata)]] == 1] [which(eligibility == 1)] <- optimal_treatments


  #### now, we want to extract the eligibility for the whole stage
  #### this means they have to BOTH have compelte cases, AND have strata1 == 1
  ### this should now be the same length as "eligibility"
  wholestage.eligibility <- as.numeric(
    ## within stage 10 variables, which are complete cases?
    ifelse(apply(get_all_vars(updated_formula, data), 1, function(x) all(!is.na(x))) &

                                                 ## within stage 10, which patients are in strata 1?
                                              long_data %>% filter(stage == nDP)%>% pull(!!sym(paste0("strata", strata))) == 1 ,
           ## if both conditions are true, put 1; otherwise, put 0
                                              1, 0) )

  ### next, create a matrix with a column for EACH patient, and merge these into the pooled matrix
  #### for pooled matrix, should have dim of wholestage.eligibility x nDP
  ### this way, we can sequence through these probs and each patient will have same # of columns



  ##########################################
  ## putting optimal predictions together ##
  ##########################################

  # Initialize the shifted probability matrix. each patient should have one column (instead of just eligible stages)
  ## we overwrite the eligible patients

  shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(wholestage.eligibility))

  ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
  ## do this only for eligible patients within the stage who were both in strata 1 and have complete cases
  shiftedprob1[, which(wholestage.eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)

  ## calculate the change in survival probability at each consecutive pair of time points
  ## subtracting each row of survShifted from the row above it
  ## append a row of 0s at the end to align with matrix dimensions
  ## probability mass vector representing the change in survival probabilities at each time point
  shiftedprob_done <- shiftedprob1 - rbind(shiftedprob1[-1L,], 0.0)

  ## sets very small values in pr to 0

  shiftedprob_done[abs(shiftedprob_done) < 1e-8] <- 0.0


  ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage
  ## add the shifted optimals into the pooled optimals together
  #### each patient will have the same number of stages by design, but ones that aren't included in the strata will just have 0s


  # Define the index to insert columns from the previous optimal
  ## this should be for the last stage
  insert_index <- seq(from = nDP, to = ncol(pr_pooled), by = nDP)

  # Initialize a counter for columns from append1_pr
  col_counter <- 1

  # Loop over each index to insert columns from append1_pr into combined_results
  for (w in seq_along(insert_index)) {
    # Insert columns from append1_pr into combined_results
    pr_pooled[,insert_index[w]] <-
      shiftedprob_done[, w]

    # Increment the counter for columns from append1_pr
    col_counter <- col_counter + 1
  }

  }


  for (i in (nDP - 1):(1)) {

    ## start at stage nDP:
    message("Convergence Strata ", strata, " starting prediction from stage ", i)



    ### first, we predict for the previous stage's optimal survival probability by feeding it through the forest
    ### if the current iteration is nDP - 2, we want to predict for nDP - 1
    ## then, for stage nDP - 2 (stage 23), predict the optimal survial and treatment for stage 24

    ####
    #### run for strata == 1
    ####

    ## then, use this previous predicted survival probability and append nDP - 2 observed info onto stage 24 predicted
    response_var <- as.character(formula(prev.iteration@FinalForest@model)[[2]][-1])
    response_with_stage <- paste0(response_var, "_", (i+1))
    terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", (i+1))
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = " , "), ") ~ ", paste(terms_with_stage, collapse = " + "))


  ### NOTE: we only predict for those patients who have the same next strata, and complete cases in that strata

    x = get_all_vars(updated_formula, data %>% filter(!!sym(paste0("strata", strata,"_",(i+1) )) == 1
                                                      & !is.na(!!sym(paste0("T","_",(i+1) )))
                     ))



    ### if there are no observations in x, we skip this whole thing
    ### there should automatically be 0's in pr_pooled for this
    ### meaning, we only run the rest of this loop if there are non-0 values for x



    if (dim(x)[1] != 0){

    ## remove stage suffixto use in prediction
    new_col_names <- gsub(paste0("_", (i+1), "$"), "", colnames(x))
    colnames(x) <- new_col_names
    last.stage.pred <- .Predict(object = prev.iteration@FinalForest,
                                newdata = x,
                                params = params,
                                findOptimal = T)

    ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
    # Iterate over each stage of interest

    # Extract optimal treatments for the current stage
    optimal_treatments <- last.stage.pred$optimal@optimalTx

    # We need an index of eligibility for the current stage at i+1 -- AKA the patients who are present in the last stage get 1, otherwise they get a 0
    ### patients who are eligible are ones who have a complete case for x
    eligibility <- as.numeric( ifelse(apply(x, 1, function(x) all(!is.na(x))), 1, 0) )

    # also update "A" column


    long_data$A[
      long_data$stage == (i + 1) &
        long_data[[paste0("strata", strata)]] == 1 &
        !is.na(long_data$T)][which(eligibility == 1)] <- optimal_treatments



## the dimensions of this are for all patients

    wholestage.eligibility <- as.numeric(
      ## within stage 10 variables, which are complete cases?
      ifelse(apply(get_all_vars(updated_formula, data), 1, function(x) all(!is.na(x))) &

               ## within stage 10, which patients are in strata 1?
               long_data %>% filter(stage == (i+1))%>% pull(!!sym(paste0("strata", strata))) == 1 ,
             ## if both conditions are true, put 1; otherwise, put 0
             1, 0) )


    ##########################################
    ## putting optimal predictions together ##
    ##########################################


    # Initialize the shifted probability matrix
    ## we overwrite the eligible patients

    shiftedprob1 <- matrix(0, nrow = nTimes, ncol = length(wholestage.eligibility))

    ### we retrieve the predicted survival curves for the patients in last.stage.pred@optimal@optimalY
    ## do this only for eligible patients
    shiftedprob1[, which(wholestage.eligibility == 1) ] <- t(last.stage.pred$optimal@optimalY)

    ###
    ### end strata == 1
    ###




    ###################################################
    ################################################# NEXT START HERE
    ###############################################################

    ## now, we want insert these predictions into the first column of a combined probability with optimal predictions from each stage

    #################### appending the previous stage
    ##### we want to make sure the appended stage is also in strata 1


    x_append1 <- stats::model.frame(formula = models,
                                    ## ## we want to exclude the A.opt.HC column and A.pool1 column, and only consider data from the prev timepoint
                                    data = long_data %>% filter(stage == i) %>% dplyr::select(-matches("^A\\.")),
                                    na.action = na.pass)

    ## current stage's eligibility based on being in strata1 and having complete cases
    ## this dimension should be the same as the total number of pts


    elig_append1 <- long_data %>%
      # Filter the data to include only rows where stage is equal to i
      filter(stage == i) %>%

      # Mutate the data to add a new column 'eligibility'
      mutate(
        eligibility = as.numeric(
          # Use ifelse to check two conditions for assigning eligibility
          ifelse(
            # Condition 1: Check if the strata column (constructed dynamically) is equal to 1
            !!sym(paste0("strata", strata)) == 1 &
              # Condition 2: Check if the row has complete cases excluding columns starting with "A"
              complete.cases(dplyr::select(., -matches("^A\\.|^gamma"))),
            # If both conditions are TRUE, assign 1
            1,
            # Otherwise, assign 0
            0
          )
        )
      ) %>%

      # Extract the 'eligibility' column as a vector
      pull(eligibility)

    # extract response and delta from model frame

    ## extract survival response for all 300 patients
    ### however, those that are not in elig_append1 = TRUE will receive an NA
    response_append1 <- stats::model.response(data = x_append1)
    ## if there are any that are not eligible for this stage, we turn into NA
    response_append1[!elig_append1, ] <- NA

    ## extract censoring indicator (delta) from the second column of the "response" data
    ## "L" is used to indicate that 2 is an integer
    delta_append1 <- response_append1[, 2L]
    delta_append1[!elig_append1] <- NA

    ## updates the "response" variable to only include the first column of the original "response" data which represents survival times
    response_append1 <- response_append1[, 1L]
    response_append1[!elig_append1] <- NA

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

    ## if there are no eligible cases, skip the following loop

    ## displays message indicating the number of cases that are still eligible for analysis at this stage

    message("cases in stage: ", sum(elig_append1))



    if (sum(elig_append1) != 0L){


    ## retrieve the number of timepoints from the params object
    ## defined in class_TimeInfo.R

    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases in the same strata

    ## NOTE: transpose is not needed bc shiftedprob has nrows = timepoints, ncols = patients
    ## only for the patients who were still eligible in the prior stage, we take their optimal shifted probabilities and overwrite them in the current matrix
    ## with columns with values that are non-0
    ### NOTE: we don't have the same number of patients in stage 9 as stage 10,
    ## so, we need to create an eligibility the same length as the patients in stage 10, with TRUE if they were in stage 9 too
    ## dimensions are same as the total number of patients

    prev.stag.elig <- long_data %>%
      filter(stage == (i+1)) %>%
      mutate(eligibility = as.numeric(ifelse(
        !is.na(T) & !!sym(paste0("strata", strata)) == 1,
        1,
        0
      ))) %>%
      pull(eligibility)


    ## initializes a survival matrix called "survMatrix" with 0's


    survMatrix <- matrix(data = 0.0,
                         nrow = nTimes,
                         ncol = length(prev.stag.elig))

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



    survMatrix[, which(prev.stag.elig == 1)] <-  t(last.stage.pred$optimal@optimalY)


    # shift the survival function down in time (T_i - Tq) and
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    ## transforms the survival functions in survMatrix to probability mass format for the current stage
    ## shifts the survival function based on the observed survival times
    ## this uses the survival matrix updated with the eligible patients & their optimal times from the last stage
    ## NOTE, we only shift patients who HAD previous stages


    append1_pr_1 <- .shiftMat(
      timePoints = .TimePoints(object = params),

      ## extracts columns from survMatrix corresponding to cases that are eligible
      ## this is a matrix matrix where each column represents survival function for an individual
      survMatrix = survMatrix[, elig_append1, drop = FALSE],

      ## extracts survival times corresponding to eligible cases
      ## this is how much to shift survival function for each individual
      shiftVector = as.matrix(response_append1[elig_append1]),

      ## probably transforming survival times into probabilities?
      surv2prob = TRUE
    )



    #############################
    #############################
    #############################
    ####### VISUALIZATION #######

    ## shifted probability output after adding new point for nDP - 3


    #    ## nDP
    #    y1 <- .shiftMat(
    #      timePoints = .TimePoints(object = params),
    #
    #      ## extracts columns from survMatrix corresponding to cases that are eligible
    #      ## this is a matrix matrix where each column represents survival function for an individual
    #      survMatrix = survMatrix[, elig_append1, drop = FALSE],
    #
    #      ## extracts survival times corresponding to eligible cases
    #      ## this is how much to shift survival function for each individual
    #      shiftVector = response_append1[elig_append1],
    #
    #      ## probably transforming survival times into probabilities?
    #      surv2prob = FALSE
    #    )[, 1]
    #
    #
    #
    #    plot(xaxis, y1)
    #    title("Convergence 2: stage 1 appending to stage 2 prediction")
    #

    #############################
    #############################
    #############################
    #############################

    ## we also have new patients, where this is their last stage. so for patients who are eligible in stage 24, but not stage 25, we need to predict
    ## if the new stages eligibility (elig_append1) doesn't match the previous stage eligibility (prev.stag.elig), we put a TRUE. Otherwise if they match, put FALSE
    ## this means if new stage eligibility is TRUE, and previous stage is FALSE, we put a value of TRUE
    newpt_elig <- ifelse(elig_append1 == T & prev.stag.elig == F, TRUE, FALSE)

    ## for these new patients, we predict their survival curves to get a stub
    response_with_stage <- paste0(response_var, "_", i)
    terms <- attr(terms(prev.iteration@FinalForest@model), "term.labels")
    terms_with_stage <- paste0(terms, "_", i)
    updated_formula <- paste("Surv(", paste(response_with_stage, collapse = " , "), ") ~ ", paste(terms_with_stage, collapse = " + "))

    x = get_all_vars(updated_formula, data)[which(newpt_elig == 1), ]

## if there are variables to predict, we move forward

    if (dim(x)[1] != 0) {
      ## remove stage suffixto use in prediction
      new_col_names <- gsub(paste0("_", i, "$"), "", colnames(x))
      colnames(x) <- new_col_names
      last.stage.pred <- .Predict(object = prev.iteration@FinalForest,
                                  newdata = x,
                                  params = params,
                                  findOptimal = T)

      ## the stage results: use this to assign optimal treatment to A.opt.HC based on the final stage prediction
      # Extract optimal treatments for the current stage
      optimal_treatments <- last.stage.pred$optimal@optimalTx

      # also update "A" column
      ## dimensions wise, this should include all the patients in the stage
      long_data$A[long_data$stage == i][which(newpt_elig == 1)] <- optimal_treatments

      ### now, we also insert these predicted probs into the shifted probabilities
      ## we want the indexes of the current stage's eligible patients, and find the indexes in the previous stage's and find where these don't match (0)
      ## for these columns, we insert the just predicted optimal for the patients who start in this stage

      newpt_shiftprob <- t(last.stage.pred$optimal@optimalY)

      ## convert into survival probabilities
      newpt_shiftprob <- newpt_shiftprob - rbind( as.matrix(newpt_shiftprob[-1L,]), 0.0)

      ## sets very small values in pr to 0

      newpt_shiftprob[abs(newpt_shiftprob) < 1e-8] <- 0.0


      ## now substitute into the append1_pr_1 matrix
      ### this now includes both shifted probabilities (double stubs) from previous stage as well as stubs from patients where this is their first stage
      ## select columns in append1_pr_1 where prev.stag.elig is 0 (didn't have a previous stage) for TRUE indices of elig_append1
      ## meaning, where the patient is eligible in the current stage, but this is their stub (no previous visit)
      append1_pr_1[, prev.stag.elig[which(elig_append1)] == 0] <- newpt_shiftprob

    }


    ## sets very small values in pr to 0

    append1_pr_1[abs(append1_pr_1) < 1e-8] <- 0.0

    #### now, we want to extract the eligibility for the whole stage
    #### this means they have to BOTH have compelte cases at stage i, AND have strata1 == 1
    ### this should now be the same length as "eligibility"
    ##### NOTE, the updated_formula gets overwritten from stage i + 1 to stage i variables

    wholestage.eligibility <- as.numeric(
      ## within stage 10 variables, which are complete cases?
      ifelse(apply(get_all_vars(updated_formula, data), 1, function(x) all(!is.na(x))) &

               ## within stage 10, which patients are in strata 1?
               long_data %>% filter(stage == i)%>% pull(!!sym(paste0("strata", strata))) == 1 ,
             ## if both conditions are true, put 1; otherwise, put 0
             1, 0) )

    ### now, create an overall matrix for all patients, inserting these appended probabilities only for patients who are eligible at the appended stage

    ## now, we put the appended probabilities in context of the full matrix of patients, with 0s in all the other locations
    ## AKA we overwrite these shifted probabilities for eligible patients
    ## number of columns is the number of patients for this stage (nDP - 2)
    ### NOTE: we need to make sure pts have the same number of stages so that it's easier to index by stage

    append1_pr <- matrix(0, nrow = nTimes, ncol = length(wholestage.eligibility))
    ## we overwrite the patients who are eligible in new stage (nDP - 3)
    append1_pr[,which(wholestage.eligibility == 1)] <-append1_pr_1


### append1_pr now contains all the probabilities (stubs + double stubs) for EACH patient in stage i


    # Define the index to insert columns from the previous optimal
    ## this should be for the last stage
    insert_index <- seq(from = i, to = ncol(pr_pooled), by = nDP)

    # Initialize a counter for columns from append1_pr
    col_counter <- 1

    # Loop over each index to insert columns from append1_pr into combined_results
    for (j in seq_along(insert_index)) {
      # Insert columns from append1_pr into combined_results
      pr_pooled[,insert_index[j]] <-
        append1_pr[, j]

      # Increment the counter for columns from append1_pr
      col_counter <- col_counter + 1
    }

    }
    }
  }

  ## for stage 1, we need to include the optimal predicted actions



  ## after we predict for all the stages, we fit a random forest to pool all stages

  # Define the dynamic variable name
  forest.name <- paste0("s2.strata", strata)

  # Perform the desired operation and assign the result to the dynamic variable
  assign(forest.name, .dtrSurvStep(
    model = models,
    data = long_data %>% filter(!!sym(paste0("strata", strata)) == 1 & !is.na(T)),
    priorStep = NULL,
    params = params,
    txName = "A",
    mTry = mTry[[nDP]],
    sampleSize = 1,

    ### NOTE: POOL1 = F in the main script IH.dtrSurv.R.... let's see if anythign changes if we set this to F
    pool1 = FALSE,

    ## we are inputting survival probabilities
    appendstep1 = TRUE,
    inputpr = pr_pooled[, colSums(pr_pooled) != 0]
  ))

  # Set the column name in long_data
  long_data$A.final <- NA

  # Extract eligibility for the current forest
  eligibility_final <- get(forest.name)@eligibility

  # Step 3: Insert pool1 results into A.pool1 for stages 3, 4, and 5
  ## also update the "A" column

  ###################
  ################### CHECK-- this might not match dimensions as eligibility is only for the pts in the strata
  ################### might need to subset the pts in the strata to update the actions
  ###################

  ## we want to filter for all the patients who are in the strata and don't have an NA


  long_data$A.final[long_data[[paste0("strata", strata)]] == 1 &
      !is.na(long_data$T)][which(eligibility_final == 1)] <- get(forest.name)@optimal@optimalTx

  long_data$A[long_data[[paste0("strata", strata)]] == 1 &
      !is.na(long_data$T)][which(eligibility_final == 1)] <- get(forest.name)@optimal@optimalTx


  ####################
  ####################
  ####################



  ######## now, we want to output this into a grid, so that all patients have all stages present
  ### for stages that are not in this strata, they will just take on values of 0
  #### for this, we need to have an overall eligibility vector across ALL stages: complete case & strata criteria
  #### we count all variables that begin with the variables in the model and count the complete cases
  # first, we update the formula to get rid of any underscores
  overall.form <- gsub("_.*?(?=\\s|\\)|$)", "", updated_formula, perl = TRUE)

## this is an eligibility across ALLL pts & time points, and should be an eligibility for pts in the strata
  allstage.eligibility <- as.numeric(

    ifelse(apply(get_all_vars(overall.form, long_data), 1, function(x) all(!is.na(x))) &

             ## which patients are in strata 1?
             long_data %>% pull(!!sym(paste0("strata", strata))) == 1 ,
           ## if both conditions are true, put 1; otherwise, put 0
           1, 0) )

  # Initialize the shifted probability matrix
  ## ## each patient will have k + 1 --> k stages; nrow(data) is the number of patients
  ## we overwrite the eligible patients
  shiftedprobfinal <- matrix(0, nrow = nTimes, ncol = length(allstage.eligibility) )
  shiftedprobfinal[,which(allstage.eligibility == T)] <-t(get(forest.name)@optimal@optimalY)



  ## get the final stage's area under the curve: we only look at the columns of the matrix in the last stage
  #finalstagepr <- shiftedprobfinal[, seq(from = 1, to = ncol(shiftedprobfinal), by = nDP)]

  ### for areas, we need to select only the columns & stages where the patient was eligible


  ## the result after appyling the area function to each column of the matrix of survival probabilities
  areas <- apply(shiftedprobfinal[, colSums(shiftedprobfinal) != 0], 2, function(surv_prob_col) {
    area_under_curve(surv_prob_col, params@timePoints)
  })



  ## captures current function call, including function name and all arguments passed to it
  cl <- match.call()

  ## ensures name of called function is set to "dtrSurv"
  cl[[1L]] <- as.name("IHdtrSurv")

  conv_forest <- new(
    Class = "DTRSurv",
    "stageResults" = list(),
    "IHstageResults" = list(),
    "FinalForest" = get(forest.name),
    "value" = NULL,
    "call" = cl,
    "params" = params,
    "integral_KM" = areas,
    "n_it" = NA,
    "avgKM_diff" = matrix(nrow = 2, ncol = 2),
    "valueTrain_list" = list(),
    "long_data" = long_data,
    "prev_probs" = shiftedprobfinal
  )

  return(conv_forest)


}
