

set.seed(24)


# Function to simulate data for n different patients
simulate_patients <- function(n) {
  # Initialize an empty list to store patient data frames
  all_patient_data <- list()
  
  # Loop through each patient
  # Initialize empty data frame to store simulation results
  for (patient_index in 1:n) {
    patient_data <- data.frame(
      Stage = integer(),
      Ak = integer(),
      Tk = numeric(),
      Uk = numeric(),
      Sk = numeric(),
      Censoring = integer(),
      Stage_Length = numeric(),
      Cumulative_Length = numeric(),
      nstages = integer()
    )
    
    # Parameters for the failure time distribution
    a1 <- 1.2
    b1 <- -0.05
    z1 <- -2.5
    p1 <- 0.1
    g1 <- -2
    r1 <- -1
    
    # Parameters for the time to next visit distribution
    a2 <- -0.3
    b2 <- 0.1
    z2<- -0.3
    p2 <- -1
    g2 <- -0.2
    r2 <- -0.8
    
    # Fixed variable Tau
    Tau <- 1000
    
    # Initialize stage counter
    stage <- 1
    
    # Initialize cumulative length variable
    cumulative_length <- 0
    
    # Initialize a vector to store stage lengths
    stage_lengths <- c()
    
    # Initialize nstages to 0: number of previous visits (part of state variable)
    nstages <- 0
    
    
    # Simulate patient data until the censoring condition is met
    while (TRUE) {
      
      
      # Initialize Sk, Ak, haz_tk, and haz_uk as NAs
      ### the hazard rates cannot be 0 or negative, so we will regenerate Ak until this is non-negative
      Sk <- NA
      Ak <- NA
      haz_tk <- NA
      haz_uk <- NA
      
      
      # we repeat this re-generation of the hazards and therefore Tk and Uk, until the values of Tk and Uk are 1 or more
      # we also add an upper bound to make sure Tk and Uk do not exceed 100
      
      repeat{
  
        # Loop to generate Ak, haz_tk, and haz_uk without NAs and non-positive values
        ## in order to have a lower bound, we also want to rates to be less than or equal to 1 in order to have generated days 1 or more
        ## we will keep generating if the haz_tk and haz_uk are greater than 1
        while (is.na(Ak) || is.na(haz_tk) || is.na(haz_uk) || haz_tk <= 0 || haz_uk <= 0 || haz_tk > 15 || haz_uk > 15) {
          
          Ak <- sample(c(-1, 1), size = 1, replace = TRUE, prob = c(0.5, 0.5))
          
          # Generating State variable from N(1, 2) 
          Sk <- rnorm(1, mean = 2, sd = 1)
          
          
          ### stuff for propensity score generation
          #if (stage == 1) {
            
            #browser()
            
            # For the first stage, randomly assign treatment
            #Ak <- sample(c(0, 1), size = 1, replace = TRUE, prob = c(0.5, 0.5))
          #} else {
            # For subsequent stages, estimate propensity scores
            
            #browser()
            
            
            #propensity_model <- glm(Ak ~ Sk + nstages + cumulative_length, family = binomial(link = "logit"), data = patient_data)
            
            # Use Propensity Score Model to Predict Treatment using the current data
            #propensity_scores <- predict(propensity_model, newdata = data.frame(Sk, nstages, cumulative_length), type = "response")
            
            # Assign Binary Treatment Based on Propensity Scores
            #threshold <- 0.5
            #Ak <- ifelse(propensity_scores > threshold, 1, 0)
         # }
          # Calculate hazard rates
          haz_tk <- exp(a1 + b1 * Sk + z1 * nstages + p1 * cumulative_length + g1 * Ak + r1 * Ak * Sk * nstages * cumulative_length)
          haz_uk <- exp(a2 + b2 * Sk + z2 * nstages + p2 * cumulative_length + g2 * Ak + r2 * Ak * Sk * nstages * cumulative_length)

          
          #browser()  # Add browser() for debugging
          
          }
        
        # generate failure time Tk and time to next visit Uk
        Tk <- rexp(n = 1, rate = haz_tk) * 100
        Uk <- rexp(n = 1, rate = haz_uk) * 100
        
        # Check if Tk and Uk meet the specified criteria
        if (Tk >= 1 && Tk <= 100 && Uk >= 1 && Uk <= 100) break
      }
      
      
      # # Generate failure time (Tk) and time to next visit (Uk)
      #  repeat {
      #    # Generate failure time (Tk) using the specified model
      #    haz_tk <- exp(a1 + b1 * Sk + z1*nstages + p1*cumulative_length + g1 * Ak + r1 * Ak *nstages*cumulative_length)
      #    Tk <- rexp(n = 1, rate = haz_tk)*10
      # 
      #    # Generate time to next visit (Uk) using the specified model
      #    haz_uk <- exp(a2 + b2 * Sk + z2*nstages + p2*cumulative_length + g2 * Ak + r2 * Ak *nstages*cumulative_length)
      #    Uk <- rexp(n = 1, rate = haz_uk)*10
      # 
      #    # Check if Tk and Uk are both at least 1
      #    if (Tk >= 1 && Uk >= 1) break
      # }
      
      # Update nstages to reflect the number of previous stages
      nstages <- nstages + 1
      
      # Calculate stage length as the minimum of Tk and Uk, this is the number of days the stage occurred across
      stage_length <- min(Tk, Uk)
      
      # Update cumulative length (cumulative sum of stage lengths): only up to the beginning of stage k
      # does not include time spent in stage k
      cumulative_length <- cumulative_length + stage_length
      

      
      if (cumulative_length > Tau) {
        # If cumulative length exceeds Tau, stop without including data for the current stage
        break
      }
      
      # Create a data frame with the variables
      stage_data <- data.frame(
        Stage = stage,
        Ak = Ak,
        Tk = Tk,
        Uk = Uk,
        Sk = Sk,
        Censoring = ifelse(Tk < Uk, 1, 0),
        Stage_Length = stage_length,
        Cumulative_Length = cumulative_length,
        nstages = nstages - 1
      )
      
      # Add the data for the current stage to the patient_data
      patient_data <- rbind(patient_data, stage_data)

      
      # Stop if the minimum is Tau or if the minimum is less than Uk
      if (Tk < Uk ) {
        break
      }
      
      # Increment the stage counter
      stage <- stage + 1
    }
    # Add patient data to the list
    all_patient_data[[patient_index]] <- patient_data
  }
  
  return(all_patient_data)
}

# Simulate data for 50 patients
simulated_data <- simulate_patients(n =5)

write.csv(simulated_data, "/nas/longleaf/home/js9gt/survrf/Outputs/pt15.csv", row.names=FALSE)


print(sapply(simulated_data, function(patient_data) max(patient_data$Stage)))