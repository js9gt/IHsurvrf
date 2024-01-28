

set.seed(24)


# Function to simulate data for n different patients
simulate_patients_vectorized <- function(n, max_stages = Inf) {
  # Initialize an empty list to store patient data frames
  all_patient_data <- vector("list", n)
  
  # Parameters for the failure time distribution
  a1 <- -0.3
  b1 <- 0.1
  z1<- -0.3
  p1 <- -1
  g1 <- -0.2
  r1 <- -0.8
  
  # Parameters for the time to next visit distribution
  a2 <- 1.2
  b2 <- -0.05
  z2 <- -2.5
  p2 <- 0.1
  g2 <- -2
  r2 <- -1
  
  # Fixed variable Tau
  Tau <- 1000
  
  # Function to generate Ak, haz_tk, and haz_uk without NAs and non-positive values
  generate_hazards <- function() {
    Ak <- sample(c(-1, 1), size = 1, replace = TRUE, prob = c(0.5, 0.5))
    Sk <- rnorm(1, mean = 2, sd = 1)
    haz_tk <- exp(a1 + b1 * Sk + z1 * nstages + p1 * cumulative_length + g1 * Ak + r1 * Ak * Sk * nstages * cumulative_length)/10
    haz_uk <- exp(a2 + b2 * Sk + z2 * nstages + p2 * cumulative_length + g2 * Ak + r2 * Ak * Sk * nstages * cumulative_length)/10
    
    return(list(Ak = Ak, Sk = Sk, haz_tk = haz_tk, haz_uk = haz_uk))
  }
  
  # Simulate data for each patient
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
    
    cumulative_length <- 0
    stage_lengths <- numeric()
    nstages <- 0
    stage <- 1
    
    while (TRUE) {
      hazards <- generate_hazards()
      
      # Generate failure time Tk and time to next visit Uk
      Tk <- rexp(n = 1, rate = hazards$haz_tk)
      Uk <- rexp(n = 1, rate = hazards$haz_uk)
      
      # Check if Tk and Uk are non-missing and meet the specified criteria
      if (!is.na(Tk) && !is.na(Uk) && Tk > 1 && 
          Tk <= 100 && 
          Uk > 1 && Uk <= 100
          ) {
        nstages <- nstages + 1
        stage_length <- min(Tk, Uk)
        cumulative_length <- cumulative_length + stage_length
        
        if (cumulative_length > Tau || nstages > max_stages) {
          break
        }
        
        stage_data <- data.frame(
          Stage = stage,
          Ak = hazards$Ak,
          Tk = Tk,
          Uk = Uk,
          Sk = hazards$Sk,
          Censoring = ifelse(Tk < Uk, 1, 0),
          Stage_Length = stage_length,
          Cumulative_Length = cumulative_length,
          nstages = nstages - 1
        )
        
        patient_data <- rbind(patient_data, stage_data)
        
        if (Tk < Uk) {
          break
        }
        
        stage <- stage + 1
      }
    }
    
    all_patient_data[[patient_index]] <- patient_data
  }
  
  return(all_patient_data)
}

# Simulate data for 50 patients with a maximum of 5 stages
simulated_data <- simulate_patients_vectorized(n = 6, max_stages = 20)

