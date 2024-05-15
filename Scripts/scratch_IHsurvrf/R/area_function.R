
### this function is used to find the area under the curve to track the integral of the survival functions for patients to see if they are converging
### this is used in IH.dtrSurv.R for the first pooling and last pooling
### also used in IH.pool1.R for each complete step of the pooling algorithm

# Define a function to calculate the area under the curve using trapezoidal rule
area_under_curve <- function(surv_prob, time_points) {
  # Sort time points and corresponding survival probabilities
  sorted_indices <- order(time_points)
  sorted_time_points <- time_points[sorted_indices]
  sorted_surv_prob <- surv_prob[sorted_indices]
  
  # Calculate the width of each trapezoid
  widths <- diff(sorted_time_points)
  
  # Calculate the average height of each trapezoid (average of adjacent survival probabilities)
  heights <- (sorted_surv_prob[-1] + sorted_surv_prob[-length(sorted_surv_prob)]) / 2
  
  # Calculate the area of each trapezoid
  areas <- widths * heights
  
  # Sum up the areas to get the total area under the curve
  total_area <- sum(areas)
  
  return(total_area)
}


