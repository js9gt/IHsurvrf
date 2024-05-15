



# Function to compute the predicted propensity based on the baseline state variable generated

## beta propensity has 0 for intercept and -1
## p will be the dimensions of the initial state vector that we use to fit the model
beta.propensity = function(p) c(0, -rep(0.5, p))

propensity_function <- function(state, p = 1) {
  beta_values <- beta.propensity(p)  # Get the beta values as a numeric vector
  ##prepare baseline state data to turn into a matrix
  #state <- data.matrix(state, rownames.force = NA)

  ## unname the dataframe

  #state <- unname(state)
  # Compute the predicted probability using the logistic function
  plogis(cbind(1, state) %*% beta.propensity(p))
}
