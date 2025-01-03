neuralNetworkLBFGSB <- function(p, df){
  # p: number of nodes
  # df: data frame with X and Y columns
  
  X <- df$X
  Y <- df$Y

  # GeLu function and its derivative
  GeLu <- function(x) {
    return(x * pnorm(x)) 
  }
  
  GeLu_derivative <- function(x) {
    return(pnorm(x) + x * dnorm(x)) # product rule, deriv of cdf is density 
  }
  
  calculate_residuals <- function(alpha) {
    alpha0 <- alpha[1]
    alpha_j <- alpha[2:(p + 1)]
    alpha_pj <- alpha[(p + 2):(2 * p + 1)]
    alpha_2pj <- alpha[(2 * p + 2):(3 * p + 1)]
    
    # Calculate the predicted values for each observation
    Y_pred <- alpha0 + rowSums(
      sapply(1:p, function(j) {
        alpha_j[j] * GeLu(alpha_pj[j] + alpha_2pj[j] * X)
      })
    )
    
    # Compute residuals
    residuals <- Y - Y_pred
    return(residuals)
  }
  

  # define the objective function to minimize from the neural net in the assignment
  objective_function <- function(alpha) {
    residuals <- calculate_residuals(alpha)
    
    # Compute mean squared error
    mse <- mean((residuals) ^ 2)
    return(mse)
  }
  
  # define the gradient function for the objective function
  gradient_function <- function(alpha) {
    alpha0 <- alpha[1]
    alpha_j <- alpha[2:(p + 1)]
    alpha_pj <- alpha[(p + 2):(2 * p + 1)]
    alpha_2pj <- alpha[(2 * p + 2):(3 * p + 1)]
    
    residuals <- calculate_residuals(alpha)
    
    # Initialize gradient vector
    grad <- numeric(3 * p + 1)
    
    # Compute gradients
    grad[1] <- -2 * mean(residuals)  # Gradient with respect to alpha0
    
    # Gradients with respect to alpha_j, alpha_pj, and alpha_2pj
    for (j in 1:p) {
      # pre-compute GeLU and its derivative
      input <- alpha_pj[j] + alpha_2pj[j] * X
      GeLu_val <- GeLu(input)
      GeLu_deriv <- GeLu_derivative(input)
      
      grad[1 + j] <- -2 * mean(residuals * GeLu_val)  # Gradient wrt alpha_j
      grad[1 + p + j] <- -2 * mean(residuals * alpha_j[j] * GeLu_deriv)  # Gradient wrt alpha_pj
      grad[1 + 2 * p + j] <- -2 * mean(residuals * alpha_j[j] * GeLu_deriv * X)  # Gradient wrt alpha_2pj
    }
    
    return(grad)
  }
  
  
  # initialize alpha vector randomly Uniform[-1,1]
  alpha_init <- runif(3 * p + 1, min = -0.5, max = 0.5) # copilot generated
  lower_bounds <- rep(-10, 3 * p + 1)
  upper_bounds <- rep(10, 3 * p + 1)
  
  # run the L-BFGS-B algorithm using optim()
  result <- optim(
    par = alpha_init,
    f = objective_function,
    gr = gradient_function,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(maxit = 5000)
  )
  
  return(result)
}


# Test the neuralNetworkLBFGSB function
# pathname <- "/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/hw7_data/test.1.tsv"
# df = read.table(pathname,header=TRUE)
# neuralNetworkLBFGSB(10,df)
