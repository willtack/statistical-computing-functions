## source('kernelRidgeRegression.R')
## Paste the function kernelRidgeRegression(df, lambda, rho, bandwidth) here
kernelRidgeRegression <- function(df, lambda, rho, bw){
  ### partially copilot generated ###
  # Kernel Ridge Regression
  # df: data frame that contains 4 columns (x_train, y_train, x_test, y_test)
  # lambda: regularization parameter
  # rho: kernel parameter
  # bw: bandwidth parameter
  # return: predictive mean square errors
  
  # Ensure parameter space
  stopifnot(rho > 0)
  stopifnot(lambda > 0)
  
  x_train <- df$x_train
  y_train <- df$y_train
  x_test <- df$x_test
  y_test <- df$y_test
  
  # Kernel function: returns exp(-rho * (x - y)^2) if |x - y| <= bw, else 0
  kernel_matrix <- function(x1, x2, rho, bw) {
    # Initialize an empty matrix of size length(x1) x length(x2)
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(0, n1, n2)  # Initialize with zeros
    
    # create sparse matrix
    for (i in 1:n1) {
      for (j in 1:n2) {
        distance<- abs(x1[i] - x2[j])
        if (distance <= bw) {
          K[i, j] <- exp(-rho * distance^2) 
        }
      }
    }
    
    return(K)
  }
  
  
  # Construct the kernel matrix K
  K <- kernel_matrix(x_train, x_train, rho, bw)
  
  # Solve (K + λI) β = Y
  n <- length(y_train)
  alpha <- solve(K + lambda * diag(n), y_train)
  
  # Predict
  K_test <- kernel_matrix(x_test, x_train, rho, bw)  # construct test matrix
  y_pred <- K_test %*% alpha
  
  # Compute predictive mean square error (PMSE)
  pmse <- mean((y_test - y_pred)^2)
  
  return(pmse)
}