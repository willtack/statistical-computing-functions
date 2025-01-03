constrainedPolynomialRegression <- function(p, y) {
  
  n <- length(y)
  X <- (1:n) / n  # Predictor values as Xi = i/n
  z <- (1:n) / n  
  beta <- numeric(p + 1)
  
  # Build z vector for least squares fitting of beta1
    # Compute factorial values outside loop first 
    facts <- sapply(1:p, function(i) factorial(i))
    xpows <- sapply(2:p, function(i) X^i)
    
    for (j in 2:p) {
      # Compute powers of X for the higher-order terms
      z <- z + xpows[, j-1] / facts[j]  # Scaled by the factorial of j
    }
    
  # Calculate beta1 using least squares formula
  beta[2] <- sum((y - mean(y)) * (z - mean(z))) / sum((z - mean(z))^2)
  
  # Calculate beta0
  beta[1] <- mean(y) - beta[2] * mean(z)
  
  # Compute other beta j terms
  for (j in 3:(p + 1)) {
    beta[j] <- beta[j - 1] / (j-1)
  }
  
  return(beta)
}

# p <- 3
# y <- c(0.5, 1, 1.5)
# constrainedPolynomialRegression(p, y)

