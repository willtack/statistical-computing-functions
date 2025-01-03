rmvnorm <- function(n, p, rho) {
  # Create covariance matrix Σ(ρ)
  i <- 1:p
  j <- 1:p
  Sigma <- outer(i, j, function(i, j) {
    exp(p * log(rho) * abs((i - j) / p)^1.99 - abs(cos(i / p)) - abs(cos(j / p)))
  })
  
  # Cholesky decomposition
  U <- chol(Sigma) # Lower triangular matrix
  
  # Generate samples in a single matrix operation
  Z <- matrix(rnorm(n * p), nrow = p, ncol = n) # Z ~ N(0, I)
  X <- crossprod(U,Z) # X ~ N(0, Σ)
  X <- t(X)
  
  # Calculate statistics
  # E[max_j X_j]
  row_max <- apply(X, 1, max) # Maximum value for each row
  E_max <- mean(row_max) # Average max across all samples
  
  # E[sqrt(sum(X_j^2))]
  E_norm <- mean(sqrt(rowSums(X^2)))
  
  # Pr(X1 * X2 > 0.5ρ)
  condition <- X[, 1] * X[, 2] > 0.5 * rho
  prob <- mean(condition) # Fraction satisfying the condition
  
  # Return the result as a vector
  return(c(E_max, E_norm, prob))
}

# print(rmvnorm(100000, 80, 0.42), digits = 1)
# print(rmvnorm(100000, 103, 0.63), digits = 1)


