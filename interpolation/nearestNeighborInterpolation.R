nearestNeighborInterpolation <- function(X, y, Z) {
  m <- nrow(Z)        # num rows of Z matrix
  p <- ncol(X)        # num cols in X matrix
  result <- numeric(m) # where to store results (interpoalted vals)
  
  X_norms <- rowSums(X^2)
  
  for (j in 1:m) {
    z <- Z[j, ]
    z_norm <- sum(z^2)
    
    # euclidean distance between z and each Xi
    distances <- X_norms + z_norm - 2 * X %*% z
    
    # Find the indices of the p+1 nearest neighbors
    nearest_indices <- order(distances)[1:(p + 1)]
    
    # construct matrix A and vector b
    A <- X[nearest_indices, ]
    b <- y[nearest_indices]
    A <- cbind(1, A) # add col of ones to A p+1 x p
    
    # wsolve the linear system eq (2)
    # beta <- solve(t(A) %*% A) %*% t(A) %*% b
    qr_decomp <- qr(A)
    beta <- qr.coef(qr_decomp, b)
    
    # Interpolated value at z, from eq (1)
    result[j] <- beta[1] + z %*% beta[-1]
  }
  
  return(result)
}

# X <- readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/hw6_data/test.1.X.rds")
# y <- readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/hw6_data/test.1.Y.rds")
# Z <- readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/hw6_data/test.1.Z.rds")
# gZ <- readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/hw6_data/test.1.gZ.rds")
# rst <- nearestNeighborInterpolation (X, y, Z)
# cor(rst, gZ)
