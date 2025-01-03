## source('fastRidgeRegression.R')
## Paste the function fastRidgeRegression(X, Y, lambda) here
fastRidgeRegression <- function(X,Y,lambda) {
  
  # Computes the integer regression coeﬀicient estimates given X, y and λ.
  # X is the design matrix, Y is the response vector, and λ is the tuning parameter.
  # returns a data frame containing two variables index and beta
  # index is the index of the non-zero coeﬀicients and beta is the value of the non-zero integer coeﬀicients.
  
  ## copilot generated ##
  n = nrow(X) 
  p = ncol(X)
  #######################
  
  thr <- 0.1 * p

  # ensure lamba is greater than 0 
  stopifnot(lambda > 0 )
  
  # precompute ?
  XtY <- crossprod(X,Y)   # t(x) %*% y
  
  if ((p - n) > thr) {
    # Use the transformation to avoid direct computation of X^T X
    # beta = X^T * (XX^T + lambda * I)^-1 * Y
    XXt <- X %*% t(X)
    A <- XXt + lambda * diag(n)
    z <- solve(A, Y)
    beta <- t(X) %*% z
    
  } else if ( n < p ) {
    XtX <- crossprod(X)     # t(x) %*% x
    # Use SVD-based ridge regression
    res_svd <- svd(X)
    d <- res_svd$d

    # Compute beta
    UtY <- t(res_svd$u) %*% Y # t(U) %*% Y
    d_inv <- d / (d^2 + lambda) # Ridge regularization on singular values
    beta <- res_svd$v %*% (d_inv * UtY)
  }
  
  else {
    stopifnot(nrow(X) == length(Y))
    # compute the cholesky decomposition 
    XtX <- crossprod(X)     # t(x) %*% x  
    A <- XtX + lambda * diag(p)
    L <- chol(A)
    z <- forwardsolve(L,XtY, upper.tri=TRUE,transpose=TRUE) 
    beta <- backsolve(L,z)  # solve U %*% beta = z for beta
  }
  
  # round beta vector to nearest integer
  ## copilot generated ##
  beta <- round(beta)
   
  # create a dataframe with the index and beta
  ## copilot generated ##
  df <- data.frame(index = which(beta != 0), beta = beta[which(beta != 0)])
  return(df)
}


# Test 
# X = readRDS('/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/test.1.X.rds')
# Y = readRDS('/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/test.1.Y.rds')
# rst = fastRidgeRegression (X,Y,0.1)
# print(rst, row.names = FALSE )
# rst2 = fastRidgeRegression (X,Y,1000)
# print(rst2, row.names = FALSE )
