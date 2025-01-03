euclideanDistanceThreshold <- function(X,Y,thres){
  # returns the number of column pairs with Euclidean
  # distance less than or equal to thres from two matrices X and Y with the same number of rows
  # X and Y are matrices with the same number of rows
  # thres is a positive number
  # returns the number of column pairs with Euclidean distance less than or equal to thres
  # between the columns of X and Y
  # if X and Y have different number of rows, return NA
  stopifnot(nrow(X) == nrow(Y))

  sqNorm_X <- colSums(X^2) # square every element in each col of X and sum 
  sqNorm_Y <- colSums(Y^2) # square every element in each col of Y and sum
  
  # compute squared distances to avoid square rooting
  # foil out original formula
  distance_squared <- outer(sqNorm_X, sqNorm_Y, "+") - 2 * t(X) %*% Y
  
  # sum the number of distances that are less than or equal to the *squared* threshold
  res <- sum(distance_squared <= thres^2)

  return(res)
  
}

# Test the function
set.seed (1234)
X = matrix(rnorm (1e5 ), 100, 1000)
Y = matrix(rnorm (2e5 ), 100, 2000)
system.time(rst <- euclideanDistanceThreshold (X, Y, 15))