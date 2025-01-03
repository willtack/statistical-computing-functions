stationaryDistribution <- function(P) {
  
  n <- nrow(P) # number of rows
  p <- diag(n) - P # in parentheses
  # add in constraint that sum pi = 1
  A <- rbind(p, rep(1, ncol(P))) # add row of 1s to 
  b <- c(rep(0, n), 1) # add 1 to last element of zero vector on right hand side
  # solve for pi
  result <- qr.solve(A, b)

  return(result)
}

# P = readRDS ('/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/learning/data/stationaryDistribution/test.1.rds')
# dim(P)
# stationaryDistribution(P)
