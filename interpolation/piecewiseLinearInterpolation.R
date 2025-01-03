piecewiseLinearInterpolation <- function(X,Y,Z) {
  # X: x values of the data points
  # Y: y values of the data points, same size as X
  # Z: x values of the points to be interpolated
  # returns: data frame containing the following two attributes
  #           f.hat : A vector containing ˆf (zj ) for j ∈ {1, . . . , m}
  #           f.true : A vector containing f (zj ) for j ∈ {1, . . . , m} for the quadratic function f
  
  ## Calculate fhat
  coefs <- piecewise.linear.interpolation(X, Y)  # Calculate the coefficients of the piecewise linear function
  f.hat <- eval.piecewise.poly(Z, coefs)  # Evaluate the piecewise linear function at Z
  
  ## Calculate ftrue
  # Calculate the coefficients of the quadratic function
  A <- cbind(X^2, X, 1)  # Design matrix with X^2, X, and a constant term
  coefs1 <- solve(t(A) %*% A, t(A) %*% Y)  # Least squares solution
  
  a <- coefs1[1]  # Coefficient of x^2
  b <- coefs1[2]  # Coefficient of x
  c <- coefs1[3]  # Constant term

  # Compute the true quadratic values at Z
  f.true <- a * Z^2 + b * Z + c

  # data frame containing fhat and ftrue
  return(data.frame(f.hat = f.hat, f.true = f.true))

}


#' piecewise.linear.interpolation()
#' @param x A vector of observed x-coordinates
#' @param y A vector of observed y-coordinates
#' @return A dataframe containing the following attributes
#'    * x : x-values interpolated (2nd to last)
#'    * b : y-intercepts for each interval line
#'    * m : slope for each interval
piecewise.linear.interpolation  = function(x,y){
  n = length(x)
  od_x = order(x)
  y = y[od_x] # reorder y based on x values
  x = x[od_x] # reorder x in the same way
  m = diff(y)/diff(x) # slopes between two nearby points
  b = y[-1] - m*x[-1] # y-intercepts from 2nd to last points
  return(cbind(x=x[-1],b,m))
}

#' evaluate polynomial function
#' @param x - A vector of x-coordinates to evaluate
#' @param coefs - Polynomial coefficient c(intercept,1st,2nd,...)
#' @return A vector of predicted values from the polynomial
eval.poly <- function(x,coefs){
  y <- rep(0, length (x))    ## y0 = 0
  for(i in length(coefs):1L)
    y <- coefs[i] + x * y    
  return(y)
}

#' eval.piecewise.poly()
#' @param x     A vector of x-coordinates to be evaluated
#' @param coefs Data frames containing (x, b, m)
#' @return A vector of y-coordinates linearly interpolated
eval.piecewise.poly <- function(x,coefs){
  n = nrow(coefs)                     # n is number coefficients, 1 less than observed points
  x_bound = c(-Inf,coefs[-n,"x"],Inf) # Make n+1 values (-Inf,x_1,x_2,x_{n-1},Inf)
  y = rep(NA,length=length(x))        # initialized interpolated values as NA
  for(i in 1:n){
    # select x points that fall into i-th interval
    idx = which((x <= x_bound[i+1]) & (x > x_bound[i]))
    y[idx] = eval.poly(x[idx],coefs[i,c("b","m")]) # apply interpolation on the subsets
  }
  return(y)
}

# setwd("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/learning/data/piecewiseLinearInterpolation")
# X = readRDS('test.1.X.rds')
# Y = readRDS('test.1.Y.rds')
# Z = readRDS('test.1.Z.rds')
# rst = piecewiseLinearInterpolation(X,Y,Z)
# dim(rst)
# head(rst)