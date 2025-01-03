nonlinearLogisticRegression <-function(df){
  # implements the secant method to find the maximum likelihood estimates for the following non-
  # linear logistic regression model:
  # logit(p) = α*Xi -α^2*Zi^2 
  # where Yi is a binary response variable taking values 0 and 1; and Xi and Zi are two predictors. α is
  # the parameter of interest and suppose that we know the true value of α is in (−5, 5)
  
  # df: a data frame with columns X, Z, and Y
  # returns estimate of α with five digit precision
  
  # f <- function(alpha){
  #   p <- exp(alpha*df$X - alpha^2*df$Z^2)/(1+exp(alpha*df$X - alpha^2*df$Z^2))
  #   return(sum(df$Y*log(p) + (1-df$Y)*log(1-p)))
  # }
  
  # score function
  f <- function(alpha) {
    eta <- alpha * df$X - alpha^2 * df$Z^2
    p <- 1 / (1 + exp(-eta))
    score <- sum((df$Y - p) * (df$X - 2 * alpha * df$Z^2))
    return(score)
  }
  
  # heuristics for choosing initial values
  
  # heuristic 1: fixed values
  # x0 <- -1; x1 <- 2
  
  # heuristic 2: mean of covariates if fall in range, 0,1 if not
  x0 <- mean(df$X) - 2.5
  x1 <- mean(df$X) + 1
  if(x0 < -5 | x0 > 5){x0 <- -2.5}
  if(x1 < -5 | x1 > 5){x1 <- 1}
  
  # heuristic 3: avg of estimated effects
  # k = 1
  # x0 = mean(df$Y)/mean(df$X) - k * mean(df$Z^2)
  # x1 = x0 + 2*sign(x0)
  # 
  # print(x0)
  # print(x1)

  # perform secant method
  max_iter <- 1000
  tol <- 1e-10
  sec_results <- secant(f, x0, x1, tol, max_iter)
  
  if (sec_results$convergence == 1){
    print("Convergence not reached")
  } else 
    print(sec_results$root)
  }

#' secant() - secant method for root finding
#' @param f : objective function to find root
#' @param x0, x1 : initial starting points
#' @param tol : absolute difference of stepwise difference in x for convergence
#' @param max_iter : maximum number of iterations
#' @return A list containing the following attributes:
#'    * root - x value with f(x) close to zero
#'    * f_root - f(root)
#'    * iter - number of iterations to reach the solution
#'    * convergence - 0 if the root was found successfully, 1 if not found
secant <- function(f,x0,x1,tol=1e-10,max_iter=1000){
  convergence = 1
  f0 = f(x0); f1 = f(x1)
  if(abs(f0-f1)<tol){ ## if x1 and x0 are too close, approx derivative is near-zero.
    warning("Expect a huge jump!")
    break
  }
  x12 <- -f1/(f1-f0)*(x1-x0) ## change in x
  x2 <- x1 + x12             ## interpolation
  for(iter in 1:max_iter){
    if(abs(x12)<tol){
      convergence = 0 ## convergence has reached
      break
    }
    f0 <- f1
    x1 <- x2
    f1 <- f(x2)
    f01 <- f1 - f0 ## difference of f(x)
    if(abs(f01)<tol){
      warning("Expect a huge jump!")
      break
    }
    x12 <- -f1/f01*x12 ## change in x - see equation in slide
    x2 <- x1 + x12     ## update rule.
  }
  return(list(root=x2, f_root = f(x2), iter=iter, convergence=convergence))
}


# pathname <- "/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/learning/data/nonlinearLogisticRegression/test.1.tsv"
# df = read.table(pathname, header=TRUE)
# print(head(df), row.names=FALSE)
# nonlinearLogisticRegression(df)
