## source('streamedCov.R')
## Paste the function streamedCov(x, y) here
#' streamedCov (x, y): calculate covariance between vectors accurately while
#' streaming . The function must read each element of
#' x and y without reusing them after an iteration .
#' @param x - A numeric vector
#' @param y - Another numeric vector
#' @return covariance between x and y (with n - 1 as denominator )
#' https://chatgpt.com/share/ac28f8fa-d363-4f13-b04e-dd0d53123e0f
#' 
streamedCov <- function (x, y) {
  ## ensure that two input vectors are the same size
  stopifnot ( length (x) == length (y))
  
  n <- 0
  mx <- 0 # mean of x
  my <- 0 # mean of y
  cov.xy <- 0 # covariance of x and y
  
  for(i in 1: length (x)) {
    
    n <- n + 1 # increment n
    
    # obtain ith element of x and y vectors
    xi <- x[i]
    yi <- y[i]
    
    # calculate difference between xi and mx, yi and my
    # and update mx and my per west algo
    deltax <- (xi - mx)
    mx <- mx + deltax / n
    deltay <- (yi - my)
    my <- my + deltay / n
    
    adj <- (n-1) / n # west algo adjustment term 
    term <- adj * deltax * deltay # use covariance formula, eg https://en.wikipedia.org/wiki/Sample_mean_and_covariance
    cov.xy <- cov.xy + term
    
  }
  cov.xy <- cov.xy / (n - 1)
  return(cov.xy)
}