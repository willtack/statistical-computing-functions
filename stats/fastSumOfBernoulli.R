fastSumOfBernoulli <- function(p) {
  # computes the probability mass function of sum of independent Bernoulli random variables using a
  # divide-and-conquer algorithm within O (n(log n)2) time complexity

  n <- length(p)
  
  # for single Bernoulli random variable
  if (n == 1) {
    return(c(1-p, p))
  }
  
  # divide vector in half
  m <- floor(n / 2)
  p1 <- p[1:m]
  p2 <- p[(m + 1):n]
  
  lPMF <- fastSumOfBernoulli(p1)
  rPMF <- fastSumOfBernoulli(p2)
  
  convolved <- convolve(lPMF, rev(rPMF), type = "open")
  
  # convolved[abs(convolved) < 1e-12] <- 0
  
  convolved <- pmax(pmin(convolved, 1), 0)

  
  return(convolved)
}

p.test <- c(1,1,1,1,1,1,1,1,1,1)
fastSumOfBernoulli(p.test)