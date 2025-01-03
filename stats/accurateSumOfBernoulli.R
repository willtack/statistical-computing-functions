accurateSumOfBernoulli <- function(p) {
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
  
  lPMF <- accurateSumOfBernoulli(p1)
  rPMF <- accurateSumOfBernoulli(p2)
  
  # convolution of two probability mass functions without using convolve function
  n1 <- length(lPMF)
  n2 <- length(rPMF)
  n <- n1 + n2 - 1
  PMF <- rep(0, n)
  for (k in 1:n) {
    for (i in 1:n1) {
      j <- k - i + 1
      if (j < 1 || j > n2) {
        next
      }
      PMF[k] <- PMF[k] + lPMF[i] * rPMF[j]
    }
  }
  
  PMF <- pmax(pmin(PMF, 1), 0)
  
  return(PMF)
  
}

p.test <- c(0.1,0.2,0.3,0.4)
accurateSumOfBernoulli(p.test)