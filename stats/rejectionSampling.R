rejectionSampling <- function(n, a, b) {
  # Define the target density π(x) up to a normalizing constant
  target_density <- function(x) {
    exp(-x^2) * x^(a - 1) * (1 - x)^(b - 1)
  }
  
  # Use the Beta(a, b) distribution as the envelope density
  envelope_density <- function(x) {
    dbeta(x, a, b)
  }
  
  # Precompute the constant M to ensure the envelope dominates the target
  M <- optimize(function(x) target_density(x) / envelope_density(x), interval = c(0, 1), maximum = TRUE)$objective
  
  attempted <- 0
  accepted <- 0
  values <- numeric(n)
  
  while (accepted < n) {
    # Sample x from the Beta(a, b) envelope density
    x <- rbeta(1, a, b)
    # Sample u uniformly from (0, 1)
    u <- runif(1)
    # Acceptance criterion
    if (u < target_density(x) / (M * envelope_density(x))) {
      accepted <- accepted + 1
      values[accepted] <- x
    }
    attempted <- attempted + 1
  }
  
  # Return the results
  list(
    attempted = attempted,
    accepted = accepted,
    values = values
  )
}


####
system.time({samp <- rejectionSampling(5000000, 49, 10)})
print(samp$accepted)
print(samp$attempted)
print(samp$accepted/samp$attempted, digits=1) # print acceptance ratio
print(quantile(samp$values, probs=c(0.01, 0.25, 0.5, 0.75, 0.99), na.rm=TRUE), digits=1)
