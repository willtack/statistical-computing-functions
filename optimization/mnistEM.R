## Implement your own E-M clustering algorithm
#' @param X : (n * 784) matrix of images (you may consider binarizing the matrix at X[i,j] == 128 if needed)
#' @param nclust : Number of clusters to generate
#' @param max.iter : Maximum number of iteration, may end early based on convergence criteria
#' @param tol : Parameters for determining convergence.
#' @return A list containing at least the following attributes:
#'    * mus - A [ ncol(X) * nclust ] matrix of mean pixel strength (0-1 scale)  for each cluster.
#'    * est - The best guess cluster of each image to one of the labels (i.e. in `1:nclust`).
mnistEM = function(X, nclust, max.iter = 100, tol=1e-6) {
  
  # Binarize the data
  X_bin <- ifelse(X > 127, 1, 0) # binarize the data
  n <- nrow(X_bin)
  p <- ncol(X_bin)
  
  # Initialize parameters
  #set.seed(123) 
  # z <- sample(1:nclust, n, replace = TRUE) # random cluster assignments
  # pi_k <- table(z) / n # mixing proportions
  # theta <- matrix(runif(nclust * p, 0.3, 0.7), nrow = nclust) # random initial probabilities
  
  # Initialize parameters using kmeans
  kmeans_init <- kmeans(X_bin, centers = nclust, nstart = 10)
  z <- kmeans_init$cluster
  pi_k <- table(z) / n
  theta <- tapply(1:n, z, function(idx) colMeans(X_bin[idx, ]))
  theta <- do.call(rbind, theta)
  
  log_likelihood <- -Inf #initialize log-likelihood
  
  for (iter in 1:max.iter) {
    ##### E-step: Compute responsibilities #####
    resp <- matrix(0, n, nclust) # Responsibility matrix
    for (k in 1:nclust) {
      prob <- X_bin %*% log(theta[k, ] + 1e-10) + (1 - X_bin) %*% log(1 - theta[k, ] + 1e-10)
      resp[, k] <- log(pi_k[k]) + prob
    }
    
    # Numerical stability adjustment
    max_log <- apply(resp, 1, max)
    resp <- exp(resp - max_log) # Convert back to probabilities
    resp <- resp / rowSums(resp) # Normalize responsibilities
    
    ##### M-STEP #######
    nk <- colSums(resp) # Effective number of points in each cluster
    pi_k <- nk / n # Update mixing proportions
    theta <- t(resp) %*% X_bin / nk # Update theta (cluster-wise probabilities)
    
    # Compute log-likelihood ---
    new_log_likelihood <- sum(log(rowSums(resp)))
    
    # Check for convergence
    if (abs(new_log_likelihood - log_likelihood) < tol) {
      break
    }
    # print(log_likelihood)
    log_likelihood <- new_log_likelihood
  }
  
  # Assign clusters based on maximum responsibility
  z <- apply(resp, 1, which.max)
  
  # Return results
  list(est = z, log_likelihood = log_likelihood, iter = iter, theta = theta, pi_k = pi_k)
  
}