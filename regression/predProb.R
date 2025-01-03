predProb<- function(x) {
  # x: a vector of predicted probabilities
  # return: a vector of predicted probabilities
  
  # Compute the predictive probabilities of a fitted five-category logistic regression model

  p <- length(x) # Initialize the vector of predicted probabilities
  categories <- 5
  linear_pred <- length(categories) # Initialize the linear predictor vector

  # Compute the predicted probabilities for each category
  for (k in 1:categories) {
    intercept <- 2^(-k)  # Intercept for category k (2^-k)
    sum_term <- sum(2^abs(k - (1:p)) * x)  # Sum of the predictors terms
    linear_pred[k] <- intercept + sum_term  # Linear predictor for category k
  }
  
  
  # log-sum-exp trick 
  max_logit <- max(linear_pred) 
  shifted_logit <- linear_pred - max_logit  
  
  exp_preds <- numeric(categories)
  for (k in 1:categories) {
    exp_preds[k] <- exp(shifted_logit[k])  
  }
  
  probs <- exp_preds / sum(exp_preds) # Normalize the predicted probabilities
  
  return(probs)
}

# x <- c(0,0.1)
# predProb(x)
  
  
  
