neuralNetworkNelderMead <- function(p, df) {

  alpha_init <- rep(0, 3 * p + 1)
  
  # run the modified Nelder-Mead algo
  result <- Nelder.Mead(
    f = objective_function,
    x0 = alpha_init,
    tol = 1e-5,
    max_iter = 10000,
    p = p,
    X = df$X,
    Y = df$Y
  )
  
  return(c(result$fmin, result$iter, result$eval_count))
}

# GeLU activation function
GeLu <- function(x) {
  x * pnorm(x)
}

# Objective function for the neural network model
objective_function <- function(alpha, p, X, Y) {
  n <- length(Y)
  predictions <- alpha[1] +
    rowSums(sapply(1:p, function(j) {
      alpha[j + 1] * GeLu(alpha[p + j + 1] + alpha[2 * p + j + 1] * X)
    }))
  sum((Y - predictions)^2) / n
}

# Modified Nelder-Mead function with tracking for function evaluations
Nelder.Mead <- function(f, x0, tol = 1e-10, max_iter = 1000, ...) {
  eval_count <<- 0
  d <- length(x0)   # d: dimension of the simplex
  X <- matrix(x0, nrow = d, ncol = d + 1)    # set d+1 simplex points
  X[,-(d+1)] <- X[,-(d+1)] + diag(d) # create a simplex
  Y <- apply(X, 2, function(x) { eval_count <<- eval_count + 1; f(x, ...) })
  
  # Initialize variables
  idx_max <- NULL; idx_min <- NULL; idx_2ndmax <- NULL
  mid_point <- NULL; tru_line <- NULL
  
  update.extremes <- function() {
    if(Y[1] > Y[2]) {
      idx_max <<- 1; idx_min <<- 2; idx_2ndmax <<- 2
    } else {
      idx_max <<- 2; idx_2ndmax <<- 1; idx_min <<- 1
    }
    if(d > 1) {
      for(i in 3:(d + 1)) {
        if(Y[i] <= Y[idx_min]) {
          idx_min <<- i
        } else if(Y[i] > Y[idx_max]) {
          idx_2ndmax <<- idx_max; idx_max <<- i
        } else if(Y[i] > Y[idx_2ndmax]) {
          idx_2ndmax <<- i
        }
      }
    }
  }
  
  update.mid.point <- function() {
    mid_point <<- apply(X[,-idx_max, drop = FALSE], 1, mean)
    tru_line <<- X[, idx_max] - mid_point
  }
  
  update.next.point <- function(step_scale) {
    next_point <- mid_point + tru_line * step_scale
    Y_next <- f(next_point, ...)
    eval_count <<- eval_count + 1
    if(Y_next < Y[idx_max]) {
      X[, idx_max] <<- next_point
      Y[idx_max] <<- Y_next
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  contract.simplex <- function() {
    X[,-idx_min] <<- 0.5 * (X[,-idx_min] + X[,idx_min])
    Y[-idx_min] <<- apply(X[,-idx_min], 2, function(x) { eval_count <<- eval_count + 1; f(x, ...) })
  }
  
  convergence <- 1
  for(iter in 1:max_iter) {
    update.extremes()
    
    if(abs(Y[idx_max] - Y[idx_min]) <= tol * (abs(Y[idx_max]) + abs(Y[idx_min]) + tol)) {
      convergence <- 0
      break
    }
    update.mid.point()
    
    update.next.point(-1.0)
    if(Y[idx_max] < Y[idx_min]) {
      update.next.point(-2.0)
    } else if(Y[idx_max] >= Y[idx_2ndmax]) {
      if(!update.next.point(0.5)) {
        contract.simplex()
      }
    }
  }
  
  list(
    xmin = X[, idx_min],
    fmin = Y[idx_min],
    convergence = convergence,
    iter = iter,
    eval_count = eval_count
  )
}


# pathname <- "/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/learning/data/neuralNetworkNelderMead/ex10testcase1/test.1.tsv"
# df = read.table(pathname,header=TRUE)
# head(df)
# neuralNetworkNelderMead (10, df)