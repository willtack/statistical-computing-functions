solveMatrixEquation <- function(A, B, C) {
  # Solve the matrix equation AX + XB = C for X
  # Inputs: A, B, C are sparse matrices
  # Output: X is a sparse matrix
  #
  # https://en.wikipedia.org/wiki/Sylvester_equation#Existence_and_uniqueness_of_the_solutions
  #
  # Rewrite as {\displaystyle (I_{m}\otimes A+B^{T}\otimes I_{n})\operatorname {vec} X=\operatorname {vec} C,}
  # where m = n in our case
  
  n <- nrow(A) # get dimensions
  I <- diag(n) # initialize identity matrix

  # transform equation into linear system with vectorized forms of matrices using kronecker products
  K1 <- kronecker(I, A)
  K2 <- kronecker(t(B), I)
  K <- K1 + K2
  
  # test if is positive definite
  # https://chatgpt.com/share/6717c866-7af8-8008-97c7-58ddcd0c2d3f
  is_pos_def <- tryCatch({
    chol(K) 
    TRUE
  }, error = function(e) {
    FALSE
  })
  

  if (is_pos_def) {
    # Cholesky solve if positive definite
    L <- chol(K)
    vec_X <- backsolve(L, forwardsolve(t(L), as.vector(C)))
  } else {
    # regular solve if not 
    vec_X <- solve(K, as.vector(C))
  }

  # reshape to matrix, round, and convert to sparse matrix
  X <- matrix(vec_X, nrow = n, ncol = n)
  X <- round(X)
  X_sparse <- as(X, "TsparseMatrix")

  return(X_sparse)
}

### TESTING ####
# library(Matrix)
# 
# A = readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/data/test.1.A.rds")
# B = readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/data/test.1.B.rds")
# C = readRDS("/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/mastery/data/test.1.C.rds")
# rst = solveMatrixEquation(A, B, C)
# print(data.frame(i=rst@i +1, j=rst@j +1, x=rst@x),row.names=FALSE)

