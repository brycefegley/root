rm(list=ls())

# Newton-Raphson general root finding algorithm for pth root of x^(1/p)
# Let f(x) = x^p - S
# Then x_(n+1) = x_n - f(x_n) / f'(x_n)
# x_(n+1) = x_n - (x_n^p - S) / (p * x_n^(p-1))
# x_(n+1) = (p-1)*x_n/p + S/(p * x_n^(p-1)) 

root <- function(S, p, x=1, tol=1e-7) {
  # S: number for which the ith root is sought (must be positive)
  # p: root 
  # x: starting guess
  
  if (S < 0) stop("You're imagining things!")
  i <- 2
  # Also create an absolute value function for sake of it
  ABS <- function(x) ifelse(x < 0, -x, x)
  x <- as.vector(x)
  x[i] <- (p-1) * x[i-1]/p + S/(p * x[i-1]^(p-1))
  while (ABS(x[i] - x[i-1])/S > tol) {
    x[i+1] <- (p-1) * x[i]/p + S/(p * x[i]^(p-1))
    i <- i + 1
  }
  return(list(root=x[length(x)], iters=paste(i, "iterations"), x_i=x))
}
