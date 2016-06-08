rm(list=ls())

# Newton-Raphson general root finding algorithm for pth root of x^(1/p)
# Let f(x) = x^p - S
# Then x_(n+1) = x_n - f(x) / f'(x)
# x_(n+1) = x_n - (x^p - S) / (p * x^(p-1))
# x_(n+1) = (p-1)*x/p + S/(p * x^(p-1)) 

root <- function(S, p, x=1, tol=1e-7) {
  # S: number for which the ith root is sought
  # i: root 
  # x: starting guess
  
  i <- 2
  ABS <- function(x) if (x < 0) x * -1 else x
  x <- as.vector(x)
  x[i] <- (p-1) * x[i-1]/p + S/(p * x[i-1]^(p-1))
  while (ABS(x[i] - x[i-1])/S > tol) {
    x[i+1] <- (p-1) * x[i]/p + S/(p * x[i]^(p-1))
    i <- i + 1
  }
  return(list(root=x[length(x)], iters=paste(i, "iterations"), x_i=x))
}