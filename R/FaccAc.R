

FaccAc <- function(y, X, d, coef, lambda, th) {
  
  a = dim(X)[1]
  b = dim(X)[2]
  p = dim(X)[3]
  
  coef = matrix(coef)
  
  if (dim(lambda)[1] != a || dim(lambda)[2] != b) {
    stop("Dimension of y, X, d or lambda is incorrect.")
  }
  
  if (length(coef) != p) {
    stop("Dimension of coef is incorrect.")
  }
  
  La = matrix(0, a, b)
  Ypre = matrix(0, a, b)
  Xabs = matrix(0, a, b)
  D = rowSums(d)
  
  Z = array(0, c(a, b, p))
  U = array(0, c(a, b, p))
  W = matrix(0, a, b)
  
  for (i in 1:a) {
    for (j in 1:b) {
      La[i, j] = sum((y <= y[i, j]) * lambda)   
    }
  }
  
  for (k in 1:b) {
    Ypre[, k] = X[, k, ] %*% coef
  }
  
  for (k in 1:p) {
    Xabs[, ] = Xabs[, ] + abs(X[,, k])
  }
  
  A = 1/th + D
  C = 1/th + rowSums(La * exp(Ypre))
  YpreExp = exp(Ypre)
  
  for (k in 1:p) {
    temp1 = X[,, k] * YpreExp
    temp2 = abs(X[,, k]) * YpreExp * Xabs
    
    Z[,,k] = FaccCal(y, temp1, Z[,,k], A, C)
    U[,,k] = FaccCal(y, temp2, U[,,k], A, C)
  }
  
  W = FaccCal(y, YpreExp, W, A, C)

  # Update parameters for Gamma Frailty
  Q1 = a*(digamma(1/th) + log(th)-1)/(th^2) + sum(A/C - digamma(A) + log(C))/(th^2) 
  Q2 = a*(3 - 2*digamma(1/th) - 2*log(th))/(th^3) + 2*sum(digamma(A)-log(C)-A/C)/(th^3) - a*trigamma(1/th)/(th^4)
  
  th1 = th - Q1/Q2
  if (th1 > 0) 
    th <- th1
  
  # Updating Coefficients 
  for (k in 1:p) {
    temp1 = sum(d*X[,, k]) - sum(d*Z[,, k]/W)
    temp2 = -sum(d*U[,, k]/W)
    coef[k, 1] = coef[k, 1] - temp1 / temp2
  }
 
  return(list(coef = coef, lambda = lambda, th = th))
}