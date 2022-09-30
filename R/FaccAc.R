

FaccAc <- function(y, X, d, coef, lambda, alp, th) {
  
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
  
  A = alp + D
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
  th =  sum(A/C) / (a*alp)
  Q1 = a*(-log(th) - digamma(alp)) + sum(digamma(A) - log(C)) 
  Q2 = -a*trigamma(alp)
  alp1 = alp - Q1/Q2
  
  if (alp1 > 0) 
    alp <- alp1
  
  # Updating Coefficients 
  for (k in 1:p) {
    temp1 = sum(d*X[,, k]) - sum(d*Z[,, k]/W)
    temp2 = -sum(d*U[,, k]/W)
    coef[k, 1] = coef[k, 1] - temp1 / temp2
  }
 
  return(list(coef = coef, lambda = lambda, alp = alp, th = th))
}