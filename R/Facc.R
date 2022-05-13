

Facc <- function(y, X, d, coef, latent, theta) {
  
  a = dim(X)[1]
  b = dim(X)[2]
  p = dim(X)[3]
  
  coef = matrix(coef)
  
  if (dim(latent)[1] != a || dim(latent)[2] != b) {
    stop("Dimension of y, X, d or latent is incorrect.")
  }
  
  if (length(coef) != p) {
    stop("Dimension of coef is incorrect.")
  }
  
  coef = matrix(coef)
  
  La = matrix(0, a, b)
  Ypre = matrix(0, a, b)
  Xabs = matrix(0, a, b)
  D = rowSums(d)
  
  Z = array(0, c(a, b, p))
  U = array(0, c(a, b, p))
  W = matrix(0, a, b)
  
  for(i in 1:a) {
    for(j in 1:b) {
      La[i, j] = sum((y <= y[i, j]) * latent)   
    }
  }
  
  for (k in 1:b) {
    Ypre[, k] = X[, k, ] %*% coef
  }
  
  for (k in 1:p) {
    Xabs[, ] = Xabs[, ] + abs(X[,, k])
  }
  
  A = 1/theta + D
  C = 1/theta + rowSums(La * exp(Ypre))
  
  for(i in 1:a) {
    for(j in 1:b) { 
      for (k in 1:p) {
        temp1 = rowSums((y >= y[i, j]) * X[,, k] * exp(Ypre))
        temp2 = rowSums((y >= y[i, j]) * abs(X[,, k]) * exp(Ypre) * Xabs)
        Z[i, j, k] = sum(A*temp1/C)
        U[i, j, k] = sum(A*temp2/C)
      }
      temp = rowSums((y >= y[i, j]) * exp(Ypre))
      W[i, j] = sum(A*temp/C)
    }
  }
  
  # Update theta for Gamma Frailty 
  Q1 = a*(digamma(1/theta) + log(theta)-1)/(theta^2) + sum(A/C - digamma(A) + log(C))/(theta^2) 
  Q2 = a*(3 - 2*digamma(1/theta) - 2*log(theta))/(theta^3) + 2*sum(digamma(A)-log(C)-A/C)/(theta^3) - a*trigamma(1/theta)/(theta^4)
  
  theta1 = theta - Q1/Q2
  if (theta1 > 0) 
    theta <- theta1
  
  # Updating Coefficients 
  for (k in 1:p) {
    temp1 = sum(d*X[,, k]) - sum(d*Z[,, k]/W)
    temp2 = -sum(d*U[,, k]/W)
    coef[k, 1] = coef[k, 1] - temp1 / temp2
  }
 
  return(list(coef = coef, latent = latent, theta = theta))
}