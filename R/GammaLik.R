

GammaLik <- function(y, X, d, coef, lambda, th) {

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
  D = rowSums(d)
  
  for (k in 1:b) {
    Ypre[, k] = X[, k, ] %*% coef
  }
  
  l3 = 0
  for(i in 1:a) {
    for(j in 1:b) {
      La[i, j] = sum((y <= y[i, j]) * lambda)
      if(lambda[i, j] > 0) {
        l3 = l3 + log(lambda[i,j])
      }
    }
  }

  A = rowSums(La * exp(Ypre))
  l1 = sum(lgamma(D + 1/th)) - a*(lgamma(1/th) + log(th)/th)
  l2 = sum(d*Ypre) - sum((1/th + D)*log(1/th + A))

  return(l1 + l2 + l3)
  
}