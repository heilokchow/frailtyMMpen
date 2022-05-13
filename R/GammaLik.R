

GammaLik <- function(y, X, d, coef, latent, theta) {

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
  
  La = matrix(0, a, b)
  Ypre = matrix(0, a, b)
  D = rowSums(d)
  
  for (k in 1:b) {
    Ypre[, k] = X[, k, ] %*% coef
  }
  
  l3 = 0
  for(i in 1:a) {
    for(j in 1:b) {
      La[i, j] = sum((y <= y[i, j]) * latent)
      if(latent[i, j] > 0) {
        l3 = l3 + log(latent[i,j])
      }
    }
  }

  A = rowSums(La * exp(Ypre))
  l1 = sum(lgamma(D + 1/theta)) - a*(lgamma(1/theta) + log(theta)/theta)
  l2 = sum(d*Ypre) - sum((1/theta + D)*log(1/theta + A))

  return(l1 + l2 + l3)
  
}