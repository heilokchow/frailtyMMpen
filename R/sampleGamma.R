

sampleGamma <- function(coef = matrix(c(-2,3,-4,5,-6,7,6,5,4,3)), lambda = 5, alp = 0.1, th = 0.1, a = 100, b = 5, cen = 2, seed = 0) {
  
  p = length(coef)
  
  set.seed(seed)
  u = rgamma(a, alp, scale = th)
  X = array(runif(a*b*p, min = 0, max = 0.5), c(a,b,p))
  U = matrix(runif(a*b,0,1),a,b)
  
  T = matrix(0,a,b)  
  
  for (k in 1:b) {
    T[, k] = -log(U[,k])/(lambda * u * exp(X[,k,] %*% coef))
  }
  
  d = 1*(T <= cen)
  y = pmin(T, cen)
  return(list(y = y, X = X, d = d))
}