

sampleGamma <- function(coef = matrix(c(-2,3,-4,5,-6,7,6,5,4,3)), lambda = 5, th = 10, a = 100, b = 5, cen = 2, seed = 0) {
  
  p = length(coef)
  
  set.seed(seed)
  u = rgamma(a, 1/th, scale = 1/th)
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

sampleGamma2E <- function(coef = matrix(c(-2,3,-4,5,-6,7,6,5,4,3)), th = 10, a1 = 3, a2 = 5, n = 100, cen1 = 3, seed = 0) {
  
  q = length(coef)
  
  set.seed(seed)
  u = rgamma(a, 1/th, scale = 1/th)
  X = array(runif(2*n*q, min = 0, max = 0.5), dim=c(2,n,q))
  
  T = matrix(0,2,n)  
  cen = runif(n,0,cen1)
  
  U1 = runif(n,0,1)    
  U2 = runif(n,0,1)
  
  
  T[1,] <- -log(U1)/(a1*u*exp(apply(x[1,,]*Be,1,sum))) 
  T[2,] <- (exp(-log(U2)/(u*exp(apply(x[2,,]*Be,1,sum))))-1)/a2
  
  d = 1*(T <= cen)
  y = pmin(T, cen)
  
  la1 = a1*y[1,]
  la2 = log(1 + a2*y[2,])
  return(list(y = y, X = X, d = d, la1 = la1, la2 = la2))
}