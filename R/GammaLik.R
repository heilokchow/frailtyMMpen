

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



GammaLik2E = function(y, X, d, coef, la1, la2, th)
{
  n = dim(X)[2]
  q = dim(X)[3]
  
  coef = matrix(coef)
  
  La1 = La2 = rep(0,n)
  for(i in 1:n){
    La1[i] = sum(la1*(y[1,] <= y[1,i]))
    La2[i] = sum(la2*(y[2,] <= y[2,i]))
  }
  
  C = 1/th + La1*exp(X[1,,] %*% coef) + La2*exp(X[2,,] %*% coef)
  
  A = 1/th+colSums(d)
  AC = A/C
  
  d1 = d[1,]
  d2 = d[2,]
  
  l1 = sum(lgamma(A)) - n*(lgamma(1/th) + log(th)/th) - sum(A*log(C)) 
  l2 = sum(log(la1[d1 != 0])) + sum(log(la2[d2 != 0])) 
  l3 = sum(d[1,]*(X[1,,] %*% coef) + d[2,]*(X[2,,] %*% coef))
  
  return(l1 + l2 + l3)
  
}