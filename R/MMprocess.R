MMprocess <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN") {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  La = (cumsum(lambda[order(vy)]))[rank(vy)]
  La = matrix(La, a, b)
  
  YpreExp = matrix(0, a, b)
  for (i in 1:a) {
    YpreExp[i,] = exp(X[i,,] %*% coef)
  }
  
  A = rowSums(La*YpreExp)
  B = apply((lambda*YpreExp)^d, 1, prod)
  D = rowSums(d) 
  
  int0 <- vector("numeric", length = a)
  int1 <- vector("numeric", length = a)
  int2 <- vector("numeric", length = a)
  
  for (i in 1:a) {  
    int0[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = A, B = B, D = D, frailty = frailty, mode = 0)$value
  }
  
  for (i in 1:a) {  
    int1[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, mode = 1)$value
    int2[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, mode = 2)$value
  }
  
  ME = matrix(int1, a, b)
  E_0 = as.vector(ME*YpreExp)
  SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
  SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
  
  # Update coefficients
  AVE_X = apply(abs(X), c(1,2), sum)
  for(k in 1:p) {
    E_1 = as.vector(ME*X[,,k]*YpreExp)
    E_2 = as.vector(ME*AVE_X*abs(X[,,k])*YpreExp)
    
    SUM_1 = cumsum((E_1[order(vy)])[seq(N, 1, -1)])
    SUM_1 = (SUM_1[seq(N, 1, -1)])[rank(vy)] 
    SUM_2 = cumsum((E_2[order(vy)])[seq(N, 1, -1)])
    SUM_2 = (SUM_2[seq(N, 1, -1)])[rank(vy)] 
    
    DE_1 = sum(d*X[,,k]) - sum(vd*SUM_1/SUM_0) 
    DE_2 = -sum(vd*SUM_2/SUM_0) 
    coef[k] = coef[k] - DE_1/DE_2
  }
  
  # Update frailty parameters
  if (frailty == "LogN") {
    est.tht = sum(int2)/a
  }
  
  if (frailty == "InvGauss") {
    q1 = a/est.tht + a - 2*est.tht*sum(int2)
    q2 = -a/(est.tht^2) - 2*sum(int2)
    est.tht = est.tht - q1/q2
  }
  
  if (est.tht < 0) {
    est.tht = 1
  }
  
  # Update lambda Variables (Use Updated coef and est.tht)
  
  YpreExp = matrix(0, a, b)
  for (i in 1:a) {
    YpreExp[i,] = exp(X[i,,] %*% coef)
  }
  
  A = rowSums(La*YpreExp)
  B = apply((lambda*YpreExp)^d, 1, prod)
  D = rowSums(d) 
  
  int0 <- vector("numeric", length = a)
  int1 <- vector("numeric", length = a)
  
  for (i in 1:a) {  
    int0[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = A, B = B, D = D, frailty = frailty, mode = 0)$value
  }
  
  for (i in 1:a) {  
    int1[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, mode = 1)$value
  }
  
  ME = matrix(int1, a, b)
  E_0 = as.vector(ME*YpreExp)
  SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
  SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
  
  lambda = vd/SUM_0
  
  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}