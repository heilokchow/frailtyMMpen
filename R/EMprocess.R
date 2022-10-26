EMprocess <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN") {
  
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
  D = rowSums(d) 
  
  itao <- vector("numeric", length = a)
  itao1 <- vector("numeric", length = a)
  itao2 <- vector("numeric", length = a)
  
  for (i in 1:a) {  
    itao[i]  = integrate(int_tao, lower = 0, upper = Inf, i = i, est.tht = est.tht, A = A, D = D, frailty = frailty, mode = 0)$value
    itao1[i] = integrate(int_tao, lower = 0, upper = Inf, i = i, est.tht = est.tht, A = A, D = D, frailty = frailty, mode = 1)$value
    itao2[i] = integrate(int_tao, lower = 0, upper = Inf, i = i, est.tht = est.tht, A = A, D = D, frailty = frailty, mode = 2)$value
  }
  
  rat1 = itao1/itao
  rat2 = itao2/itao    
  ME = matrix(rat2, a, b)
  E_0 = as.vector(ME*YpreExp)
  SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
  SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
  
  # Update frailty parameters
  if (frailty == "LogN") {
    est.tht = sum(rat1)/a
  }
  
  if (frailty == "InvGauss") {
    DETA = a^2 + 4*a*sum(rat1)
    est.tht = (a + sqrt(DETA))/(2*sum(rat1))
  }
  
  # Update lambda Variables
  lambda = vd/SUM_0
  
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
  
  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}