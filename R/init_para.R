init_para <- function(y, X, d, frailty = "LogN") {
  
  p = dim(X)[3]
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  # Initialize Parameters
  newX = matrix(0, nrow = 0, ncol = p)
  for (i in 1:a) {
    newX = rbind(newX, X[i,,])
  }
  newY <- Surv(rep(0, N), as.vector(t(y)), as.vector(t(d)))
  
  coxret = survival::agreg.fit(x = newX, y = newY, strata = NULL, offset = NULL, init = NULL, control = survival::coxph.control(),
                               weights = NULL, method = "breslow", rownames = NULL)  
  
  coef = coxret$coefficients
  
  if (frailty == "Gamma") {
    est.tht = 1
    lambda = matrix(rep(1/N/10, N), a, b)
    return(list(coef = coef, est.tht = est.tht, lambda = lambda))
  }
  
  YpreExp = matrix(0, a, b)
  for (i in 1:a) {
    YpreExp[i,] = exp(X[i,,] %*% coef)
  }
  
  lambda = rep(1/N/10, N)
  est.tht = 1
  
  La = (cumsum(lambda[order(vy)]))[rank(vy)]
  La = matrix(La, a, b)
  
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
  
  if (frailty == "LogN" || frailty == "InvGauss") {
    est.tht = sum(int2)/a
  }
  
  if (est.tht < 0) {
    est.tht = 1
  }
  
  ME = matrix(int1, a, b)
  E_0 = as.vector(ME*YpreExp)
  SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
  SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
  
  lambda = vd/SUM_0
  
  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}