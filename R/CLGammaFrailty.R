

CLGammaFrailty <- function(y, X, d, coef, lambda, th, penalty = NULL, tune = NULL) {
  
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
  
  l0 = GammaLik(y, X, d, coef, lambda, th)
  error = 3
  num = 0
  
  while (error > 0.000001 && num < 1000) {
    
    coef0 = coef
    th0 = th
    lambda0 = lambda
    
    temp1 = FaccAc(y, X, d, coef, lambda, th, penalty, tune)
    
    coef1 = temp1$coef
    th1 = temp1$th
    u0_be = coef1 - coef
    
    temp2 = FaccAc(y, X, d, coef1, lambda, th1, penalty, tune)
    
    coef2 = temp2$coef
    th2 = temp2$th
    v0_be = coef2 - 2*coef1 +  coef
    
    al0_be = sum(u0_be*v0_be) / sum(v0_be^2)
    if (al0_be > -1) {
      al0_be = -1
    }
    
    # Update Coefficients
    coef = coef - 2*al0_be*u0_be + al0_be^2*v0_be
    th = th2
    
    Ypre = matrix(0, a, b)
    La = matrix(0, a, b)
    
    for (i in 1:a) {
      for (j in 1:b) {
        La[i, j] = sum((y <= y[i, j]) * lambda)
      }
    }
    
    for (k in 1:b) {
      Ypre[, k] = X[, k, ] %*% coef
    }
    
    D = rowSums(d)
    A = 1/th + D
    C = 1/th + rowSums(La * exp(Ypre))
    YpreExp = exp(Ypre)
    
    # Update lambda Variables
    lambda1 = matrix(0, a, b)
    lambda1 = FaccCal(y, YpreExp, lambda1, A, C)
    lambda = d / lambda1

    mLambda = sum((y <= 1) * lambda) 
    
    l1 = GammaLik(y, X, d, coef, lambda, th)
    if (is.null(penalty)) {
      error = abs(l1 - l0)/(1 + abs(l0))
      l0 = l1
    } else {
      error = sum(abs(coef - coef0)) + sum(abs(th - th0)) + sum(abs(lambda - lambda0))
    }
    
    num = num + 1
    cat(error, '\n')
  }
  
  return(list(coef = coef, est.tht = th, lambda = lambda, likelihood = l1))
}