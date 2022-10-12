

CLGammaFrailty <- function(y, X, d, coef, lambda, th) {
  
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
  
  
  ell = rep(0, 10000)
  iter = 1
  ell[iter] = GammaLik(y, X, d, coef, lambda, th)
  error = 3
  
  while (error > 0.000001) {
    
    temp1 = FaccAc(y, X, d, coef, lambda, th)
    
    coef1 = temp1$coef
    th1 = temp1$th
    u0_be = coef1 - coef
    
    temp2 = FaccAc(y, X, d, coef1, lambda, th1)
    
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
    
    ell[iter + 1] = GammaLik(y, X, d, coef, lambda, th)
  
    error = abs(ell[iter + 1] - ell[iter]) / (1 + abs(ell[iter]))
    iter = iter + 1
    cat(error, '\n')
  }
  
  return(list(coef = coef, th = th, likelihood = ell[iter], iteration = iter, mLambda = mLambda))
}