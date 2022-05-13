

MEGammaFrailty <- function(y, X, d, coef, latent, theta) {
  
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
  
  
  ell = rep(0, 10000)
  iter = 1
  ell[iter] = GammaLik(y, X, d, coef, latent, theta)
  error = 3
  
  while (error > 0.000001) {
    
    temp1 = Facc(y, X, d, coef, latent, theta)
    
    coef1 = temp1$coef
    theta1 = temp1$theta
    u0_be = coef1 - coef
    
    temp2 = Facc(y, X, d, coef1, latent, theta1)
    
    coef2 = temp2$coef
    theta2 = temp2$theta
    v0_be = coef2 - 2*coef1 +  coef
    
    al0_be = sum(u0_be*v0_be) / sum(v0_be^2)
    if (al0_be > -1) {
      al0_be = -1
    }
    
    # Update Coefficients
    coef = coef - 2*al0_be*u0_be + al0_be^2*v0_be
    theta = theta2
    
    Ypre = matrix(0, a, b)
    La = matrix(0, a, b)
    
    for (i in 1:a) {
      for (j in 1:b) {
        La[i, j] = sum((y <= y[i, j]) * latent)
      }
    }
    
    for (k in 1:b) {
      Ypre[, k] = X[, k, ] %*% coef
    }
    
    D = rowSums(d)
    A = 1/theta + D
    C = 1/theta + rowSums(La * exp(Ypre))
    
    # Update Latent Variables
    for (i in 1:a) {
      for (j in 1:b) { 
        w = rowSums((y >= y[i, j]) * exp(Ypre)) 
        latent[i, j] = d[i, j] / sum(A*w/C)     
      }
    }
    
    mLambda = sum((y <= 1) * latent) 
    
    ell[iter + 1] = GammaLik(y, X, d, coef, latent, theta)
    error = abs(ell[iter + 1] - ell[iter]) / (1 + abs(ell[iter]))
    iter = iter + 1
    cat(error, '\n')
  }
  
  return(list(coef = coef, theta = theta, likelihood = ell[iter], iteration = iter, mLambda = mLambda))
}