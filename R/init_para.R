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
  lambda = rep(1/N/10, N)
  est.tht = 1
  
  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}