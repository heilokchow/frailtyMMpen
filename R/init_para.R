init_para <- function(y, X, d, frailty = "LogN") {
  
  p = dim(X)[2]
  N = length(y)

  newY <- Surv(y, d)
  
  coxret = survival::coxph.fit(x = X, y = newY, strata = NULL, offset = NULL, init = NULL, control = survival::coxph.control(),
                               weights = NULL, method = "breslow", rownames = NULL)  
  
  coef = coxret$coefficients
  lambda = rep(1/N/10, N)
  est.tht = 1
  
  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}