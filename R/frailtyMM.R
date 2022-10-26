frailtyMM <- function(y, X, d, frailty = "LogN") {
  
  p = dim(X)[3]
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  coef = rep(0.5, p)
  est.tht = 1
  lambda = rep(1/N, N)
  
  ell = rep(0,1000000)
  k = 1
  
  ell[k]= logLikihood(y, X, d, coef, lambda, est.tht, frailty = "LogN")
  error = 3

  while(error > 0.000001) {
    
    rs1 = EMprocess(y, X, d, coef, lambda, est.tht, frailty = "LogN") 
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda = rs1$lambda
    
    u_be = coef1 - coef
    
    rs2 = EMprocess(y, X, d, coef1, lambda, est.tht, frailty = "LogN") 
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda = rs2$lambda
    
    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}
    
    coef = coef - 2*al_be*u_be + al_be^2*v_be
    
    rs = EMprocess(y, X, d, coef, lambda, est.tht, frailty = "LogN") 
    coef = rs$coef
    est.tht = rs$est.tht
    lambda = rs$lambda
    
    ell[k+1] = logLikihood(y, X, d, coef, lambda, est.tht, frailty = "LogN")
    
    error = abs(ell[k+1] - ell[k])/(1 + abs(ell[k]))
    cat(error, '\n')
    k = k + 1
  }

  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}