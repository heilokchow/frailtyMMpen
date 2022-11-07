frailtyMM <- function(y, X, d, frailty = "LogN") {
  
  p = dim(X)[3]
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  coef = rep(0.5, p)
  est.tht = 1
  lambda = rep(1/N/10, N)
  
  ell = rep(0,1000000)
  k = 1
  
  ell[k]= logLikihood(y, X, d, coef, lambda, est.tht, frailty = frailty)
  error = 3

  while(error > 0.000001) {

    rs1 = MMprocess(y, X, d, coef, lambda, est.tht, frailty = frailty)
    coef1 = rs1$coef
    est.tht1 = rs1$est.tht
    lambda = rs1$lambda

    u_be = coef1 - coef

    rs2 = MMprocess(y, X, d, coef1, lambda, est.tht1, frailty = frailty)
    coef2 = rs2$coef
    est.tht2 = rs2$est.tht
    lambda = rs2$lambda

    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}

    coef = coef - 2*al_be*u_be + al_be^2*v_be
    est.tht = est.tht2

    rs = MMprocess(y, X, d, coef, lambda, est.tht, frailty = frailty) 
    coef = rs$coef
    est.tht = rs$est.tht
    lambda = rs$lambda
    
    ell[k+1] = logLikihood(y, X, d, coef, lambda, est.tht, frailty = frailty)
    
    error = abs(ell[k+1] - ell[k])/(1 + abs(ell[k]))
    cat(error, " ", est.tht, '\n')
    k = k + 1
  }

  return(list(coef = coef, est.tht = est.tht, lambda = lambda, likelihood = ell[k]))
}