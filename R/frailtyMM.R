frailtyMM_CL <- function(y, X, d, frailty = "LogN") {
  
  p = dim(X)[3]
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  # Initialize Parameters
  init = init_para(y, X, d, frailty = frailty)
  coef = init$coef
  est.tht = init$est.tht
  lambda = init$lambda

  if (frailty == "Gamma") {
    return(CLGammaFrailty(y, X, d, coef, lambda, est.tht))
  }
  
  ell = rep(0,1000000)
  k = 1
  
  ell[k]= logLikihood_CL(y, X, d, coef, lambda, est.tht, frailty = frailty)
  error = 3

  while(error > 0.000001) {

    rs1 = MMprocess_CL(y, X, d, coef, lambda, est.tht, frailty = frailty)
    coef1 = rs1$coef
    est.tht1 = rs1$est.tht
    lambda = rs1$lambda

    u_be = coef1 - coef

    rs2 = MMprocess_CL(y, X, d, coef1, lambda, est.tht1, frailty = frailty)
    coef2 = rs2$coef
    est.tht2 = rs2$est.tht
    lambda = rs2$lambda

    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}

    coef = coef - 2*al_be*u_be + al_be^2*v_be
    est.tht = est.tht2

    rs = MMprocess_CL(y, X, d, coef, lambda, est.tht, frailty = frailty) 
    coef = rs$coef
    est.tht = rs$est.tht
    lambda = rs$lambda
    
    ell[k+1] = logLikihood_CL(y, X, d, coef, lambda, est.tht, frailty = frailty)
    
    error = abs(ell[k+1] - ell[k])/(1 + abs(ell[k]))
    cat(error, " ", est.tht, '\n')
    k = k + 1
  }

  return(list(coef = coef, est.tht = est.tht, lambda = lambda, likelihood = ell[k]))
}


frailtyMM_ME <- function(y, X, d, frailty = "LogN") {
  
  p = dim(X)[3]
  n = ncol(y)
  
  # Initialize Parameters
  coef = rep(0.5, p)
  est.tht = 1
  lambda1 = rep(1/n, n)
  lambda2 = rep(1/n, n)
  
  ell = rep(0,1000000)
  k = 1
  
  ell[k]= logLikihood_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty)
  error = 3
  
  while(error > 0.000001) {
    
    rs1 = MMprocess_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty)
    coef1 = rs1$coef
    est.tht1 = rs1$est.tht
    lambda1 = rs1$lambda1
    lambda2 = rs1$lambda2
    
    u_be = coef1 - coef
    
    rs2 = MMprocess_ME(y, X, d, coef1, lambda1, lambda2, est.tht1, frailty = frailty)
    coef2 = rs2$coef
    est.tht2 = rs2$est.tht
    lambda1 = rs2$lambda1
    lambda2 = rs2$lambda2
    
    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}
    
    coef = coef - 2*al_be*u_be + al_be^2*v_be
    est.tht = est.tht2
    
    rs = MMprocess_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty) 
    coef = rs$coef
    est.tht = rs$est.tht
    lambda1 = rs$lambda1
    lambda2 = rs$lambda2
    
    ell[k+1] = logLikihood_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty)
    
    error = abs(ell[k+1] - ell[k])/(1 + abs(ell[k]))
    cat(error, " ", est.tht, '\n')
    k = k + 1
  }
  
  return(list(coef = coef, est.tht = est.tht, lambda1 = lambda1, lambda2 = lambda2, likelihood = ell[k]))
}
