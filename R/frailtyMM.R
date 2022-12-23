frailtyMM_CL <- function(y, X, d, coef.ini = NULL, est.tht.ini = NULL, lambda.ini = NULL, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL, maxit = 1000, threshold = 1e-6) {
  
  p = dim(X)[3]
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  # Initialize Parameters
  
  if (!is.null(lambda.ini) && !is.null(coef.ini) && !is.null(est.tht.ini)) {
    
    coef = coef.ini
    est.tht = est.tht.ini
    lambda = lambda.ini
    
  } else {
    
    init = init_para(y, X, d, frailty = frailty)
    coef = init$coef
    est.tht = init$est.tht
    lambda = init$lambda
    
    if (!is.null(lambda.ini)) {
      lambda = lambda.ini
    }
    
    if (!is.null(coef.ini)) {
      coef = coef.ini
    }
    
    if (!is.null(est.tht.ini)) {
      est.tht = est.tht.ini
    }
  }
  
  if (frailty == "Gamma") {
    return(CLGammaFrailty(y, X, d, coef, lambda, est.tht, penalty = penalty, tune = tune, maxit = maxit, threshold = threshold))
  }
  
  l0 = logLikihood_CL(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power)
  error = 3

  num = 0
  while(error > threshold && num < maxit) {
    
    coef0 = coef
    est.tht0 = est.tht
    lambda0 = lambda
    
    rs1 = MMprocess_CL(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda = rs1$lambda
    
    if (sum(abs(coef1)) < threshold) {
      break
    }

    u_be = coef1 - coef

    rs2 = MMprocess_CL(y, X, d, coef1, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda = rs2$lambda

    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}

    for (k1 in 0:4) {
      dc = (2*al_be*u_be - al_be^2*v_be) * 2^(-k1)
      coef.temp = coef - dc
      rs = MMprocess_CL(y, X, d, coef.temp, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune) 
      if (!backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda)) {
        break
      } 
    }
    
    coef = rs$coef
    est.tht = rs$est.tht
    lambda = rs$lambda
      
    l1 = logLikihood_CL(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power)
    
    if (is.null(penalty)) {
      error = abs(l1 - l0)/(1 + abs(l0))
      l0 = l1
    } else {
      error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0)) + sum(abs(lambda - lambda0))
    }
    
    num = num + 1
    cat(error, " ", est.tht, " ", k1, " ", al_be, " ", num, '\n')
  }

  return(list(coef = coef, est.tht = est.tht, lambda = lambda, likelihood = l1))
}


frailtyMM_ME <- function(y, X, d, coef.ini = NULL, est.tht.ini = NULL, lambda1.ini = NULL, lambda2.ini = NULL, frailty = "LogN", penalty = NULL, tune = NULL, maxit = 1000, threshold = 1e-6) {
  
  p = dim(X)[3]
  n = ncol(y)
  
  # Initialize Parameters
  coef = rep(0.5, p)
  est.tht = 1
  lambda1 = rep(1/n, n)
  lambda2 = rep(1/n, n)
  
  if (!is.null(lambda1.ini)) {
    lambda1 = lambda1.ini
  }
  
  if (!is.null(lambda2.ini)) {
    lambda2 = lambda2.ini
  }
  
  if (!is.null(coef.ini)) {
    coef = coef.ini
  }
  
  if (!is.null(est.tht.ini)) {
    est.tht = est.tht.ini
  }
  
  l0 = logLikihood_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty)
  error = 3
  
  num = 0
  while(error > threshold && num < maxit) {
    
    coef0 = coef
    est.tht0 = est.tht
    lambda01 = lambda1
    lambda02 = lambda2
    
    rs1 = MMprocess_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty, penalty = penalty, tune = tune)
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda1 = rs1$lambda1
    lambda2 = rs1$lambda2
    
    if (sum(abs(coef1)) < threshold) {
      break
    }
    
    u_be = coef1 - coef
    
    rs2 = MMprocess_ME(y, X, d, coef1, lambda1, lambda2, est.tht, frailty = frailty, penalty = penalty, tune = tune)
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda1 = rs2$lambda1
    lambda2 = rs2$lambda2
    
    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}
    
    for (k1 in 0:4) {
      dc = (2*al_be*u_be - al_be^2*v_be) * 2^(-k1)
      coef.temp = coef - dc
      rs = MMprocess_ME(y, X, d, coef.temp, lambda1, lambda2, est.tht, frailty = frailty, penalty = penalty, tune = tune) 
      if (!backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda1, lambda2 = lambda2)) {
        break
      } 
    }
    
    coef = rs$coef
    est.tht = rs$est.tht
    lambda1 = rs$lambda1
    lambda2 = rs$lambda2

    l1 = logLikihood_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty)
    
    if (is.null(penalty)) {
      error = abs(l1 - l0)/(1 + abs(l0))
      l0 = l1
    } else {
      error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0)) + sum(abs(lambda1 - lambda01)) + sum(abs(lambda2 - lambda02))
    }
    
    num = num + 1
    cat(error, " ", est.tht, " ", k1, " ", al_be, " ", num, '\n')
  }
  
  return(list(coef = coef, est.tht = est.tht, lambda1 = lambda1, lambda2 = lambda2, likelihood = l1))
}
