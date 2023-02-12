frailtyMM_CL <- function(y, X, d, coef.ini = NULL, est.tht.ini = NULL, lambda.ini = NULL, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL, maxit = 200, threshold = 1e-5) {
  
  p = dim(X)[3]
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  if (frailty == "Gamma") {
    frailtyc = 0
  }
  
  if (frailty == "LogN") {
    frailtyc = 1
  }
  
  if (frailty == "InvGauss") {
    frailtyc = 2
  }
  
  if (frailty == "PVF") {
    frailtyc = 3
  }
  
  if (is.null(penalty)) {
    penaltyc = 0
    tune = 0
  }
  
  if (!is.null(penalty)) {
    
    if (penalty == "LASSO") {
      penaltyc = 1
    }
    
    if (penalty == "MCP") {
      penaltyc = 2
    }
    
    if (penalty == "SCAD") {
      penaltyc = 3
    }
    
  }
  
  if (is.null(power)) {
    power = 0.0
  }
  
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
    
    # coef = rep(1/p, p)
    # est.tht = 1
    # lambda = rep(1/N/10, N)

    if (!is.null(coef.ini)) {
      coef = coef.ini
    }
    
    if (!is.null(est.tht.ini)) {
      est.tht = est.tht.ini
    }
    
    if (!is.null(lambda.ini)) {
      lambda = lambda.ini
    } 
    
  }
  
  coef.ini = coef
  lambda.ini = lambda
  est.tht.ini = est.tht
  
  l0 = LogLikCL(y, X, d, coef, lambda, est.tht, frailtyc, a, b, p, power)
  error = 3

  ## MM iteration
  num = 0
  SQS1 = 1
  p0 = 0
  threshold = threshold * (p+1)
  while(error > threshold && num < maxit) {
    
    coef0 = coef
    est.tht0 = est.tht
    lambda0 = lambda

    p1 = proc.time()[1]
    # rs1 = MMprocess_CL(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    rs1 = MMCL_TEST(y, X, d, coef, lambda, est.tht, frailtyc, penaltyc, tune, a, b, p, power)
    p2 = proc.time()[1]
    
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda = rs1$lambda
    
    if (sum(abs(coef1)) < threshold) {
      break
    }

    u_be = coef1 - coef

    p3 = proc.time()[1]
    # rs2 = MMprocess_CL(y, X, d, coef1, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    rs2 = MMCL_TEST(y, X, d, coef1, lambda, est.tht, frailtyc, penaltyc, tune, a, b, p, power)
    p4 = proc.time()[1]
    
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda = rs2$lambda

    if (!SQS1) {
      coef = coef2
    }
  
    if (SQS1) {
      
      v_be = coef2 - 2*coef1 + coef
      al_be = sum(u_be*v_be)/sum(v_be^2)
      if (al_be > -1) {al_be = -1}
      
      for (k1 in 0:7) {
        dc = (2*al_be*u_be - al_be^2*v_be) * 2^(-k1)
  
        if (sum(abs(dc)) > p) {
          next
        }
  
        coef.temp = coef - dc
  
        p5 = proc.time()[1]
        # rs = MMprocess_CL(y, X, d, coef.temp, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
        rs = MMCL_TEST(y, X, d, coef.temp, lambda, est.tht, frailtyc, penaltyc, tune, a, b, p, power)
        p6 = proc.time()[1]
  
        if (rs$error) {
          coef = coef.ini
          est.tht = est.tht.ini
          lambda = lambda.ini
          SQS1 = 0
          num = 0
          break
        }
        
        if (!backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda)) {
          break
        }
      }
      
      if (!rs$error) {
        coef = rs$coef
        est.tht = rs$est.tht
        lambda = rs$lambda
      }
    }
 
    error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0))
    
    num = num + 1
    cat(error, " ", est.tht, " ", " ", al_be, " ", num, '\n')
    p0 = p0 + (p2 - p1) + (p4 - p3) + (p6 - p5)
  }
  
  l1 = LogLikCL(y, X, d, coef, lambda, est.tht, frailtyc, a, b, p, power)
  
  output = list(coef = coef,
                est.tht = est.tht,
                lambda = lambda,
                likelihood = l1,
                input = list(y = y, X = X, d = d),
                iter = num,
                convergence = error)
  
  cat("->", p0)
  return(output)
}


frailtyMM_ME <- function(y, X, d, coef.ini = NULL, est.tht.ini = NULL, lambda1.ini = NULL, lambda2.ini = NULL, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL, maxit = 200, threshold = 1e-6) {
  
  p = dim(X)[2]
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
  
  l0 = logLikihood_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty, power = power)
  l1 = l0
  error = 3
  
  num = 0
  # MM iteration
  while(error > threshold && num < maxit) {
    
    coef0 = coef
    est.tht0 = est.tht
    lambda01 = lambda1
    lambda02 = lambda2
    
    rs1 = MMprocess_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda1 = rs1$lambda1
    lambda2 = rs1$lambda2
    
    if (sum(abs(coef1)) < threshold) {
      break
    }
    
    u_be = coef1 - coef
    
    rs2 = MMprocess_ME(y, X, d, coef1, lambda1, lambda2, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
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
      rs = MMprocess_ME(y, X, d, coef.temp, lambda1, lambda2, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune) 
      if (!backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda1, lambda2 = lambda2)) {
        break
      } 
    }
    
    coef = rs$coef
    est.tht = rs$est.tht
    lambda1 = rs$lambda1
    lambda2 = rs$lambda2

    l1 = logLikihood_ME(y, X, d, coef, lambda1, lambda2, est.tht, frailty = frailty, power = power)
    
    if (is.null(penalty)) {
      error = abs(l1 - l0)/(1 + abs(l0))
      l0 = l1
    } else {
      error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0)) + sum(abs(lambda1 - lambda01)) + sum(abs(lambda2 - lambda02))
    }
    
    num = num + 1
    cat(error, " ", est.tht, " ", k1, " ", al_be, " ", num, '\n')
  }
  
  output = list(coef = coef,
                est.tht = est.tht,
                lambda1 = lambda1,
                lambda2 = lambda2,
                likelihood = l1,
                input = list(y = y, X = X, d = d),
                iter = num,
                convergence = error)
  
  return(output)
}




frailtyMM_RE <- function(y, X, d, coef.ini = NULL, est.tht.ini = NULL, lambda.ini = NULL, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL, maxit = 200, threshold = 1e-6) {
  
  p = dim(X[[1]])[2]
  n = length(y)
  
  # Initialize Parameters
  coef = rep(0.5, p)
  est.tht = 1
  lambda = lapply(y, function(x) {x^0*1/n})
  
  if (!is.null(lambda.ini)) {
    lambda = lambda.ini
  }
  
  if (!is.null(coef.ini)) {
    coef = coef.ini
  }
  
  if (!is.null(est.tht.ini)) {
    est.tht = est.tht.ini
  }
  
  l0 = logLikihood_RE(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power)
  l1 = l0
  error = 3
  
  num = 0
  # MM iteration
  while(error > threshold && num < maxit) {
    
    coef0 = coef
    est.tht0 = est.tht
    lambda0 = lambda
    
    rs1 = MMprocess_RE(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda = rs1$lambda

    if (sum(abs(coef1)) < threshold) {
      break
    }
    
    u_be = coef1 - coef
    
    rs2 = MMprocess_RE(y, X, d, coef1, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda = rs1$lambda
    
    v_be = coef2 - 2*coef1 + coef
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if (al_be > -1) {al_be = -1}
    
    for (k1 in 0:4) {
      dc = (2*al_be*u_be - al_be^2*v_be) * 2^(-k1)
      coef.temp = coef - dc
      rs = MMprocess_RE(y, X, d, coef.temp, lambda, est.tht, frailty = frailty, power = power, penalty = penalty, tune = tune)
      if (!backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda)) {
        break
      }
    }

    coef = rs$coef
    est.tht = rs$est.tht
    lambda = rs$lambda
    
    l1 = logLikihood_RE(y, X, d, coef, lambda, est.tht, frailty = frailty, power = power)
    
    if (is.null(penalty)) {
      error = abs(l1 - l0)/(1 + abs(l0))
      l0 = l1
    } else {
      error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0)) + sum(abs(unlist(lambda) - unlist(lambda0)))
    }
    
    num = num + 1
    cat(error, " ", est.tht, " ", l1, " ", num, '\n')
  }
  
  output = list(coef = coef,
                est.tht = est.tht,
                lambda = lambda,
                likelihood = l1,
                input = list(y = y, X = X, d = d),
                iter = num,
                convergence = error)
  
  return(output)
}
