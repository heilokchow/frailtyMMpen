MMroutine <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type) {
  
  TEST = NULL
  
  if (type == 1) {
    TEST = MMCL(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power)
  }
  
  if (type == 2) {
    TEST = MMME(y, X, d, coef, lambda, tht, frailty, penalty, tune, N, a, p, power)
  }
  
  if (type == 3) {
    TEST = MMRE(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power)
  }
  
  return(TEST)
}

MMRE_TEST <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type) {
  
  TEST = NULL
  
  if (type == 1) {
    TEST = MMRE(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power)
  }
  
  return(TEST)
}

logLikcal <- function(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type) {
  
  s = 0
  
  if (type == 1) {
    s = LogLikCL(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power)
  }
  
  if (type == 2) {
    s = LogLikME(y, X, d, coef, lambda, est.tht, frailtyc, N, a, p, power)
  }
  
  if (type == 3) {
    s = LogLikRE(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power)
  }
  
  return(s)
  
}

logLik <- function(x1, data, lambda, frailtyc, id, N, a, p, power, type) {
  
  est.tht = x1[1]
  coef = x1[-1]
  y = data$y
  X = data$X
  d = data$d
  
  return(logLikcal(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type))
  
}

#'@export
cluster <- function(x) {
  x
}

#'@export
event <- function(x) {
  x
}