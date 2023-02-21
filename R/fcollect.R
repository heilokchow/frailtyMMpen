#'@export
MMroutine <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type) {
  
  TEST = NULL
  
  if (type == 1) {
    TEST = MMCL(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power)
  }
  
  if (type == 2) {
    TEST = MMME(y, X, d, coef, lambda, tht, frailty, penalty, tune, N, a, p, power)
  }
  
  if (type == 3) {
    TEST = MMRE(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type)
  }
  
  return(TEST)
}

#'@export
MMRE_TEST <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type) {
  
  TEST = NULL
  
  if (type == 1) {
    TEST = MMRE(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power)
  }
  
  return(TEST)
}

#'@export
logLikcal <- function(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type) {
  
  s = 0
  
  if (type == 1) {
    s = LogLikCL(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power)
  }
  
  if (type == 2) {
    s = LogLikME(y, X, d, coef, lambda, est.tht, frailtyc, N, a, p, power)
  }
  
  return(s)
  
}

#'@export
cluster <- function(x) {
  x
}

#'@export
event <- function(x) {
  x
}