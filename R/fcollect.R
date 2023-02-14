#'@export
MMroutine <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type) {
  
  TEST = NULL
  
  if (type == 1) {
    TEST = MMCL(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power)
  }
  
  return(TEST)
}

logLikcal <- function(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type) {
  
  s = 0
  
  if (type == 1) {
    s = LogLikCL(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power)
  }
  
  return(s)
  
}