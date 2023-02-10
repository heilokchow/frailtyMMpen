#'@export
MMCL_TEST <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, a, b, p, power) {
  TEST = MMCL(y, X, d, coef, lambda, tht, frailty, penalty, tune, a, b, p, power)
  return(TEST)
}