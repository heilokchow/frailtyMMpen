#'@export
MMCL_TEST <- function(y, X, d, coef, lambda, tht, frailty, penalty, tune, a, b, p) {
  TEST = MMCL(y, X, d, coef, lambda, tht, frailty, penalty, tune, a, b, p)
  return(TEST)
}