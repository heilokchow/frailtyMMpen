backtrackerror <- function(model, coef, est.tht, lambda, lambda2 = NULL) {
  coefn = model$coef
  est.thtn = model$est.tht
  
  c1 = abs((sum(coefn) - sum(coef)) / sum(coef))
  c2 = abs((est.thtn - est.tht) / est.tht)
  
  if (is.null(lambda2)) {
    lambdan = model$lambda
    c3 = abs((sum(lambdan) - sum(lambda)) / sum(lambda))
  } else {
    lambdan1 = model$lambda1
    lambdan2 = model$lambda2
    c3 = max(abs((sum(lambdan1) - sum(lambda)) / sum(lambda)), abs((sum(lambdan2) - sum(lambda2)) / sum(lambda2)))
  }

  if (c1 > 10 || c2 > 5 || c3 > 10) {
    return(1)
  }
  return(0)
}