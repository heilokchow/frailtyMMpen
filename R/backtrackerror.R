backtrackerror <- function(model, coef, est.tht, lambda) {
  coefn = model$coef
  est.thtn = model$est.tht
  lambdan = model$lambda
  
  c1 = abs((sum(coefn) - sum(coef)) / sum(coef))
  c2 = abs((est.thtn - est.tht) / est.tht)
  c3 = abs((sum(lambdan) - sum(lambda)) / sum(lambda))

  if (c1 > 10 || c2 > 5 || c3 > 10) {
    return(1)
  }
  return(0)
}