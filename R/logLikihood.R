
logLikihood <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN") {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  La = (cumsum(lambda[order(vy)]))[rank(vy)]
  La = matrix(La, a, b)
  BE = array(rep(coef, each = a*b),c(a, b, length(coef)))
  A = rowSums(La*exp(apply(X*BE, c(1,2), sum)))
  D = rowSums(d)
  
  l1 = sum(log(lambda[vd != 0])) 
  l2 = sum(d*(apply(X*BE, c(1,2), sum))) 
  itao <- vector("numeric", length = a)
  for (i in 1:a) {  
    itao[i] = integrate(int_tao, lower = 0, upper = Inf, i = i, est.tht = est.tht, A = A, D = D, frailty = frailty, mode = 0)$value
  }
  l3 = sum(log(itao[]))
  
  if (frailty == "LogN") {
    l0 = -a*log(est.tht)/2 + l1 + l2 + l3
  }
  if (frailty == "InvGauss") {
    l0 = a*log(est.tht) + a*est.tht + l1 + l2 + l3
  }
  
  return(l0)
}