
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
  B = apply((lambda*exp(apply(X*BE, c(1,2), sum)))^d, 1, prod)
  D = rowSums(d)

  int0 <- vector("numeric", length = a)
  for (i in 1:a) {  
    int0[i] = integrate(int_tao, lower = 0, upper = Inf,
                        i = i, est.tht = est.tht, A = A, B = B, D = D, frailty = frailty, mode = 0)$value
  }

  return(sum(log(int0)))
}