
logLikihood_CL <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN", power = NULL) {
  
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
                        i = i, est.tht = est.tht, A = A, B = B, D = D, frailty = frailty, power = power, mode = 0)$value
  }

  return(sum(log(int0)))
}


logLikihood_ME <- function(y, X, d, coef, lambda1, lambda2, est.tht, frailty = "LogN") {
  
  p = length(coef)
  coef = as.matrix(coef)
  n = ncol(y)
  
  La1 = La2 = rep(0,n)
  for(i in 1:n){
    La1[i] = sum(lambda1*(y[1,] <= y[1,i]))
    La2[i] = sum(lambda2*(y[2,] <= y[2,i]))
  }
  
  YpreExp = matrix(0, 2, n)
  for (i in 1:2) {
    YpreExp[i,] = exp(X[i,,] %*% coef)
  }
  
  CC =  La1*YpreExp[1,] + La2*YpreExp[2,]
  D = colSums(d)
  d1 = d[1,]
  d2 = d[2,]
  
  AA = (lambda1*YpreExp[1,])^(d1)*(lambda2*YpreExp[2,])^(d2)
  
  if (frailty == "Gamma") {
    C = 1/est.tht + CC
    A = 1/est.tht + D
    AC = A/C
    
    l1 = sum(lgamma(A)) - n*(lgamma(1/est.tht)+log(est.tht)/est.tht) - sum(A*log(C)) 
    l2 = sum(log(lambda1[d1 != 0])) + sum(log(lambda2[d2 != 0])) 
    l3 = sum(d1*rowSums(X[1,,] %*% coef) + d2*rowSums(X[2,,] %*% coef))
    
    return(l1+l2+l3) 
  }
  
  int0 <- vector("numeric", length = n)
  
  for (i in 1:n) {  
    int0[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = CC, B = AA, D = D, frailty = frailty, mode = 0)$value
  }
  
  return(sum(log(int0)))
}