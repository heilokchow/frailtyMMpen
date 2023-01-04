
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


logLikihood_RE <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN") {
  
  p = length(coef)
  coef = as.matrix(coef)
  n = length(y)
  
  vy = unlist(y)
  yend = unlist(lapply(y, max))
  vd = unlist(d)
  vl = unlist(lambda)

  tall = sort(vy)
  ntime = length(tall)
  Xmat = array(0, c(n, p, ntime))
  
  for (i in 1:n) {
    for (j in 1:p) {
      cont = 1
      for (z in 1:ntime) {
        it = length(y[[i]])
        Xmat[i, j, z] = X[[i]][cont, j]
        if (cont < it && y[[i]][cont] < tall[z] && y[[i]][cont + 1] >= tall[z]) {
          cont = cont + 1
        }
      }
    }
  }
  
  La = lambda
  for (i in 1:n) {
    it = length(y[[i]])
    for (j in 1:it) {
      La[[i]][j] =  sum(vl*(vy <= y[[i]][j]))
    }
  }
  
  Yexp = lapply(X, function(mx, mcoef) {c(exp(mx %*% mcoef))}, mcoef = coef)
  Ycoef = lapply(X, function(mx, mcoef) {c((mx %*% mcoef))}, mcoef = coef)
  
  AA = rep(0, n)
  for (i in 1:n) {
    it = length(y[[i]])
    L0 = 0
    for (j in 1:it) {
      AA[i] = AA[i] + (La[[i]][j] - L0) * Yexp[[i]][j]
      L0 = La[[i]][j]
    }
  }
  
  DD = unlist(lapply(d, sum))
  
  int0 <- vector("numeric", length = n)
  int1 <- vector("numeric", length = n)
  
  if (frailty == "Gamma") {
    C = 1/est.tht + AA
    A = 1/est.tht + DD
    AC = A/C
    
    l1 = sum(lgamma(A)) - n*(lgamma(1/est.tht)+log(est.tht)/est.tht) - sum(A*log(C)) 
    lambdaall = unlist(lambda)
    l2 = sum(log(lambdaall[lambdaall > 0]))
    l3 = sum(unlist(purrr::pmap(list(d, Ycoef), function(x, y) {x*y})))
    
    return(l1+l2+l3) 
  }
  
  int0 <- vector("numeric", length = n)
  BB = unlist(purrr::pmap(list(lambda, Yexp, d), function(x, y, z) {prod((x*y)^z)}))
  
  for (i in 1:n) {  
    int0[i] = integrate(int_tao, lower = 0, upper = Inf, stop.on.error = FALSE,
                        i = i, est.tht = est.tht, A = AA, B = BB, D = DD, frailty = frailty, mode = 0)$value
  }
  
  return(sum(log(int0)))
}