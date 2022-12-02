

FaccAc <- function(y, X, d, coef, lambda, th, penalty = NULL, tune = NULL) {
  
  a = dim(X)[1]
  b = dim(X)[2]
  p = dim(X)[3]
  N = a*b
  
  coef = matrix(coef)
  
  if (dim(lambda)[1] != a || dim(lambda)[2] != b) {
    stop("Dimension of y, X, d or lambda is incorrect.")
  }
  
  if (length(coef) != p) {
    stop("Dimension of coef is incorrect.")
  }
  
  La = matrix(0, a, b)
  Ypre = matrix(0, a, b)
  Xabs = matrix(0, a, b)
  D = rowSums(d)
  
  Z = array(0, c(a, b, p))
  U = array(0, c(a, b, p))
  W = matrix(0, a, b)
  
  for (i in 1:a) {
    for (j in 1:b) {
      La[i, j] = sum((y <= y[i, j]) * lambda)   
    }
  }
  
  for (k in 1:b) {
    Ypre[, k] = X[, k, ] %*% coef
  }
  
  for (k in 1:p) {
    Xabs[, ] = Xabs[, ] + abs(X[,, k])
  }
  
  A = 1/th + D
  C = 1/th + rowSums(La * exp(Ypre))
  YpreExp = exp(Ypre)
  
  for (k in 1:p) {
    temp1 = X[,, k] * YpreExp
    temp2 = abs(X[,, k]) * YpreExp * Xabs
    
    Z[,,k] = FaccCal(y, temp1, Z[,,k], A, C)
    U[,,k] = FaccCal(y, temp2, U[,,k], A, C)
  }
  
  W = FaccCal(y, YpreExp, W, A, C)

  # Update parameters for Gamma Frailty
  Q1 = a*(digamma(1/th) + log(th)-1)/(th^2) + sum(A/C - digamma(A) + log(C))/(th^2) 
  Q2 = a*(3 - 2*digamma(1/th) - 2*log(th))/(th^3) + 2*sum(digamma(A)-log(C)-A/C)/(th^3) - a*trigamma(1/th)/(th^4)
  
  th1 = th - Q1/Q2
  if (th1 > 0) 
    th <- th1
  
  # Updating Coefficients 
  for (k in 1:p) {
    temp1 = sum(d*X[,, k]) - sum(d*Z[,, k]/W)
    temp2 = -sum(d*U[,, k]/W)
    
    if (!is.null(penalty)) {
      if (penalty == "LASSO") {
        temp1 = temp1 - N*sign(coef[k])*tune
        temp2 = temp2 - N*tune/abs(coef[k])
      }
      if (penalty == "MCP") {
        temp1 = temp1 - N*sign(coef[k])*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)
        temp2 = temp2 - N*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)/abs(coef[k])
      }
      if (penalty == "SCAD") {
        temp1 = temp1 - N*sign(coef[k])*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)
        temp2 = temp2 - N*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)/abs(coef[k])
      }
    }
    
    coef[k, 1] = coef[k, 1] - temp1 / temp2
  }
 
  return(list(coef = coef, lambda = lambda, th = th))
}


FaccAc2E <- function(y, X, d, coef, la1, la2, th) {
  
  n = dim(X)[2]
  q = dim(X)[3]
  
  coef = matrix(coef)
  
  La1 = La2 = rep(0,n)
  for(i in 1:n){
    La1[i] = sum(la1*(y[1,] <= y[1,i]))
    La2[i] = sum(la2*(y[2,] <= y[2,i]))
  }
  
  # Updating la1 and la2
  Ypre1 = X[1,,] %*% coef
  Ypre2 = X[2,,] %*% coef
  
  YpreExp1 = exp(Ypre1)
  YpreExp2 = exp(Ypre2)
  
  C0 = La1*YpreExp1 + La2*YpreExp2
  C = 1/th + C0
  
  A = 1/th + colSums(d)
  AC = A/C
  
  E_01 = AC*YpreExp1
  SUM_01 = cumsum((E_01[order(y[1,])])[seq(n, 1, -1)])
  SUM_01 = (SUM_01[seq(n, 1, -1)])[rank(y[1,])] 
  
  la1 = d[1,]/SUM_01
  if(any(is.na(la1))) la1 = rep(1/n, n)
  
  E_02 = AC*YpreExp2
  SUM_02 = cumsum((E_02[order(y[2,])])[seq(n, 1, -1)])
  SUM_02 = (SUM_02[seq(n, 1, -1)])[rank(y[2,])] 
  
  la2 = d[2,]/SUM_02
  if(any(is.na(la2))) la2 = rep(1/n, n)
  
  # Updating Coefficients 
  for(p in 1:q){
    E1 =  La1*X[1,,p]*YpreExp1 + La2*X[2,,p]*YpreExp2
    DE_1 = sum(d*X[,,p]) - sum(AC*E1)  
    
    AVE_X = apply(abs(X), c(1,2), sum)
    E2 = La1*abs(X[1,,p])*AVE_X[1,]*YpreExp1 + La2*abs(X[2,,p])*AVE_X[2,]*YpreExp2
    DE_2 = -2*sum((colSums(d) + 2/th)*E2/C)   
    
    coef[p] = coef[p] - DE_1/DE_2
  }    
  
  # Updating th 
  Q01 = n*(digamma(1/th) + log(th) - 2)/(th^2) + sum(log(C) + colSums(d)/C - digamma(A) + 2/(C*th))/(th^2) + sum(C0/C)/(th^2)
  Q02 = n*(4 - 2*digamma(1/th) - trigamma(1/th)/th - log(th))/(th^3) +
    2*sum(trigamma(A)/(2*th) + digamma(A) - log(C) - colSums(d)/C - 3/(C*th))/(th^3) -
    5*sum( C0/C )/(th^3)
  
  th1 = th - Q01/Q02
  if(th > 0) {th = th1}
  
  return(list(coef = coef, la1 = la1, la2 = la2, th = th))
}