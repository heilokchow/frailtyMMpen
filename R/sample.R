
sample_CL <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 26))), lambda = 5, frailty = "LogN", power = NULL, init.var = 10, a = 50, b = 10, cen = 2) {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  cen = matrix(runif(a*b, 0.001 * cen, 0.002 * cen), a, b)
  
  if (frailty == "Gamma") {
    th = init.var
    u = rgamma(a, 1/th, scale = 1/th)
  }
  
  if (frailty == "LogN") {
    si = sqrt(init.var)
    u = rlnorm(a, 0, si)  
  }
  
  if (frailty == "InvGauss") {
    alp = 1/init.var
    u = SuppDists::rinvGauss(a, 1, alp)
  }
  
  if (frailty == "PVF") {
    u = rtweedie(a, mu = 1, phi = init.var, power = power)
  }
  
  X = array(runif(a*b*p, 0, 0.5), c(a, b, p))
  T = matrix(0, a, b)  
  
  for(i in seq_len(a)) {
    U = runif(b, 0, 1)
    T[i,] <- -log(U)/(lambda*u[i]*exp(X[i,,] %*% coef))    
  }
  
  d = 1*(T <= cen)
  y = pmin(T, cen)
  return(list(y = y, X = X, d = d))
}



sample_ME <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 26))), lambda1 = 3, lambda2 = 5, frailty = "LogN", init.var = 10, n = 200, cen = 2) {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  if (frailty == "Gamma") {
    th = init.var
    u = rgamma(n, 1/th, scale = 1/th)
  }
  
  if (frailty == "LogN") {
    si = sqrt(init.var)
    u = rlnorm(n, 0, si)  
  }
  
  if (frailty == "InvGauss") {
    alp = 1/init.var
    u = SuppDists::rinvGauss(n, 1, alp)
  }
  
  X = array(runif(2*n*p, 0, 0.5), c(2, n, p))
  T = matrix(0, 2, n)  
  cen1 = matrix(0, 2, n)  
  
  cen1[1,] = runif(n, 0, cen)
  cen1[2,] = runif(n, 0, cen)
  
  U1 = runif(n, 0, 1)
  U2 = runif(n, 0, 1)
  
  T[1,] <- -log(U1)/(lambda1*u*exp(X[1,,] %*% coef)) 
  T[2,] <- (exp(-log(U2)/(u*exp(X[2,,] %*% coef))) - 1)/lambda2
  
  d = 1*(T <= cen1)
  y = pmin(T, cen1)
  return(list(y = y, X = X, d = d))
}


sample_RE <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 26))), lambda = 5, frailty = "LogN", power = NULL, init.var = 10, a = 50, b = 10, cen = 500) {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  cen = matrix(runif(a, 0.001 * cen, 0.002 * cen), a)
  
  if (frailty == "Gamma") {
    th = init.var
    u = rgamma(a, 1/th, scale = 1/th)
  }
  
  if (frailty == "LogN") {
    si = sqrt(init.var)
    u = rlnorm(a, 0, si)  
  }
  
  if (frailty == "InvGauss") {
    alp = 1/init.var
    u = SuppDists::rinvGauss(a, 1, alp)
  }
  
  if (frailty == "PVF") {
    u = rtweedie(a, mu = 1, phi = init.var, power = power)
  }
  
  X1 = array(runif(a*b*p, 0, 0.5), c(a, b, p))
  T1 = matrix(0, a, b)  
  
  for(i in seq_len(a)) {
    U = runif(b, 0, 1)
    T1[i,] <- -log(U)/(lambda*u[i]*exp(X1[i,,] %*% coef))    
  }
  
  if (b == 1) {
    T2 = T1
  } else {
    T2 = t(apply(T1, 1, cumsum))
  }
  
  y = list()
  X = list()
  d = list()
  
  for(i in seq_len(a)) {
    lst = b
    fl = 0
    pos = which(T2[i,] > cen[i,])
    if (length(pos) != 0) {
      lst = pos[1]
      if (lst == 1) {
        lstlen = cen[i,]
      } else {
        lstlen = cen[i,]
      }
    } else {
      fl = 1
      lst = b
      lstlen = T2[i, lst]
    }
    
    y[[i]] = c(T2[i, seq_len(lst-1)], lstlen)
    d[[i]] = c(rep(1, lst-1), fl)
    X[[i]] = matrix(X1[i,seq_len(lst),], ncol = p)
  }
  
  return(list(y = y, X = X, d = d))
}


