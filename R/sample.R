sample_CL <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 26))), lambda = 5, frailty = "LogN", power = NULL, init.var = 10, a = 50, b = 10, cen = 2) {
  
  p = length(coef)
  coef = as.matrix(coef)
  N = a * b
  
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
  
  X1 = aperm(X, c(2,1,3))
  outX = as.data.frame(matrix(c(X1), nrow = N, ncol = p))
  outX$time = as.vector(t(y))
  outX$status = as.vector(t(d))
  outX$id = sort(rep(seq(1, a, 1), b))
  
  return(list(data = outX, y = y, X = X, d = d))
}


sample_ME <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 26))), lambda1 = 3, lambda2 = 5, frailty = "LogN", power = NULL, init.var = 10, n = 200, cen = 2) {
  
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
  
  if (frailty == "PVF") {
    u = rtweedie(a, mu = 1, phi = init.var, power = power)
  }
  
  X = array(runif(2*n*p, 0, 0.5), c(n, p, 2))
  T = matrix(0, n, 2)  
  cen1 = matrix(0, n, 2)  
  
  cen1[,1] = runif(n, 0, cen)
  cen1[,2] = runif(n, 0, cen)
  
  U1 = runif(n, 0, 1)
  U2 = runif(n, 0, 1)
  
  T[,1] <- -log(U1)/(lambda1*u*exp(X[,,1] %*% coef)) 
  T[,2] <- (exp(-log(U2)/(u*exp(X[,,2] %*% coef))) - 1)/lambda2
  
  d = 1*(T <= cen1)
  y = pmin(T, cen1)
  
  outX = as.data.frame(rbind(X[,,1], X[,,2]))
  outX$time = as.vector(y)
  outX$status = as.vector(d)
  outX$id = c(rep(1, n), rep(2, n))
  
  return(list(data = outX, y = y, X = X, d = d))
}


sample_RE <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 26))), lambda = 5, frailty = "LogN", power = NULL, init.var = 1, n = 50, cen = 400) {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  cen = matrix(runif(n, 0.001 * cen, 0.002 * cen), n)
  
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
  
  if (frailty == "PVF") {
    u = rtweedie(n, mu = 1, phi = init.var, power = power)
  }
  
  y = list()
  y1 = list()
  X = list()
  d = list()
  id = list()
  
  for(i in seq_len(n)) {
    temps = 0
    y[[i]] = c(0)
    d[[i]] = c(0)
    X[[i]] = matrix(rep(0, p), nrow = 1)
    while (temps < cen[i,]) {
      tempX = runif(p, 0, 0.5)
      U = runif(1, 0, 1)
      tempT = -log(U)/(lambda*u[i]*exp(tempX %*% coef)) 
      temps = temps + tempT
      y[[i]] = c(y[[i]], temps)
      d[[i]] = c(d[[i]], 1)
      X[[i]] = rbind(X[[i]], tempX)
    }
    li = length(y[[i]])
    y[[i]][li] = cen[i,]
    d[[i]][li] = 0
    
    y1[[i]] = y[[i]][1:(li-1)]
    y[[i]] = y[[i]][2:li]
    d[[i]] = d[[i]][2:li]
    X[[i]] = X[[i]][2:li,, drop = FALSE]
    id[[i]] = rep(i, li-1)
  }
  
  out = unlist(y)
  N = length(out)
  
  lambda = y
  for (i in 1:n) {
    for (j in seq_len(length(y[[i]]))) {
      lambda[[i]][j] = 1/N
    }
  }
  
  outX = matrix(0, nrow = 0, ncol = p)
  for (i in seq_len(n)) {
    outX = rbind(outX, X[[i]])
  }
  outX = as.data.frame(outX)
  outX$start = unlist(y1)
  outX$end = unlist(y)
  outX$status = unlist(d)
  outX$id = unlist(id)
  
  return(list(outX, y, X, d, lambda))
}


