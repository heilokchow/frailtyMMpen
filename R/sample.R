
sample_CL <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 46))), lambda = 5, frailty = "LogN", power = NULL, init.var = 10, a = 50, b = 10, cen = 2) {
  
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



sample_ME <- function(coef = matrix(c(1, 2, 3, 4, rep(0, 46))), lambda1 = 3, lambda2 = 5, frailty = "LogN", init.var = 10, n = 300, cen = 2) {
  
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
  cen = matrix(0, 2, n)  
  
  cen[1,] = runif(n, 0, 1)
  cen[2,] = runif(n, 0, 1)
  
  U1 = runif(n, 0, 1)
  U2 = runif(n, 0, 1)
  
  T[1,] <- -log(U1)/(lambda1*u*exp(X[1,,] %*% coef)) 
  T[2,] <- (exp(-log(U2)/(u*exp(X[2,,] %*% coef))) - 1)/lambda2
  
  d = 1*(T <= cen)
  y = pmin(T, cen)
  return(list(y = y, X = X, d = d))
}