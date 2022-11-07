
sample <- function(coef = matrix(c(rep(-2,6), rep(-1,6), rep(1,6), rep(2,6), rep(3,6))), lambda = 5, frailty = "LogN", init.var = 10, a = 50, b = 10, cen = 2) {
  
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
    u = rinvGauss(a, 1, alp)
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