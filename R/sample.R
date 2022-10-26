
sample <- function(coef = matrix(c(-2,3,-4,5,-6,7,6,5,4,3)), lambda = 5, frailty = "LogN", rho = 0.2, init.var = 10, a = 40, b = 20, cen = 2) {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  C_0 = matrix(rho, p, p)
  Row = matrix(c(1:p), p, p)
  Col = t(Row)
  Sigma = C_0^(abs(Row - Col))
  
  if (frailty == "Gamma") {
    th = init.var
    u = rgamma(a, 1/th, scale = 1/th)
  }
  
  if (frailty == "LogN") {
    si = sqrt(init.var)
    u = rlnorm(a, 0, si)  
  }
  
  if (frailty == "InvGauss") {
    mu = sqrt(init.var)
    u = rinvGauss(a, mu, mu^2)
  }
  
  X = array(0, c(a, b, p))
  T = matrix(0, a, b)  
  
  for(i in seq_len(a)) {
    X[i,,] = mvrnorm(b, rep(0, p), Sigma)
    U = runif(b, 0, 1)
    T[i,] <- -log(U)/(lambda*u[i]*exp(X[i,,] %*% coef))    
  }
  
  d = 1*(T <= cen)
  y = pmin(T, cen)
  return(list(y = y, X = X, d = d))
}