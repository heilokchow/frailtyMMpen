MMprocess_CL <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL) {
  
  p = length(coef)
  coef = as.matrix(coef)
  
  a = nrow(y)
  b = ncol(y)
  N = a*b
  vy = as.vector(y)
  vd = as.vector(d)
  
  La = (cumsum(lambda[order(vy)]))[rank(vy)]
  La = matrix(La, a, b)
  
  YpreExp = matrix(0, a, b)
  for (i in 1:a) {
    YpreExp[i,] = exp(X[i,,] %*% coef)
  }
  
  A = rowSums(La*YpreExp)
  B = apply((lambda*YpreExp)^d, 1, prod)
  D = rowSums(d) 
  
  int0 <- vector("numeric", length = a)
  int1 <- vector("numeric", length = a)
  
  if (frailty == "LogN" || frailty == "InvGauss" || frailty == "PVF") {
    for (i in 1:a) {  
      int0[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = A, B = B, D = D, frailty = frailty, power = power, mode = 0)$value
    }
    
    if (any(int0 == 0)) {
      return(list(coef = rep(1/p, p), est.tht = est.tht, lambda = rep(1/N, N)))
    }
    
    for (i in 1:a) {  
      int1[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, power = power, mode = 1)$value
    }
  }
  
  if (frailty == "Gamma") {
    C2 = 1/est.tht + A
    A2 = 1/est.tht + D
    int1 = A2/C2
  }
  
  # Update lambda Variables
  ME = matrix(int1, a, b)
  E_0 = as.vector(ME*YpreExp)
  SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
  SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
  
  lambda = vd/SUM_0
  
  # Update coefficients
  AVE_X = apply(abs(X), c(1,2), sum)
  for(k in 1:p) {
    E_1 = as.vector(ME*X[,,k]*YpreExp)
    E_2 = as.vector(ME*AVE_X*abs(X[,,k])*YpreExp)
    
    SUM_1 = cumsum((E_1[order(vy)])[seq(N, 1, -1)])
    SUM_1 = (SUM_1[seq(N, 1, -1)])[rank(vy)] 
    SUM_2 = cumsum((E_2[order(vy)])[seq(N, 1, -1)])
    SUM_2 = (SUM_2[seq(N, 1, -1)])[rank(vy)] 
    
    DE_1 = sum(d*X[,,k]) - sum(vd*SUM_1/SUM_0) 
    DE_2 = -sum(vd*SUM_2/SUM_0) 
    
    if (!is.null(penalty)) {
      if (penalty == "LASSO") {
        DE_1 = DE_1 - N*sign(coef[k])*tune
        DE_2 = DE_2 - N*tune/abs(coef[k])
      }
      if (penalty == "MCP") {
        DE_1 = DE_1 - N*sign(coef[k])*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)
        DE_2 = DE_2 - N*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)/abs(coef[k])
      }
      if (penalty == "SCAD") {
        DE_1 = DE_1 - N*sign(coef[k])*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)
        DE_2 = DE_2 - N*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)/abs(coef[k])
      }
    }
    
    coef[k] = coef[k] - DE_1/DE_2
  }
  
  # Update frailty parameters
  if (frailty == "PVF") {
    gn = function(target, est.tht, power, A, B, D, int0) {
      int_f = function(x, i, target, est.tht, A, B, D, tao0) {
        temp0 = dtweedie.dldphi(y = x, mu = 1, phi = target, power = power)*
          dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)/tao0
        return(temp0)
      }
      
      s = 0
      for (i in seq_len(length(A))) {
        temp = integrate(int_f, lower = 0, upper = 20, stop.on.error = FALSE,
                         i = i, target = target, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i])$value
        s = s + temp
      }
      return(s)
    }
    
    if (gn(max(0.2, est.tht - 3), est.tht, power, A, B, D, int0) > 0) {
      est.tht = max(0.2, est.tht - 3)
    } else if (gn(min(est.tht + 3, 20), est.tht, power, A, B, D, int0) < 0) {
      est.tht = min(est.tht + 3, 20)
    } else {
      est.tht = uniroot(gn, c(max(0.2, est.tht - 3), min(est.tht + 3, 20)), est.tht = est.tht, power = power, A = A, B = B, D = D, int0 = int0, tol = 1e-3)$root
    }
  }
  
  if (frailty == "LogN" || frailty == "InvGauss") {
    int2 <- vector("numeric", length = a)
    for (i in 1:a) {  
      int2[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, mode = 2)$value
    }
    est.tht = sum(int2)/a
  }

  if (est.tht < 0) {
    est.tht = 1
  } 
  
  if (frailty == "Gamma") {
    Q01 = a*(digamma(1/est.tht)+log(est.tht)-1)/(est.tht^2) + sum(A2/C2-digamma(A2)+log(C2))/(est.tht^2) 
    Q02 = a*(3-2*digamma(1/est.tht)-2*log(est.tht))/(est.tht^3)+2*sum(digamma(A2)-log(C2)-A2/C2)/(est.tht^3) - a*trigamma(1/est.tht)/(est.tht^4)
    
    est.tht1 = est.tht - Q01/Q02
    if(est.tht1 > 0) {
      est.tht = est.tht1 
    }
  }

  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}


MMprocess_ME <- function(y, X, d, coef, lambda1, lambda2, est.tht, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL) {
  
  p = length(coef)
  coef = as.matrix(coef)
  n = nrow(y)
  
  La1 = La2 = rep(0,n)
  for(i in 1:n){
    La1[i] = sum(lambda1*(y[,1] <= y[i,1]))
    La2[i] = sum(lambda2*(y[,2] <= y[i,2]))
  }
  
  YpreExp = matrix(0, n, 2)
  for (i in 1:2) {
    YpreExp[,i] = exp(X[,,i] %*% coef)
  }
  
  CC =  La1*YpreExp[,1] + La2*YpreExp[,2]
  D = colSums(d)
  d1 = d[,1]
  d2 = d[,2]
  
  int0 <- vector("numeric", length = n)
  int1 <- vector("numeric", length = n)

  if (frailty == "LogN" || frailty == "InvGauss" || frailty == "PVF") {
    AA = (lambda1*YpreExp[,1])^(d1)*(lambda2*YpreExp[,2])^(d2)
    
    for (i in 1:n) {  
      int0[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = CC, B = AA, D = D, frailty = frailty, power = power, mode = 0)$value
    }
    
    if (any(int0 == 0)) {
      return(list(coef = rep(1/p, p), est.tht = est.tht, lambda1 = rep(1/n, n), lambda2 = rep(1/n, n)))
    }
    
    for (i in 1:n) {  
      int1[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = CC, B = AA, D = D, tao0 = int0[i], frailty = frailty, power = power, mode = 1)$value
    }
  }
  
  if (frailty == "Gamma") {
    C = 1/est.tht + CC
    A = 1/est.tht + D
    int1 = A/C
  }
  
  # Update lambda Variables
  E_01 = int1*YpreExp[,1]
  SUM_01 = cumsum((E_01[order(y[,1])])[seq(n,1,-1)])
  SUM_01 = (SUM_01[seq(n,1,-1)])[rank(y[,1])] 
  
  lambda1 = d[,1]/SUM_01
  if(any(is.na(lambda1))) lambda1 = rep(1/n, n)
  
  E_02 = int1*YpreExp[,2]
  SUM_02 = cumsum((E_02[order(y[,2])])[seq(n,1,-1)])
  SUM_02 = (SUM_02[seq(n,1,-1)])[rank(y[,2])] 
  
  lambda2 = d[,2]/SUM_02
  if(any(is.na(lambda2))) lambda2 = rep(1/n, n)
  

  # Update coefficients
  AVE_X = apply(abs(X), c(1,3), sum)
  for(k in 1:p)
  {
    E1 = La1*X[,k,1]*YpreExp[,1] + La2*X[,k,2]*YpreExp[,2]
    DE_1 = sum(d*X[,k,]) - sum(int1*E1)   
    
    E2 = La1*abs(X[,k,1])*AVE_X[,1]*YpreExp[,1] + La2*abs(X[,k,2])*AVE_X[,2]*YpreExp[,2]
    
    if (frailty == "LogN" || frailty == "InvGauss") {
      DE_2 = -2*sum(int1*E2)  
    }
    
    if (frailty == "Gamma") {
      DE_2 = -2*sum((D+2/est.tht)*E2/C) 
    }
    
    if (!is.null(penalty)) {
      if (penalty == "LASSO") {
        DE_1 = DE_1 - n*sign(coef[k])*tune
        DE_2 = DE_2 - n*tune/abs(coef[k])
      }
      if (penalty == "MCP") {
        DE_1 = DE_1 - n*sign(coef[k])*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)
        DE_2 = DE_2 - n*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)/abs(coef[k])
      }
      if (penalty == "SCAD") {
        DE_1 = DE_1 - n*sign(coef[k])*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)
        DE_2 = DE_2 - n*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)/abs(coef[k])
      }
    }
    
    coef[k] = coef[k] - DE_1/DE_2
  }     
  
  # Update frailty parameters
  if (frailty == "PVF") {
    gn = function(target, est.tht, power, A, B, D, int0) {
      int_f = function(x, i, target, est.tht, A, B, D, tao0) {
        temp0 = dtweedie.dldphi(y = x, mu = 1, phi = target, power = power)*
          dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)/tao0
        return(temp0)
      }
      
      s = 0
      for (i in seq_len(length(A))) {
        temp = integrate(int_f, lower = 0, upper = 20, stop.on.error = FALSE,
                         i = i, target = target, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i])$value
        s = s + temp
      }
      return(s)
    }
    
    if (gn(max(0.2, est.tht - 3), est.tht, power, CC, AA, D, int0) > 0) {
      est.tht = max(0.2, est.tht - 3)
    } else if (gn(min(est.tht + 3, 20), est.tht, power, CC, AA, D, int0) < 0) {
      est.tht = min(est.tht + 3, 20)
    } else {
      est.tht = uniroot(gn, c(max(0.2, est.tht - 3), min(est.tht + 3, 20)), est.tht = est.tht, power = power, A = CC, B = AA, D = D, int0 = int0, tol = 1e-3)$root
    }
  }
  
  if (frailty == "LogN" || frailty == "InvGauss") {
    int2 <- vector("numeric", length = n) 
    
    for (i in 1:n) {  
      int2[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = CC, B = AA, D = D, tao0 = int0[i], frailty = frailty, power = power, mode = 2)$value
    }
    
    est.tht = sum(int2)/n
    if (est.tht < 0) {
      est.tht = 1
    }
  }
  
  if (frailty == "Gamma") {
    Q01 = n*(digamma(1/est.tht)+log(est.tht)-2)/(est.tht^2) + 
      sum(log(C)+D/C-digamma(A)+2/(C*est.tht))/(est.tht^2) + sum(CC/C)/(est.tht^2)
    Q02 = n*(4-2*digamma(1/est.tht)-trigamma(1/est.tht)/est.tht-log(est.tht))/(est.tht^3)+
      2*sum(trigamma(A)/(2*est.tht)+digamma(A)-log(C)-D/C-3/(C*est.tht))/(est.tht^3) -
      5*sum(CC/C)/(est.tht^3)
    
    est.tht1 = est.tht - Q01/Q02
    if(est.tht1 > 0) {
      est.tht = est.tht1 
    }
  }
  
  return(list(coef = coef, est.tht = est.tht, lambda1 = lambda1, lambda2 = lambda2))
}



MMprocess_RE <- function(y, X, d, coef, lambda, est.tht, frailty = "LogN", power = NULL, penalty = NULL, tune = NULL) {
  
  p = length(coef)
  coef = as.matrix(coef)
  n = length(y)
  
  vy = unlist(y)
  yend = unlist(lapply(y, max))
  vd = unlist(d)
  vl = unlist(lambda)
  
  # Possible C Code
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
  
  if (frailty == "LogN" || frailty == "InvGauss" || frailty == "PVF") {
    
    BB = unlist(purrr::pmap(list(lambda, Yexp, d), function(x, y, z) {prod((x*y)^z)}))
    
    for (i in 1:n) {  
      int0[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = AA, B = BB, D = DD, frailty = frailty, power = power, mode = 0)$value
    }
    
    if (any(int0 == 0)) {
      return(list(coef = rep(1/p, p), est.tht = est.tht, lambda = lapply(lambda, function(x) {x^0*1/n})))
    }
    
    for (i in 1:n) {  
      int1[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = AA, B = BB, D = DD, tao0 = int0[i], frailty = frailty, power = power, mode = 1)$value
    }
  }
  
  if (frailty == "Gamma") {
    C = 1/est.tht + AA
    A = 1/est.tht + DD
    int1 = A/C
  }
  
  # Update lambda Variables
  Xmatexp = apply(Xmat, 3, function(x, coef) {exp(x %*% coef)}, coef = coef)
  SUM_0 = lambda
  
  for (i in 1:n) {
    it = length(y[[i]])
    for (j in 1:it) {
      th = which(tall == y[[i]][j])
      SUM_0[[i]][j] = sum((yend >= y[[i]][j]) * int1 * Xmatexp[,th])
    }
  }
  
  lambda = purrr::pmap(list(d, SUM_0), function(x, y) {x/y})
  
  # Update coefficients
  AVE_X = apply(abs(Xmat), c(1,3), sum)
  for(k in 1:p) {
    
    SUM_1 = lambda
    SUM_2 = lambda
    
    for (i in 1:n) {
      it = length(y[[i]])
      for (j in 1:it) {
        th = which(tall == y[[i]][j])
        SUM_1[[i]][j] = sum((yend >= y[[i]][j]) * int1 * Xmatexp[,th] * Xmat[,k,th]) 
        SUM_2[[i]][j] = sum((yend >= y[[i]][j]) * int1 * Xmatexp[,th] * abs(Xmat[,k,th]) * AVE_X[,th]) 
      }
    }
    
    DE_1 = sum(unlist(purrr::pmap(list(d, X, SUM_0, SUM_1), function(a, b, c, d) {a*b[,k] - a*d/c})))   
    DE_2 = sum(unlist(purrr::pmap(list(d, SUM_0, SUM_2), function(a, c, f) {- a*f/c})))   
    
    if (!is.null(penalty)) {
      if (penalty == "LASSO") {
        DE_1 = DE_1 - n*sign(coef[k])*tune
        DE_2 = DE_2 - n*tune/abs(coef[k])
      }
      if (penalty == "MCP") {
        DE_1 = DE_1 - n*sign(coef[k])*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)
        DE_2 = DE_2 - n*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)/abs(coef[k])
      }
      if (penalty == "SCAD") {
        DE_1 = DE_1 - n*sign(coef[k])*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)
        DE_2 = DE_2 - n*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)/abs(coef[k])
      }
    }
    
    coef[k] = coef[k] - DE_1/DE_2
  }     
  
  # Update frailty parameters
  if (frailty == "PVF") {
    gn = function(target, est.tht, power, A, B, D, int0) {
      int_f = function(x, i, target, est.tht, A, B, D, tao0) {
        temp0 = dtweedie.dldphi(y = x, mu = 1, phi = target, power = power)*
          dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)/tao0
        return(temp0)
      }
      
      s = 0
      for (i in seq_len(length(A))) {
        temp = integrate(int_f, lower = 0, upper = 20, stop.on.error = FALSE,
                         i = i, target = target, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i])$value
        s = s + temp
      }
      return(s)
    }
    
    if (gn(max(0.2, est.tht - 3), est.tht, power, AA, BB, DD, int0) > 0) {
      est.tht = max(0.2, est.tht - 3)
    } else if (gn(min(est.tht + 3, 20), est.tht, power, AA, BB, DD, int0) < 0) {
      est.tht = min(est.tht + 3, 20)
    } else {
      est.tht = uniroot(gn, c(max(0.2, est.tht - 3), min(est.tht + 3, 20)), est.tht = est.tht, power = power, A = AA, B = BB, D = DD, int0 = int0, tol = 1e-3)$root
    }
  }
  
  if (frailty == "LogN" || frailty == "InvGauss") {
    int2 <- vector("numeric", length = n) 
    
    for (i in 1:n) {  
      int2[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
                          i = i, est.tht = est.tht, A = AA, B = BB, D = DD, tao0 = int0[i], frailty = frailty, power = power, mode = 2)$value
    }
    
    est.tht = sum(int2)/n
    if (est.tht < 0) {
      est.tht = 1
    }
  }
  
  if (frailty == "Gamma") {
    Q01 = n*(digamma(1/est.tht)+log(est.tht)-2)/(est.tht^2) + 
      sum(log(C)+DD/C-digamma(A)+2/(C*est.tht))/(est.tht^2) + sum(AA/C)/(est.tht^2)
    Q02 = n*(4-2*digamma(1/est.tht)-trigamma(1/est.tht)/est.tht-log(est.tht))/(est.tht^3)+
      2*sum(trigamma(A)/(2*est.tht)+digamma(A)-log(C)-DD/C-3/(C*est.tht))/(est.tht^3) -
      5*sum(AA/C)/(est.tht^3)
    
    est.tht1 = est.tht - Q01/Q02
    if(est.tht1 > 0) {
      est.tht = est.tht1 
    }
  }
  
  return(list(coef = coef, est.tht = est.tht, lambda = lambda))
}
