

MEGammaFrailty <- function(y, X, d, coef, la1, la2, th) {
  
  n = dim(X)[2]
  q = dim(X)[3]
  
  ell=rep(0, 10000)
  iter = 1
  ell[iter] = GammaLik2E(y, X, d, coef, la1, la2, th)
  error = 3

  while(error > 0.0000001)  {
    
    temp1 = FaccAc2E(y, X, d, coef, la1, la2, th)
    la1 = temp1$la1
    la2 = temp1$la2
    coef1 = temp1$coef
    th = temp1$th
    
    u_be = coef1 - coef
    
    temp2 = FaccAc2E(y, X, d, coef1, la1, la2, th)
    la1 = temp2$la1
    la2 = temp2$la2
    coef2 = temp2$coef
    th = temp2$th
    
    v_be = coef2 - 2*coef1 + coef
    
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if(al_be > -1) {al_be = -1}
    
    coef = coef2 - 2*al_be*u_be + al_be^2*v_be
    
    La1 = La2 = rep(0, n)
    for(i in 1:n){
      La1[i] = sum(la1*(y[1,] <= y[1,i]))
      La2[i] = sum(la2*(y[2,] <= y[2,i]))
    }
    
    # Update lambda Variables
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
    
    # Update th
    Q01 = n*(digamma(1/th) + log(th) - 2)/(th^2) + sum(log(C) + colSums(d)/C - digamma(A) + 2/(C*th))/(th^2) + sum(C0/C)/(th^2)
    Q02 = n*(4 - 2*digamma(1/th) - trigamma(1/th)/th - log(th))/(th^3) + 
      2*sum(trigamma(A)/(2*th) + digamma(A) - log(C) - colSums(d)/C - 3/(C*th))/(th^3) -
      5*sum(C0/C)/(th^3)
    
    th1 = th - Q01/Q02
    if(th1 > 0) {th = th1}
    
    ell[iter + 1] = GammaLik2E(y, X, d, coef, la1, la2, th)
    
    error = abs(ell[iter + 1] - ell[iter])/(1 + abs(ell[iter]))
    
    iter = iter + 1
    cat(error, '\n')
  }
  
  return(list(coef = coef, th = th, likelihood = ell[iter], iteration = iter))
  
}