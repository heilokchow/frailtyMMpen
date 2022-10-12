library(survival)

n = 1
estMM = matrix(0,n,15)

sample <- function(la,alp,th,be1,be2,be3,be4,be5,be6,be7,be8,be9,be10,a,b)
{
  u = rgamma(a, alp, scale = th) # HH
  x1 = matrix(0,a,b)
  x2 = matrix(0,a,b)
  x3 = matrix(0,a,b)
  x4 = matrix(0,a,b)
  x5 = matrix(0,a,b)
  x6 = matrix(0,a,b)
  x7 = matrix(0,a,b)
  x8 = matrix(0,a,b)
  x9 = matrix(0,a,b)
  x10 = matrix(0,a,b)
  
  T = matrix(0,a,b)  
  cen = 10
  for(i in 1:a)
  {
    x1[i,]=runif(b, min = 0, max = 0.5)
    x2[i,]=runif(b, min = 0, max = 0.5)
    x3[i,]=runif(b, min = 0, max = 0.5)
    x4[i,]=runif(b, min = 0, max = 0.5)
    x5[i,]=runif(b, min = 0, max = 0.5)
    x6[i,]=runif(b, min = 0, max = 0.5)
    x7[i,]=runif(b, min = 0, max = 0.5)
    x8[i,]=runif(b, min = 0, max = 0.5)
    x9[i,]=runif(b, min = 0, max = 0.5)
    x10[i,]=runif(b, min = 0, max = 0.5)
    
    
    U = runif(b,0,1)
    T[i,] <- -log(U)/(la*u[i]*exp(x1[i,]*be1+x2[i,]*be2+x3[i,]*be3+x4[i,]*be4+x5[i,]*be5+x6[i,]*be6+x7[i,]*be7+x8[i,]*be8+x9[i,]*be9+x10[i,]*be10))    
  }
  d = 1*(T<=cen)
  y=pmin(T,cen)
  return(list(y=y,d=d,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10))
}


set.seed(10)
for (r in 1:n) {
  yy = sampleGamma(seed = r*1000)

  a = dim(yy$X)[1]
  b = dim(yy$X)[2]
  p = dim(yy$X)[3]
  lambda = matrix(1, a, b)/(a * b)
  
  start = proc.time()[1]
  
  # yy = sample(5,0.1,0.5,-2,3,-4,5,-6,7,6,5,4,3,500,5)
  # y = yy$y 
  # d = yy$d
  # X = array(c(yy$x1, yy$x2, yy$x3, yy$x4, yy$x5, yy$x6, yy$x7, yy$x8, yy$x9, yy$x10), dim = c(dim(y), 10))
  # 
  # a = 500
  # b = 5
  y = yy$y
  X = yy$X
  d = yy$d
  th = 1
  coef = rep(1, p)
  lambda = matrix(1, a, b)/(a*b)
  ret = CLGammaFrailty(y, X, d, coef, lambda, th)
  
  end = proc.time()[1]
  time = end-start
  
  estMM[r,] = c(ret$iteration, time, ret$likelihood, ret$th, ret$coef, ret$mLambda)
}


aveMM=c(mean(estMM[,1]),mean(estMM[,2]),mean(estMM[,3]),mean(estMM[,4]),mean(estMM[,5]),mean(estMM[,6]),mean(estMM[,7]),mean(estMM[,8]),mean(estMM[,9]),mean(estMM[,10]),mean(estMM[,11]),mean(estMM[,12]),mean(estMM[,13]),mean(estMM[,14]),mean(estMM[,15]),mean(estMM[,16]))
stdMM=c(sd(estMM[,1]),sd(estMM[,2]),sd(estMM[,3]),sd(estMM[,4]),sd(estMM[,5]),sd(estMM[,6]),sd(estMM[,7]),sd(estMM[,8]),sd(estMM[,9]),sd(estMM[,10]),sd(estMM[,11]),sd(estMM[,12]),sd(estMM[,13]),sd(estMM[,14]),sd(estMM[,15]),sd(estMM[,16]))
biasMM = aveMM - c(0,0,0,0.1,0.1,-2,3,-4,5,-6,7,6,5,4,3,5)

aveMM
stdMM
biasMM

# Dump
X[,,1] = x1
X[,,2] = x2
X[,,3] = x3
X[,,4] = x4
X[,,5] = x5
X[,,6] = x6
X[,,7] = x7
X[,,8] = x8
X[,,9] = x9
X[,,10] = x10
