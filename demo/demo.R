library(survival)

n = 1
estMM = matrix(0,n,15)

sample <- function(la,th,be1,be2,be3,be4,be5,be6,be7,be8,be9,be10,a,b)
{
  u = rgamma(a, 1/th, rate=1/th) # HH
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
  cen = 2
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


set.seed(1)
for (r in 1:n) {
  # yy = sampleGamma(seed = r*1000)
  # 
  # a = dim(yy$X)[1]
  # b = dim(yy$X)[2]
  # p = dim(yy$X)[3]
  # lambda = matrix(1, a, b)/(a * b)
  # y = yy$y
  # X = yy$X
  # d = yy$d
  
  start = proc.time()[1]
  
  yy = sample(5,10,-2,3,-4,5,-6,7,6,5,4,3,100,5)
  y = yy$y
  d = yy$d
  X = array(c(yy$x1, yy$x2, yy$x3, yy$x4, yy$x5, yy$x6, yy$x7, yy$x8, yy$x9, yy$x10), dim = c(dim(y), 10))
   
  a = 100
  b = 5
  th = 1
  coef = rep(1, 10)
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

# MEGamma

nsam <- function(a1,a2,th,be,n) {
  u = rgamma(n, 1/th, rate=1/th)
  
  q = length(be)
  
  x = array( runif(2*n*q, min=0, max=0.5), dim=c(2,n,q) )
  
  Be = matrix(rep(be,each=n),n,q)
  
  T = matrix(0,2,n)
  
  cen = runif(n,0,3)    ## (0,3)/100=10%  
  
  
  U1 = runif(n,0,1)    
  U2 = runif(n,0,1)
  
  T[1,] <- -log(U1)/(a1*u*exp(apply(x[1,,]*Be,1,sum))) 
  T[2,] <- ( exp( -log(U2)/(u*exp(apply(x[2,,]*Be,1,sum))) )-1 )/a2
  
  d = 1*(T<=cen)
  y=pmin(T,cen)
  
  la1 = a1*y[1,]
  la2 = log(1+a2*y[2,])
  
  return(list(y=y,d=d,x=x,la1=la1,la2=la2,w=u))
}

N = 1

mm.gamma = matrix(0,N,23)
set.seed(1)

for(j in 1:N){
  
  be0 = c(rep(-2,10), rep(3,10))
  q = length(be0)
  n = 100
  
  yy = nsam(3,5,1,be0,n)
  y = yy$y 
  d = yy$d
  x = yy$x

  la1 = rep(1/n, n)
  la2 = rep(1/n, n)
  th = 0.5
  be = rep(0.5, q)
  coef = be
  X = x
  
  start = proc.time()[1]
  rs = MEGammaFrailty(y, x, d, be, la1, la2, th)
  end = proc.time()[1]
  
  mm.gamma[j,] = c(rs$coef, rs$th, end - start, rs$likelihood) 
}



Time = apply(mm.gamma,2,mean)[22]
Time

# MLE 
apply(mm.gamma,2,mean)[c(1,5,10,15,20,21)]

# BIAS 
abs(apply(mm.gamma,2,mean)[c(1,5,10,15,20,21)] - c(rep(-2,3), rep(3,2),1))

# SD 
apply(mm.gamma,2,sd)[c(1,5,10,15,20,21)]



# InvGauss and LogN -------------------------------------------------------

set.seed(2)
yy = sample()
y = yy$y 
d = yy$d
X = yy$X

a = nrow(y)
b = ncol(y)
N = a*b
vy = as.vector(y)
vd = as.vector(d)

est.tht = 1
coef = rep(0.5, 10)
lambda = rep(1/N, N)

ell=rep(0,1000000)
k=1

ell[k]= logLikihood(y, X, d, coef, lambda, est.tht, frailty = "LogN")
error = 3

start = proc.time()[1]
while(error > 0.000001) {
  
  rs = EMprocess(y, X, d, coef, lambda, est.tht, frailty = "LogN") 
  coef = rs$coef
  lambda = rs$lambda
  est.tht = rs$est.tht
  
  ell[k+1] = logLikihood(y, X, d, coef, lambda, est.tht, frailty = "LogN")
  
  error = abs(ell[k+1]-ell[k])/(1+abs(ell[k]))
  cat(error, '\n')
  k = k+1
}
end = proc.time()[1]
time = end - start

# SQS3

ell=rep(0,1000000)
k=1

ell[k]= logLikihood(y, X, d, coef, lambda, est.tht, frailty = "LogN")
error = 3

start = proc.time()[1]
while(error > 0.000001) {
  
  rs1 = EMprocess(y, X, d, coef, lambda, est.tht, frailty = "LogN") 
  coef1 = rs1$coef
  est.tht = rs1$est.tht
  lambda = rs1$lambda
  
  u_be = coef1 - coef
  
  rs2 = EMprocess(y, X, d, coef1, lambda, est.tht, frailty = "LogN") 
  coef2 = rs2$coef
  est.tht = rs2$est.tht
  lambda = rs2$lambda
  
  v_be = coef2 - 2*coef1 + coef
  al_be = sum(u_be*v_be)/sum(v_be^2)
  if (al_be > -1) {al_be = -1}
  
  coef = coef - 2*al_be*u_be + al_be^2*v_be
  
  rs = EMprocess(y, X, d, coef, lambda, est.tht, frailty = "LogN") 
  coef = rs$coef
  est.tht = rs$est.tht
  lambda = rs$lambda
  
  ell[k+1] = logLikihood(y, X, d, coef, lambda, est.tht, frailty = "LogN")
  
  error = abs(ell[k+1]-ell[k])/(1+abs(ell[k]))
  cat(error, '\n')
  k = k+1
}
end = proc.time()[1]
time = end - start
