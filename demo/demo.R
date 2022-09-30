library(survival)

n = 10
estMM = matrix(0,n,16)

for (r in 1:n) {
  yy = sampleGamma(th = 0.5, a = 500, seed = r*1000)
  
  a = dim(yy$X)[1]
  b = dim(yy$X)[2]
  p = dim(yy$X)[3]
  lambda = matrix(1, a, b)/(a * b)
  
  start = proc.time()[1]
  
  y = yy$y
  X = yy$X
  d = yy$d
  alp = 0.1
  th = 0.1
  coef = rep(1, p)
  ret = MEGammaFrailty(y, X, d, coef, lambda, alp, th)
  
  end = proc.time()[1]
  time = end-start
  
  estMM[r,] = c(ret$iteration, time, ret$likelihood, ret$alp, ret$th, ret$coef, ret$mLambda)
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
