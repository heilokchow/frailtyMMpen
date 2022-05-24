library(survival)

n = 10
estMM = matrix(0,n,15)

for (r in 1:n) {
  yy = sampleGamma(seed = r*1000)
  
  a = dim(yy$X)[1]
  b = dim(yy$X)[2]
  p = dim(yy$X)[3]
  latent = matrix(1, a, b)/(a * b)
  
  start = proc.time()[1]
  
  ret = MEGammaFrailty(yy$y, yy$X, yy$d, rep(1, p), latent, 1)
  
  end = proc.time()[1]
  time = end-start
  
  estMM[r,] = c(ret$iteration, time, ret$likelihood, ret$theta, ret$coef, ret$mLambda)
}


aveMM=c(mean(estMM[,1]),mean(estMM[,2]),mean(estMM[,3]),mean(estMM[,4]),mean(estMM[,5]),mean(estMM[,6]),mean(estMM[,7]),mean(estMM[,8]),mean(estMM[,9]),mean(estMM[,10]),mean(estMM[,11]),mean(estMM[,12]),mean(estMM[,13]),mean(estMM[,14]),mean(estMM[,15]))
stdMM=c(sd(estMM[,1]),sd(estMM[,2]),sd(estMM[,3]),sd(estMM[,4]),sd(estMM[,5]),sd(estMM[,6]),sd(estMM[,7]),sd(estMM[,8]),sd(estMM[,9]),sd(estMM[,10]),sd(estMM[,11]),sd(estMM[,12]),sd(estMM[,13]),sd(estMM[,14]),sd(estMM[,15]))
biasMM = aveMM - c(0,0,0,10,-2,3,-4,5,-6,7,6,5,4,3,5)

aveMM
stdMM
biasMM