
# frailtyMM ---------------------------------------------------------------

library("survival")
library("SuppDists")

# Gamma

yy = sampleF(init.var = 1, cen = 5, frailty = "Gamma")

y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM(y, X, d, frailty = "Gamma")
end = proc.time()[1]
end - start

rs1$coef
rs1$est.tht
rs1$likelihood


# InvGauss

yy = sampleF(init.var = 1, cen = 5, frailty = "InvGauss")

y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM(y, X, d, frailty = "InvGauss")
end = proc.time()[1]
end - start

rs1$coef
rs1$est.tht
rs1$likelihood

# LogN

yy = sampleF(init.var = 1, cen = 4)
y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM(y, X, d)
end = proc.time()[1]
end - start

rs1$coef
rs1$est.tht
rs1$likelihood


# coxph -------------------------------------------------------------------

library(survival)
newX = matrix(0, nrow = 0, ncol = p)
id = c()
for (i in 1:a) {
  newX = rbind(newX, X[i,,])
  id = c(id, rep(i, b))
}


vy = as.vector(t(y))
vd = as.vector(t(d))

data = data.frame(id=id,x1=newX[,1],x2=newX[,2],x3=newX[,3],x4=newX[,4],x5=newX[,5],
                  x6=newX[,6],x7=newX[,7],x8=newX[,8],x9=newX[,9],x10=newX[,10],times=vy,status=vd)

aa = coxph(Surv(times,status) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+frailty(id,dist='gauss'), data)
aa = coxph(Surv(times,status) ~ x1, data)
summary(aa)
aa$coefficients
