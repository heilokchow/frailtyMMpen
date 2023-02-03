
# frailtyMM ---------------------------------------------------------------

library("survival")
library("SuppDists")

# Gamma

coef_test = rep(0, 30)
th_test = 0

for (i in 1:100) {
  yy = sample_CL(init.var = 1, cen = 5, frailty = "Gamma")
  
  y = yy$y 
  d = yy$d
  X = yy$X
  
  start = proc.time()[1]
  rs1 = frailtyMM_CL(y, X, d, frailty = "Gamma")
  end = proc.time()[1]
  end - start
  
  coef_test = coef_test + rs1$coef
  th_test = th_test + rs1$est.tht
}

round(rs1$coef, 2)
rs1$est.tht
rs1$likelihood


# InvGauss

yy = sample_CL(init.var = 2, cen = 20, frailty = "LogN")

y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM_CL(y, X, d, frailty = "Gamma")
rs1 = frailtyMM_CL(y, X, d, coef.ini = rs1$coef, est.tht.ini = rs1$est.tht, lambda.ini = rs1$lambda, frailty = "LogN")
rs1 = frailtyMM_CL(y, X, d, frailty = "LogN")
end = proc.time()[1]
end - start

start = proc.time()[1]
rs1 = frailtyMM_CL(y, X, d, coef.ini = rs1$coef, est.tht.ini = rs1$est.tht, lambda.ini = rs1$lambda, frailty = "LogN",
                   penalty = "SCAD", tune = 2)
end = proc.time()[1]
end - start


start = proc.time()[1]
rs2 = frailtyMMpen(y, X, d, type = "Cluster", frailty = "Gamma", penalty = "SCAD")
end = proc.time()[1]
end - start

round(rs1$coef, 2)
rs2$tune.min

plot.fpen(rs2)

start = proc.time()[1]
rs1 = frailtyMM_CL(y, X, d, frailty = "InvGauss", penalty = "MCP", tune = 3)
end = proc.time()[1]
end - start


round(rs1$coef, 2)
rs1$est.tht
rs1$likelihood

# LogN

yy = sample_CL(coef = c(1, 2, 3), a = 30, b = 15, init.var = 1, frailty = "PVF", power = 2.5, cen = 500)
y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM_CL(y, X, d)
end = proc.time()[1]
end - start



start = proc.time()[1]
rs1 = frailtyMM_CL(y, X, d, frailty = "PVF", power = 2.5)
end = proc.time()[1]
end - start

rs1$coef
rs1$est.tht
rs1$likelihood



# frailtyMM ME ------------------------------------------------------------


# LogN

yy = sample_ME(init.var = 1, cen = 0.1, frailty = "InvGauss")
y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM_ME(y, X, d, frailty = "InvGauss")
end = proc.time()[1]
end - start


start = proc.time()[1]
rs1 = frailtyMM_ME(y, X, d, frailty = "InvGauss", penalty = "MCP", tune = 2)
end = proc.time()[1]
end - start


start = proc.time()[1]
rs2 = frailtyMMpen(y, X, d, type = "Multiple", frailty = "InvGauss", penalty = "MCP")
end = proc.time()[1]
end - start

plot.fpen(rs2)

round(rs1$coef, 2)
rs1$est.tht
rs1$likelihood


# InvGauss

yy = sample_ME(init.var = 1, cen = 4, frailty = "InvGauss")
y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM_ME(y, X, d, frailty = "InvGauss", penalty = "LASSO", tune = 0.01)
end = proc.time()[1]
end - start

rs1$coef
rs1$est.tht
rs1$likelihood

# Gamma

yy = sample_ME(init.var = 1, cen = 4, frailty = "Gamma")

y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs1 = frailtyMM_ME(y, X, d, frailty = "Gamma", penalty = "SCAD", tune = 0.1)
end = proc.time()[1]
end - start

rs1$coef
rs1$est.tht
rs1$likelihood


# frailty RE --------------------------------------------------------------


yy = sample_CL(coef = c(1, 2), init.var = 1, a = 50, b = 10, cen = 50, frailty = "LogN")

y1 = yy$y 
d1 = yy$d
X1 = yy$X

y = list()
X = list()
d = list()

for (i in 1:length(y1)) {
  y[[i]] = y1[i]  
  X[[i]] = t(matrix(X1[i,,]))  
  d[[i]] = d1[i, ]  
}

rs = frailtyMM_RE(y, X, d, frailty = "LogN")
rs1 = frailtyMM_CL(y1, X1, d1, coef.ini = c(0.5, 0.5), est.tht.ini = 1, lambda.ini = rep(1/50, 50), frailty = "LogN")
rs1 = frailtyMM_CL(y1, X1, d1, frailty = "LogN")
rs$coef
rs1$coef
rs1$est.tht
rs$est.tht


coef_test = c(0,0)

for (i in 1:100) {
  df = sample_RE(coef = c(1, 2), init.var = 1, n = 50, cen = 200, frailty = "Gamma")

  start = proc.time()[1]
  rs1 = frailtyMM(Surv(start, end, status) ~ . + cluster(id), df, frailty = "Gamma")
  end = proc.time()[1]
  end - start
  
  coef_test = coef_test + rs1$coef
  
  gam <- emfrail(Surv(start, end, status) ~ V1 + V2 + cluster(id), data = df)
  summary(gam)
  
  cat(i, "--------\n")
}

yy = sample_RE(init.var = 1, n = 50, cen = 100, frailty = "InvGauss")

y = yy$y 
d = yy$d
X = yy$X

start = proc.time()[1]
rs2 = frailtyMMpen(y, X, d, type = "Recurrent", frailty = "InvGauss", penalty = "LASSO")
end = proc.time()[1]

plot.fpen(rs2)

# coxph -------------------------------------------------------------------

library(survival)

yy = sample_CL(init.var = 1, cen = 5, frailty = "Gamma")

y = yy$y 
d = yy$d
X = yy$X

a = nrow(y)
b = ncol(y)

newX = matrix(0, nrow = 0, ncol = 30)
id = c()
for (i in 1:a) {
  newX = rbind(newX, X[i,,])
  id = c(id, rep(i, b))
}


vy = as.vector(t(y))
vd = as.vector(t(d))

data = data.frame(id=id,x1=newX[,1],x2=newX[,2],x3=newX[,3],x4=newX[,4],x5=newX[,5],
                  x6=newX[,6],x7=newX[,7],x8=newX[,8],x9=newX[,9],x10=newX[,10],
                  x11=newX[,11],x12=newX[,12],x13=newX[,13],x14=newX[,14],x15=newX[,15],
                  x16=newX[,16],x17=newX[,17],x18=newX[,18],x19=newX[,19],x20=newX[,20],
                  x21=newX[,21],x22=newX[,22],x23=newX[,23],x24=newX[,24],x25=newX[,25],
                  x26=newX[,26],x27=newX[,27],x28=newX[,28],x29=newX[,29],x30=newX[,30],times=vy,status=vd)

aa = coxph(Surv(times,status) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+
             x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+frailty(id,dist='gamma'), data)
aa = coxph(Surv(times,status) ~ x1, data)
summary(aa)
aa$coefficients

start = proc.time()[1]
bb = emfrail(Surv(times,status) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+
               x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30 + cluster(id), data = data)
end = proc.time()[1]
end - start

X1 <- model.matrix(Surv(times,status) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+
                     x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30 + cluster(id), data)

pos_cluster_X1 <- grep("cluster", colnames(X1))
pos_terminal_X1 <- grep("terminal", colnames(X1))
pos_strata_X1 <- grep("strata", colnames(X1))

X2 <- X1[,-c(1, pos_cluster_X1, pos_terminal_X1, pos_strata_X1), drop=FALSE]

Y <- Surv(rep(0, length(vy)), vy, vd)
mcox <- survival::agreg.fit(x = X2, y = Y, strata = NULL, offset = NULL, init = NULL,
                            control = survival::coxph.control(),
                            weights = NULL, method = "breslow", rownames = NULL)  

summary(mcox)
mcox$coefficients

# Real Data ---------------------------------------------------------------

library(frailtyEM)

y = matrix(rats$time, nrow = 100, ncol = 3, byrow = TRUE)
d = matrix(rats$status, nrow = 100, ncol = 3, byrow = TRUE)
x1 = matrix(rats$rx, nrow = 100, ncol = 3, byrow = TRUE)
x2 = matrix(rats$sex, nrow = 100, ncol = 3, byrow = TRUE)
x3 = matrix(0, 100, 3)
for (i in 1:100) {
  for (j in 1:3) {
    if (x2[i, j] == "f") {
      x3[i, j] = 0
    } else {
      x3[i, j] = 1
    }
  }
}


X = array(0, c(100, 3, 2))
X[,,1] = x1
X[,,2] = x3


# CGD
rs1 = frailtyMM(Surv(tstart, tstop, status) ~ sex + treat + cluster(id), cgd, frailty = "InvGauss", tol = 1e-6)
rs1$coef

gam <- emfrail(Surv(tstart, tstop, status) ~ sex + treat + cluster(id), data = cgd, distribution = emfrail_dist(dist = 'pvf'))
summary(gam)


# Kidney
kidney <- kidney[c("time", "status", "id", "age", "sex" )]
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
head(kidney)

rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), kidney, frailty = "Gamma")
rs2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), kidney, frailty = "LogN")
rs1$coef
rs1$est.tht
summary(rs1)

m_gam <- emfrail(Surv(time, status) ~ age + sex + cluster(id), data = kidney, distribution = emfrail_dist(dist = 'gamma'))
summary(m_gam)

f_pack <- frailtyPenal(Surv(time, status) ~ age + sex + cluster(id), data = kidney, n.knots = 14, kappa = 10000)
summary(f_pack)
plot(rs2)


# Rat
head(rats)

rs1 = frailtyMM(Surv(time, status) ~ . + cluster(litter), rats, frailty = "Gamma")
rs1$coef
rs1$est.tht


mm_gam = emfrail(Surv(time, status) ~ rx + sex + cluster(litter), data = rats, distribution = emfrail_dist(dist = 'gamma'))
summary(mm_gam)

f_pack <- frailtyPenal(Surv(time, status) ~ rx + sex + cluster(litter), data = rats, n.knots = 14, kappa = 10000)
summary(f_pack)

rs2 = frailtyMMpen(Surv(time, status) ~ . + cluster(litter), rats, frailty = "LogN")

