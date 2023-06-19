
# frailtyMM ---------------------------------------------------------------

library("survival")
library("SuppDists")

# Gamma

coef_test = rep(0, 2)
th_test = 0

for (i in 1:100) {
  yy = sample_CL(coef = c(1,2), init.var = 1, cen = 5, frailty = "Gamma")
  
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

rs1 = frailtyMM_CL(y, X, d, frailty = "Gamma")
rs1 = frailtyMM_CL(y, X, d, coef.ini = rs1$coef, est.tht.ini = rs1$est.tht, lambda.ini = rs1$lambda, frailty = "LogN")
start = proc.time()[1]
rs1 = frailtyMM_CL(y, X, d, frailty = "LogN")
# rs1 = frailtyMM_CL(y, X, d, frailty = "LogN", penalty = "SCAD", tune = 1)
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

set.seed(1)
yy = sample_CL(init.var = 1, cen = 5, frailty = "Gamma")

test = yy[[1]]
colnames(test)[1] = "(ehe)"
simdataCL = yy[[1]]
# rownames(simdataME) = seq_len(nrow(simdataME))
save(simdataCL, file = "data/simdataCL.RData")

y = yy[[2]] 
d = yy[[4]] 
X = yy[[3]] 

rs = frailtyMM_RE(y, X, d, frailty = "Gamma")

p1 = proc.time()[1]
rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), simdataCL, frailty = "Gamma")
rs1 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), simdataCL, frailty = "Gamma", penalty = "SCAD")
p2 = proc.time()[1]

p1 = proc.time()[1]
gam <- emfrail(Surv(time, status) ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + 
                 V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + 
                 V21 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + cluster(id), data = simdataCL)
p2 = proc.time()[1]

p1 = proc.time()[1]
f_pack <- frailtyPenal(Surv(time, status) ~ V1 + V2 + cluster(id), data = simdataCL, n.knots = 14, kappa = 10000)
p2 = proc.time()[1]
p2 - p1

rs1 = frailtyMM_CL(y1, X1, d1, coef.ini = c(0.5, 0.5), est.tht.ini = 1, lambda.ini = rep(1/50, 50), frailty = "LogN")
rs1 = frailtyMM_CL(y1, X1, d1, frailty = "LogN")
rs$coef
rs1$coef
rs1$est.tht
rs$est.tht


coef_test = c(0,0)
th_test = 0

for (i in 1:100) {
  set.seed(i)
  df = sample_CL(coef = c(1, 2), init.var = 4, cen = 20, frailty = "InvGauss")

  # start = proc.time()[1]
  rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), df[[1]], frailty = "Gamma")
  # end = proc.time()[1]
  # end - start
  
  
  # gam <- emfrail(Surv(start, end, status) ~ V1 + V2 + cluster(id), data = df[[1]])
  coef_test = coef_test + rs1$coef
  th_test = th_test + rs1$est.tht
  # summary(gam)
  
  cat(i, "--------\n")
}

rs2 = emfrail(Surv(time, status) ~ V1 + V2 + cluster(id), df[[1]], distribution = emfrail_dist(dist = 'gamma'))
rs3 = frailtyPenal(Surv(time, status) ~ V1 + V2 + cluster(id), data = df[[1]], n.knots = 14, kappa = 10000)

test = predict(rs1, df[[1]][1,], surv = TRUE)
test = predict(rs1)
plot(rs2, type = "pred", newdata = df[[1]][1,])
predict(rs2, type = "marginal", newdata = df[[1]][1,], conf_int = "regular")
predict(rs2, lp = 0, type = "marginal", conf_int = "regular")

plot(rs3, type.plot = "Hazard", conf.bands=TRUE)

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
rs1 = frailtyMM(Surv(tstart, tstop, status) ~ sex + treat + cluster(id), cgd, frailty = "Gamma", tol = 1e-6)
rs1 = frailtyMMpen(Surv(tstart, tstop, status) ~ sex + treat + cluster(id), cgd, frailty = "LogN", penalty = "SCAD", tol = 1e-6)
rs1$coef
summary(rs1)

gam <- emfrail(Surv(tstart, tstop, status) ~ sex + treat + cluster(id), data = cgd, distribution = emfrail_dist(dist = 'gamma'))
summary(gam)


# Kidney
kidney <- kidney[c("time", "status", "id", "age", "sex" )]
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
head(kidney)

kidney2 = kidney[seq(1,75,2),]

p1 = proc.time()[1]
rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), kidney, frailty = "Gamma")
rs1 = frailtyMMpen(Surv(time, status) ~ age + sex + cluster(id), kidney, frailty = "Gamma")
rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), kidney1, frailty = "InvGauss")
# rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), kidney, frailty = "LogN")
# rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), kidney, frailty = "Gamma")
# rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), kidney, frailty = "PVF", power = 1.5)
rs2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), kidney, frailty = "InvGauss", penalty = "SCAD", tune = 0.2)
rs2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), kidney, frailty = "InvGauss", penalty = "SCAD")
p2 = proc.time()[1]
p2 - p1

rs1$coef
rs1$est.tht
summary(rs1)
plot(rs2)

kidney1 = kidney[c(-3,-5),]

p3 = proc.time()[1]
m_gam <- emfrail(Surv(time, status) ~ age + sex + cluster(id), data = kidney1, distribution = emfrail_dist(dist = 'pvf', pvfm = -0.5))
m_gam <- emfrail(Surv(time, status) ~ age + sex + cluster(id), data = kidney, distribution = emfrail_dist(dist = 'gamma'))
m_gam <- emfrail(Surv(time, status) ~ age + sex + cluster(id), data = kidney, distribution = emfrail_dist(dist = 'pvf', pvfm = -0.5))
p4 = proc.time()[1]
p4 - p3

summary(mm_gam)

p5 = proc.time()[1]
f_pack <- frailtyPenal(Surv(time, status) ~ age + sex + cluster(id), data = kidney, n.knots = 14, kappa = 10000)
f_pack <- frailtyPenal(Surv(time, status) ~ age + sex + cluster(id), data = kidney1, n.knots = 14, kappa = 10000)
p6 = proc.time()[1]

summary(f_pack)
plot(f_pack, type.plot = "Su", conf.bands=TRUE)
prediction(f_pack, kidney, t=seq(10,20,by=1), window=5)

# DEBUG -------------------------------------------------------------------

y = rs2$y
X = rs2$X
d = rs2$d

res = frailtyMM_CL(y, X, d, frailty = "LogN", penalty = "SCAD", tune = 3)


set.seed(5)
sdata = sample_CL(init.var = 1, cen = 5, frailty = "LogN", a = 50, b = 10)

test_coef = rep(0, 30)
test_th = 0
for (i in 1:100) {
  set.seed(i)
  sdata = sample_ME(init.var = 1, cen = 5, frailty = "LogN", n = 100)
  sdata = sample_CL(init.var = 1, cen = 5, frailty = "LogN", a = 50, b = 10)
  # y = sdata$y
  # X = sdata$X
  # d = sdata$d
  # # 
  # y = sdata$data$time
  # X = unname(unlist(sdata$data[,1:30]))
  # d = sdata$data$status
  
  # rs0 = frailtyMM_ME(y, X, d, frailty = "LogN")
  
  p3 = proc.time()[1]
  rs1 = frailtyMM(Surv(time, status) ~ . + event(id), sdata$data, frailty = "LogN", tol = 1e-6)
  rs1 = frailtyMM(Surv(time, status) ~ . + cluster(id), sdata$data, frailty = "LogN")
  rs1 = frailtyMMpen(Surv(time, status) ~ . + event(id), sdata$data, frailty = "LogN", tol = 1e-5, power = 1.5, penalty = "LASSO", maxit = 400)
  p4 = proc.time()[1]
  p4 - p3
  
  test_coef = test_coef + rs1$coef
  test_th = test_th + rs1$est.tht
}

summary(m_gam)

p1 = proc.time()[1]
m_gam <- emfrail(Surv(time, status) ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + 
                   V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + 
                   V21 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + cluster(id), 
                 data = sdata, distribution = emfrail_dist(dist = 'pvf', pvfm = -0.5))
p2 = proc.time()[1]
p2 - p1


p5 = proc.time()[1]
f_pack <- frailtyPenal(Surv(time, status) ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + 
                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + 
                         V21 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + cluster(id),
                       data = sdata, n.knots = 14, kappa = 10000)
p6 = proc.time()[1]

summary(f_pack)

# Rat
head(rats)

p1 = proc.time()[1]
rs1 = frailtyMM(Surv(time, status) ~ . + cluster(litter), rats, frailty = "InvGauss", tol = 1e-6)
p2 = proc.time()[1]
rs1$coef
rs1$est.tht
summary(rs1)

mm_gam = emfrail(Surv(time, status) ~ rx + sex + cluster(litter), data = rats, distribution = emfrail_dist(dist = 'pvf', pvfm = -0.5))
summary(mm_gam)

f_pack <- frailtyPenal(Surv(time, status) ~ rx + sex + cluster(litter), data = rats, n.knots = 14, kappa = 10000)
summary(f_pack)

rs2 = frailtyMMpen(Surv(time, status) ~ . + cluster(litter), rats, frailty = "Gamma", penalty = "SCAD")
plot(rs2)

# Glmnet

library(glmnet)
data(QuickStartExample)
x <- QuickStartExample$x
y <- QuickStartExample$y
fit <- glmnet(x, y)
plot(fit)
fit
summary(fit)


# Transfer to GSL ---------------------------------------------------------

set.seed(5)
yy = sample_RE(init.var = 1, cen = 5, frailty = "Gamma")

sdata = simdataCL
y = sdata$time
d = sdata$status
X = sdata[,1:30]

X = unname(unlist(sdata[,1:30]))
X = matrix(X, ncol = 30)
N = nrow(sdata)
lambda = c(rep(1/N, N))

y1 = yy[[2]]
X1 = yy[[3]]
d1 = yy[[4]]
lambda1 = yy[[5]]

y = yy[[2]]
X = yy[[3]]
d = yy[[4]]
lambda = yy[[5]]

vy = as.vector(y)
vd = as.vector(d)

a = 50
p = 30
coef = rep(1/p, p)

est.tht = 1
frailty = "LogN"
power = 1.5
penalty = NULL

y1 = vy
d1 = vd
X1 = matrix(as.vector(X), nrow = N, ncol = p)
newid = sdata$id-1

coef0 = coef
est.tht0 = 1
lambda0 = lambda

formula = Surv(time, status) ~ . + cluster(id)
rs1 = frailtyMM(Surv(time, status) ~ . + event(id), sdata, frailty = "InvGauss", tol = 1e-6)
rs1 = frailtyMM(Surv(start, end, status) ~ . + cluster(id), sdata, frailty = "LogN", tol = 1e-6, maxit = 100)
rs1 = frailtyMM(Surv(start, end, status) ~ . + cluster(id), sdata, frailty = "Gamma", tol = 1e-6, maxit = 200)

gam_cl2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), simdataCL, frailty = "Gamma", penalty = "SCAD", gam = 4)
gam_cl1 = frailtyMM(Surv(time, status) ~ . + cluster(id), simdataCL, frailty = "Gamma")

predict(gam_cl1)

rs2 = emfrail(Surv(start, end, status) ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + 
                V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + 
                V21 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + cluster(id), data = sdata, distribution = emfrail_dist(dist = 'gamma'))

ttemp = c()
for (i in 1:93) { 
  ttemp = c(ttemp, test1[[7]][3+(i-1)*a])
}

f1 <- function() {
  # test = MMCL_TEST(y1, X1, d1, coef, lambda, 1.0, 1, 0, 0.1, id, N, a, p, power)
  # test = MMCL_TEST(y, X, d, coef, lambda, 1.0, 1, 0, 0.1, id, N, a, p, power)
  # p1 = proc.time()[1]
  test1 = MMRE_TEST(y, X, d, coef, lambda, 1, 1, 0, 0, id, N, a, p, power, 1)
  # p2 = proc.time()[1]
  # MMME_TEST(y, X, d, coef, lambda, tht, frailty, penalty, tune, id, N, a, p, power, type)
}

f2 <- function(){
  test2 = MMprocess_RE(y1, X1, d1, coef, lambda1, 1)
}


microbenchmark(f1(), f2())

p1 = proc.time()[1]
f1()
p2 = proc.time()[1]
p2-p1


f1 <- function(x) {
  return (exp(ldTweedie(x,mu=1,p=2,phi=10)[1]))
}

integrate(Vectorize(f1), 0, Inf)

for (i in 1:20) {
  cat(ldTweedie(6,1,1.5,i)[2], "\n")
}

log(dtweedie(1,mu=1,p=2,phi=100))
f <- function(x) {
  ldTweedie(10,mu=1,p=2,phi=x)[1]
}

ldTweedie(0.001,mu=1,p=2,phi=100)
log(dgamma(0.001, 1/100, 1/100))

numDeriv::grad(f, 1)

dtweedie.dldphi.saddle(y=2,mu=1,phi=2,power=1.5)
dtweedie.dldphi(y=2,mu=1,phi=2,power=1.5)
ldTweedie(2,mu=1,p=1.5,phi=0.3)


# Example -----------------------------------------------------------------

# Kidney data fitted by Clustered Inverse Gaussian Frailty Model
InvG_real_cl = frailtyMM(Surv(time, status) ~ age + sex + cluster(id), kidney, frailty = "InvGauss")
InvG_real_cl

# Cgd data fitted by Recurrent Log-Normal Frailty Model

logN_real_re = frailtyMM(Surv(tstart, tstop, status) ~ sex + treat + cluster(id), cgd, frailty = "Gamma")
logN_real_re

gam_re = frailtyMM(Surv(tstart, tstop, status) ~  sex + treat + cluster(id), cgd, frailty = "Gamma")

# Plot the survival curve based on baseline hazard
plot(gam_re, surv = TRUE)

# Construct new data and plot the cumulative hazard based on new data
newre = c(1, 1, 2)
names(newre) = c(gam_re$coefname, "id")
plot(gam_re, newdata = newre)

# Simulated data example

data(simdataCL)
data(simdataME)
data(simdataRE)

# Parameter estimation under different model structure and frailties

# Clustered Gamma Frailty Model
gam_cl = frailtyMM(Surv(time, status) ~ . + cluster(id),
                   simdataCL, frailty = "Gamma")

gam_cl = frailtyMMpen(Surv(time, status) ~  . + cluster(id), simdataCL, frailty = "Gamma")

# Calculate the survival curve based on baseline hazard
predict(gam_cl, surv = TRUE)

# Construct new data and calculate the cumulative hazard based on new data
newcl = c(gam_cl$X[1,], 2)
names(newcl) = c(gam_cl$coefname, "id")
predict(gam_cl, newdata = newcl)

# Clustered Log-Normal Frailty Model
logn_cl = frailtyMM(Surv(time, status) ~ . + cluster(id),
                    simdataCL, frailty = "LogN")

# Clustered Inverse Gaussian Frailty Model
invg_cl = frailtyMM(Surv(time, status) ~ . + cluster(id),
                    simdataCL, frailty = "InvGauss")

# Multi-event Gamma Frailty Model
gam_me = frailtyMM(Surv(time, status) ~ . + cluster(id),
                   simdataCL, frailty = "Gamma")

# Recurrent event Gamma Frailty Model
gam_re = frailtyMM(Surv(start, end, status) ~ . + cluster(id),
                   simdataRE, frailty = "Gamma")

# Recurrent event Log-Normal Frailty Model
logn_re = frailtyMM(Surv(start, end, status) ~ . + cluster(id),
                   simdataRE, frailty = "LogN")

# Recurrent event Inverse Gaussian Frailty Model
invg_re = frailtyMM(Surv(start, end, status) ~ . + cluster(id),
                    simdataRE, frailty = "InvGauss")


# Obtain the summary statistics under fitted model

coef(gam_cl)
summary(gam_cl)

# Penalized regression under clustered frailty model

# Clustered Gamma Frailty Model

gam_cl1 = frailtyMM(Surv(time, status) ~ . + cluster(id),
                       simdataCL, frailty = "Gamma")

# Using default tuning parameter sequence
gam_cl1 = frailtyMMpen(Surv(time, status) ~ . + cluster(id),
                       simdataCL, frailty = "Gamma")

# Using given tuning parameter sequence
gam_cl2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id),
                       simdataCL, frailty = "Gamma", tune = 0.1)

# Obtain the coefficient where minimum BIC is obtained
coef(gam_cl1)

# Obtain the coefficient with tune = a.
coef(gam_cl1, tune = 0.2)

# Plot the regularization path
plot(gam_cl1)

# Get the degree of freedom and BIC for the sequence of tuning parameters provided
print(gam_cl1)



# Paper writing -----------------------------------------------------------

# HD data generate

hddata = sample_CL(coef = c(rep(4, 5), rep(0, 95)), a = 100, b = 2, init.var = 1, cen = 3, frailty = "LogN")
hdata1 = hddata[[1]]
sum(hddata$d)

rs1 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), hdata1, frailty = "LogN", penalty = "SCAD", maxit = 100)
plot(rs1)
tune = exp(seq(-1, 1.5, 0.125))
tune = seq(0.2, 3, 0.2)

library(DLBCL)

data(exprLym)
test = exprs(exprLym)
test1 = pData(exprLym)

test2 = cbind(test1, t(test)[,1:10])
test2 = test2[1:190,]

rm = which(test2$FollowUpYears < 0.3)
test2 = test2[-rm, ]
colnames(test2)[6:15] = sapply(seq(1,10,1), function (x) {paste0("x", x, sep = "")})


rs = frailtyMMpen(Surv(FollowUpYears, Status) ~ ., data = test2, tune = exp(seq(-0.1, 2, 0.1)), frailty = "Gamma")
rs = frailtyMM(Surv(FollowUpYears, Status) ~ x5, data = test2, frailty = "Gamma")

plot(test2$x3)
plot(log(test2$FollowUpYears))

cor(test2$x3, test2$x1)

data(pbc, package="survival")
pbc_data = na.omit(pbc)
death = which(pbc_data$status == 2)
trans = which(pbc_data$status == 1)
pbc_data$status[death] = 1
pbc_data = pbc_data[-trans,]
pbc_data = pbc_data[,-1]
pbc_data$trt = pbc_data$trt - 1

gamfrail = frailtyMM(Surv(time, status) ~ ., data = pbc_data, frailty = "Gamma")
summary(gamfrail)

gamfrailLASSO = frailtyMMpen(Surv(time, status) ~ ., data = pbc_data, frailty = "Gamma", penalty = "LASSO")
gamfrailSCAD = frailtyMMpen(Surv(time, status) ~ ., data = pbc_data, frailty = "Gamma", penalty = "SCAD")
gamfrailMCP = frailtyMMpen(Surv(time, status) ~ ., data = pbc_data, frailty = "Gamma", penalty = "MCP")

plot(gamfrailLASSO)
plot(gamfrailSCAD)
plot(gamfrailMCP)

ecoef = round(data.frame(coef(gamfrail), coef(gamfrailLASSO), coef(gamfrailSCAD), coef(gamfrailMCP)), 3)
colnames(ecoef) = c("ori", "lasso", "scad", "mcp")
ecoef

round(cbind(coef(gamfrail), coef(gamfrailLASSO), coef(gamfrailSCAD), coef(gamfrailMCP)), 3)

coef(gamfrailLASSO)
rowSums(pbc)

data("simdataME")
logNfrail = frailtyMM(Surv(time, status) ~ . + event(id), data = simdataME, frailty = "LogN", maxit = 1000)

data("simdataRE")
InvGfrail = frailtyMM(Surv(start, end, status) ~ . + cluster(id), data = simdataRE, frailty = "InvGauss", maxit = 1000)


# Comparison --------------------------------------------------------------

kk = seq(1, 96, 10)
kk = c(1, 31)

MSEall = matrix(0, nrow = 5, ncol = 10)
Biaall = matrix(0, nrow = 5, ncol = 10)

MSEvar = matrix(0, nrow = 5, ncol = 10)
Biavar = matrix(0, nrow = 5, ncol = 10)


Timeall = matrix(0, nrow = 5, ncol = 10)
convall = matrix(0, nrow = 5, ncol = 10)

run1 = run2 = run3 = run4 = run5 = 1

for (j in 2) {
  for (i in 1:100) {
    coefNew = matrix(c(1, 2, 3, 4, rep(0, kk[j])))
    # coefNew = matrix(c(rep(-2, 10), rep(3, 10)))
    
    p = length(coefNew)
    set.seed(i)
    # DK = sample_CL(coef = coefNew, a = 50, b = 5, cen = 10, frailty = "InvGauss", init.var = 1)
    # DK = sample_RE(coef = coefNew, init.var = 1, cen = 100, frailty = "LogN")
    DK = sample_ME(coef = coefNew, init.var = 0.4812, cen = 0.1, frailty = "LogN", n = 200)
    # DK = sample_ME(coef = coefNew, init.var = 1, cen = 0.1, frailty = "Gamma", n = 50)
    n1 = nrow(DK$data)/2
    DK$data$pid = rep(seq(1, n1/2, 1), 2)
    
    st = "V1"
    for (z in 2:p) {
      st = paste0(st, "+V", z)
    }
    
    # form = formula(paste("Surv(time, status) ~ ", st, " + cluster(id)"))
    form = formula(paste("Surv(time, status) ~ ", st, " + event(id)"))
    # form = formula(paste("Surv(start, end, status) ~ ", st, " + cluster(id)"))
    
    if (run1 > 0) {
      p1 = Sys.time()
      er = tryCatch(expr = {
        md = frailtyMM(form, data = DK$data, frailty = "LogN", maxit = 200)
      }, error = function(e){
        return(1)
      })
      p2 = Sys.time()

      ms = sum((md$coef - coefNew)^2)/p

      if (ms > 10) {
        convall[1, j] = convall[1, j] + 1
      } else {
        MSEall[1, j] = MSEall[1, j] + ms
        Biaall[1, j] = Biaall[1, j] + sum(md$coef - coefNew)/p
        MSEvar[1, j] = MSEvar[1, j] + (md$est.tht - 1)^2
        Biavar[1, j] = Biavar[1, j] + md$est.tht - 1
        diff = difftime(p2, p1, units = "secs")
        Timeall[1, j] = Timeall[1, j] + diff
        if (diff > 120) {
          run1 = run1 - 1
        }
      }

    } else {
      Timeall[1, j] = Timeall[1, j] + 120
    }
    
    cat(i, j, ms, "1\n")

    # if (i <= 100 && run2 > 0) {
    #   p1 = Sys.time()
    #   er = tryCatch(expr = {
    #     md = emfrail(form, data = DK$data, emfrail_dist(dist = "pvf", pvfm = -0.5))
    #   }, error = function(e){
    #     return(1)
    #   })
    #   p2 = Sys.time()
    # 
    #   diff = difftime(p2, p1, units = "secs")
    # 
    #   if (diff > 120) {
    #     run2 = run2 - 1
    #   }
    # 
    #   if (length(er) == 1 && er == 1) {
    #     convall[2, j] = convall[2, j] + 1
    #   } else {
    # 
    #     ms = sum((md$coefficients - coefNew)^2)/p
    #     MSEall[2, j] = MSEall[2, j] + ms
    #     Timeall[2, j] = Timeall[2, j] + diff
    #     MSEvar[2, j] = MSEvar[2, j] + (exp(md$logtheta) - 1)^2
    #     Biavar[2, j] = Biavar[2, j] + exp(md$logtheta) - 1
    #     Biaall[2, j] = Biaall[2, j] + sum(md$coefficients - coefNew)/p
    #   }
    # } else {
    #   convall[2, j] = convall[2, j] + 1
    # }
    
    # cat(i, j, "2\n")

    if (run3 > 0) {

      form2 = formula(paste("Surv(time, status) ~ ", st, " + frailty.gaussian(id)"))
      # form2 = formula(paste("Surv(start, end, status) ~ ", st, " + frailty.gaussian(id)"))
      form2 = formula(paste("Surv(time, status) ~ ", st, " + frailty.gaussian(pid) + strata(id)"))

      p1 = Sys.time()
      er = tryCatch(expr = {
        md = coxph(form2, data = DK$data)
      }, error = function(e){
        return(1)
      })
      p2 = Sys.time()

      diff = difftime(p2, p1, units = "secs")
      if (diff > 120) {
        run3 = run3 - 1
      }

      if ((length(er) == 1 && er == 1) || !prod(!is.na(md$coefficients))) {
        convall[3, j] = convall[3, j] + 1
      } else {
        ms = sum((md$coefficients - coefNew)^2)/p
        MSEall[3, j] = MSEall[3, j] + ms
        Biaall[3, j] = Biaall[3, j] + sum(md$coefficients - coefNew)/p
        MSEvar[3, j] = MSEvar[3, j] + (md[[17]]$`frailty.gaussian(pid)`$theta - 1)^2
        # MSEvar[3, j] = MSEvar[3, j] + (md[[17]]$`frailty(pid)`$theta - 1)^2
        Biavar[3, j] = Biavar[3, j] + md[[17]]$`frailty.gaussian(pid)`$theta - 1
        # Biavar[3, j] = Biavar[3, j] + md[[17]]$`frailty(pid)`$theta - 1
        Timeall[3, j] = Timeall[3, j] + diff

      }
    } else {
      Timeall[3, j] = Timeall[3, j] + 120
    }

    cat(i, j, ms, "3\n")
    #
    # if (run4 > 0) {
    #
    #   p1 = Sys.time()
    #   er = tryCatch(expr = {
    #     md = fitfrail(form, DK$data, frailty = "invgauss")
    #   }, error = function(e){
    #     return(1)
    #   })
    #   p2 = Sys.time()
    #
    #   diff = difftime(p2, p1, units = "secs")
    #   if (diff > 120) {
    #     run3 = run3 - 1
    #   }
    #
    #   if ((length(er) == 1 && er == 1) || is.null(md$beta)) {
    #     convall[4, j] = convall[4, j] + 1
    #   } else {
    #     ms = sum((md$beta - coefNew)^2)/p
    #     MSEall[4, j] = MSEall[4, j] + ms
    #     Biaall[4, j] = Biaall[4, j] + sum(md$beta - coefNew)/p
    #     MSEvar[4, j] = MSEvar[4, j] + (md$theta - 1)^2
    #     Biavar[4, j] = Biavar[4, j] + md$theta - 1
    #     Timeall[4, j] = Timeall[4, j] + diff
    #
    #   }
    # } else {
    #   Timeall[4, j] = Timeall[4, j] + 120
    # }

    # cat(i, j, "4\n")
    #
    if (i <= 100 && run5 > 0) {
      # form1 = formula(paste("Surv(time, status) ~ ", st, " + (1|id)"))
      form1 = formula(paste("Surv(time, status) ~ ", st, " + (1|pid) + (1|id)"))
      p1 = Sys.time()
      er = tryCatch(expr = {
        md = frailtyHL(form1, data = DK$data, mord = 0, dord = 1, RandDist = "Normal")
        }, error = function(e){
        return(1)
        })
      p2 = Sys.time()
      diff = difftime(p2, p1, units = "secs")

      if (diff > 120) {
        run4 = run4 - 1
      }

      if ((length(er) == 1 && er == 1)) {
        convall[5, j] = convall[5, j] + 1
      } else {
        ms = sum((md$FixCoef[,1] - coefNew)^2)/p
        MSEall[5, j] = MSEall[5, j] + ms
        Biaall[5, j] = Biaall[5, j] + sum(md$FixCoef[,1] - coefNew)/p
        MSEvar[5, j] = MSEvar[5, j] + (md$RandCoef[1] - 1)^2
        Biavar[5, j] = Biavar[5, j] + md$RandCoef[1] - 1
        Timeall[5, j] = Timeall[5, j] + diff
      }
    } else {
      convall[5, j] = convall[5, j] + 1
    }
    cat(i, j, ms, "5\n")
  }
  
  
}

Gamma_Cluster = list(MSEall = MSEall,
                     Biaall = Biaall,
                     MSEvar = MSEvar,
                     Biavar = Biavar,
                     Timeall = Timeall,
                     convall = convall)

Gamma_Cluster = list(MSEall = MSEall1,
                     Biaall = Biaall1,
                     MSEvar = MSEvar1,
                     Biavar = Biavar1,
                     Timeall = Timeall1,
                     convall = convall1)

LogN_Cluster = list(MSEall = MSEall,
                     Biaall = Biaall,
                     MSEvar = MSEvar,
                     Biavar = Biavar,
                     Timeall = Timeall,
                     convall = convall)

Gamma_Recurrent = list(MSEall = MSEall,
                    Biaall = Biaall,
                    MSEvar = MSEvar,
                    Biavar = Biavar,
                    Timeall = Timeall,
                    convall = convall)

Inverse_Cluster = list(MSEall = MSEall,
                       Biaall = Biaall,
                       MSEvar = MSEvar,
                       Biavar = Biavar,
                       Timeall = Timeall,
                       convall = convall)

Inverse_Recurrent = list(MSEall = MSEall,
                       Biaall = Biaall,
                       MSEvar = MSEvar,
                       Biavar = Biavar,
                       Timeall = Timeall,
                       convall = convall)

Normal_Recurrent = list(MSEall = MSEall,
                         Biaall = Biaall,
                         MSEvar = MSEvar,
                         Biavar = Biavar,
                         Timeall = Timeall,
                         convall = convall)

# Gamma Recurrent

MSEall = Gamma_Recurrent$MSEall[,1:7]
Biaall = Gamma_Recurrent$Biaall[,1:7]
MSEvar = Gamma_Recurrent$MSEvar[,1:7]
Biavar = Gamma_Recurrent$Biavar[,1:7]
Timeall = Gamma_Recurrent$Timeall[,1:7]
convall = Gamma_Recurrent$convall[,1:7]

dplot = data.frame(x = rep(seq(5, 66, 10), 3),
                   MSE = c(MSEall[1,], MSEall[2,], MSEall[3,]),
                   Bias = c(Biaall[1,], Biaall[2,], Biaall[3,]),
                   MSEtheta = c(MSEvar[1,], MSEvar[2,], MSEvar[3,]),
                   Biastheta = c(Biavar[1,], Biavar[2,], Biavar[3,]),
                   Timeall = c(Timeall[1,], Timeall[2,], Timeall[3,]),
                   convall = c(convall[1,], convall[2,], convall[3,]),
                   group = c(rep("frailtyMMpen", 7), rep("frailtyEM", 7),
                             rep("survival", 7)))

q1 = ggplot(dplot, aes(x = x, y = MSE/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "red", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("MSE(", hat(italic(beta)), ")"))) +
  theme_bw()

q2 = ggplot(dplot, aes(x = x, y = Bias/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "red", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("BIAS(", hat(italic(beta)), ")"))) +
  theme_bw()

q3 = ggplot(dplot, aes(x = x, y = MSEtheta/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "red", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("MSE(", hat(italic(theta)), ")"))) +
  theme_bw()

q4 = ggplot(dplot, aes(x = x, y = Biastheta/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "red", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("BIAS(", hat(italic(theta)), ")"))) +
  theme_bw()

q5 = ggplot(dplot, aes(x = x, y = Timeall/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "red", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("Time(sec)"))) +
  theme_bw() +
  theme(legend.position = "none") 

ggarrange(q1, q2, q3, q4, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")


# LogN Cluster

MSEall = LogN_Cluster$MSEall
Biaall = LogN_Cluster$Biaall
MSEvar = LogN_Cluster$MSEvar
Biavar = LogN_Cluster$Biavar
Timeall = LogN_Cluster$Timeall
convall = LogN_Cluster$convall

dplot = data.frame(x = rep(seq(5, 96, 10), 4),
                   MSE = c(MSEall[1,], MSEall[3,], MSEall[4,], MSEall[5,]),
                   Bias = c(Biaall[1,], Biaall[3,], Biaall[4,], Biaall[5,]),
                   MSEtheta = c(MSEvar[1,], MSEvar[3,], MSEvar[4,], MSEvar[5,]),
                   Biastheta = c(Biavar[1,], Biavar[3,], Biavar[4,], Biavar[5,]),
                   Timeall = c(Timeall[1,], Timeall[3,], Timeall[4,], Timeall[5,]),
                   convall = c(convall[1,], convall[3,], convall[4,], convall[5,]),
                   group = c(rep("frailtyMMpen", 10), rep("survival", 10),
                             rep("frailtySurv", 10), rep("frailtyHL", 10)))

q1 = ggplot(dplot, aes(x = x, y = MSE/(100 - convall), group = group, color = group, fill = group)) +
            geom_line(linewidth = 1) +
            scale_color_manual(values = c("#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
            scale_x_continuous(breaks=(seq(5, 96, 20))) + 
            xlab(expression(italic("p"))) +
            ylab(expression(paste("MSE(", hat(italic(beta)), ")"))) +
            theme_bw()

q2 = ggplot(dplot, aes(x = x, y = Bias/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("BIAS(", hat(italic(beta)), ")"))) +
  theme_bw()

q3 = ggplot(dplot, aes(x = x, y = MSEtheta/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("MSE(", hat(italic(theta)), ")"))) +
  theme_bw()

q4 = ggplot(dplot, aes(x = x, y = Biastheta/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("BIAS(", hat(italic(theta)), ")"))) +
  theme_bw()

q5 = ggplot(dplot, aes(x = x, y = Timeall/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("Time(sec)"))) +
  theme_bw() +
  theme(legend.position = "none") 

ggarrange(q1, q2, q3, q4, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")


# Gamma Cluster

MSEall1 = Gamma_Cluster$MSEall[,1:7]
Biaall1 = Gamma_Cluster$Biaall[,1:7]
MSEvar1 = Gamma_Cluster$MSEvar[,1:7]
Biavar1 = Gamma_Cluster$Biavar[,1:7]
Timeall1 = Gamma_Cluster$Timeall[,1:7]
convall1 = Gamma_Cluster$convall[,1:7]

MSEall1[c(2,5), ] = MSEall[c(2,5), 1:7]
Biaall1[c(2,5), ] = Biaall[c(2,5), 1:7]
MSEvar1[c(2,5), ] = MSEvar[c(2,5), 1:7]
Biavar1[c(2,5), ] = Biavar[c(2,5), 1:7]
Timeall1[c(2,5), ] = Timeall[c(2,5), 1:7]
convall1[c(2,5), ] = convall[c(2,5), 1:7]

MSEall = Gamma_Cluster$MSEall[,1:7]
Biaall = Gamma_Cluster$Biaall[,1:7]
MSEvar = Gamma_Cluster$MSEvar[,1:7]
Biavar = Gamma_Cluster$Biavar[,1:7]
Timeall = Gamma_Cluster$Timeall[,1:7]
convall = Gamma_Cluster$convall[,1:7]
Timeall[5,3] = (Timeall[5,2] + Timeall[5,4])/2

dplot = data.frame(x = rep(seq(5, 66, 10), 5),
                   MSE = c(MSEall[1,], MSEall[2,], MSEall[3,], MSEall[4,], MSEall[5,]),
                   Bias = c(Biaall[1,], Biaall[2,], Biaall[3,], Biaall[4,], Biaall[5,]),
                   MSEtheta = c(MSEvar[1,], MSEvar[2,], MSEvar[3,], MSEvar[4,], MSEvar[5,]),
                   Biastheta = c(Biavar[1,], Biavar[2,],Biavar[3,], Biavar[4,], Biavar[5,]),
                   Timeall = c(Timeall[1,], Timeall[2,], Timeall[3,], Timeall[4,], Timeall[5,]),
                   convall = c(convall[1,], convall[2,], convall[3,], convall[4,], convall[5,]),
                   group = c(rep("frailtyMMpen", 7), rep("frailtyEM", 7), rep("survival", 7),
                             rep("frailtySurv", 7), rep("frailtyHL", 7)))

q1 = ggplot(dplot, aes(x = x, y = MSE/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("MSE(", hat(italic(beta)), ")"))) +
  theme_bw()

q2 = ggplot(dplot, aes(x = x, y = Bias/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("BIAS(", hat(italic(beta)), ")"))) +
  theme_bw()

q3 = ggplot(dplot, aes(x = x, y = MSEtheta/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("MSE(", hat(italic(theta)), ")"))) +
  theme_bw()

q4 = ggplot(dplot, aes(x = x, y = Biastheta/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("BIAS(", hat(italic(theta)), ")"))) +
  theme_bw()

q5 = ggplot(dplot, aes(x = x, y = Timeall/(100 - convall), group = group, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("pink", "#56B4E9", "red", "#F0E442", "#009E73"), name=NULL) + 
  scale_x_continuous(breaks=(seq(5, 96, 20))) + 
  xlab(expression(italic("p"))) +
  ylab(expression(paste("Time(sec)"))) +
  theme_bw() +
  theme(legend.position = "none") 

ggarrange(q1, q2, q3, q4, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")