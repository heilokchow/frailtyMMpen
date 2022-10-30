int_tao <- function(x, i, est.tht, A, B, D, tao0 = 0, frailty = "LogN", mode = 0) {
  
  ta = 0
  if (frailty == "LogN") {
    si2 = est.tht
    if (mode == 0) {
      ta = 1/(x*sqrt(2*pi*si2))*exp(-(log(x))^2/(2*si2))*x^D[i]*B[i]*exp(-A[i]*x)
    }
    if (mode == 1) {
      ta = 1/(x*sqrt(2*pi*si2))*exp(-(log(x))^2/(2*si2))*x^D[i]*B[i]*exp(-A[i]*x)*x/tao0
    }
    if (mode == 2) {
      ta = (log(x))^2*1/(x*sqrt(2*pi*si2))*exp(-(log(x))^2/(2*si2))*x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
  }

  if (frailty == "InvGauss") {
    mth = est.tht
    if (mode == 0) {
      ta = mth/(sqrt(2*pi)*x^(3/2))*exp(-(x-mth)^2/(2*x))*x^D[i]*B[i]*exp(-A[i]*x)
    }
    if (mode == 1) {
      ta = mth/(sqrt(2*pi)*x^(3/2))*exp(-(x-mth)^2/(2*x))*x^D[i]*B[i]*exp(-A[i]*x)*x/tao0
    }
    if (mode == 2) {
      ta =  (1/(2*x))*mth/(sqrt(2*pi)*x^(3/2))*exp(-(x-mth)^2/(2*x))*x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
  }

  return(ta)
}