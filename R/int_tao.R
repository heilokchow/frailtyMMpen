int_tao <- function(x, i, est.tht, A, B, D, tao0 = 0, frailty = "LogN", power = NULL, mode = 0) {
  
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
    if (mode == 3) {
      ta = (-0.5*log(2*pi*si2)-(log(x))^2/(2*si2))*1/(x*sqrt(2*pi*si2))*exp(-(log(x))^2/(2*si2))*x^D[i]*B[i]*exp(-A[i]*x)
    }
  }

  if (frailty == "InvGauss") {
    alp = est.tht
    if (mode == 0) {
      ta = 1/(sqrt(2*pi*alp)*x^(3/2))*exp(-(x-1)^2/(2*x*alp))*x^D[i]*B[i]*exp(-A[i]*x)
    }
    if (mode == 1) {
      ta = 1/(sqrt(2*pi*alp)*x^(3/2))*exp(-(x-1)^2/(2*x*alp))*x^D[i]*B[i]*exp(-A[i]*x)*x/tao0
    }
    if (mode == 2) {
      ta = ((x-1)^2/x)/(sqrt(2*pi*alp)*x^(3/2))*exp(-(x-1)^2/(2*x*alp))*x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
  }
  
  if (frailty == "PVF") {
    if (mode == 0) {
      ta = dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)
    }
    if (mode == 1) {
      ta = dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)*x/tao0
    }
    if (mode == 2) {
      ta = ldTweedie(x, mu = 1, phi = est.tht, p = power)[2]*exp(ldTweedie(x, mu = 1, phi = est.tht, p = power)[1])*
        x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
    if (mode == 3) {
      ta = ldTweedie(x, mu = 1, phi = est.tht, p = power)[3]*exp(ldTweedie(x, mu = 1, phi = est.tht, p = power)[1])*
        x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
    if (mode == 4) {
      ta = ldTweedie(x, mu = 1, phi = est.tht, p = power)[1]*exp(ldTweedie(x, mu = 1, phi = est.tht, p = power)[1])*
        x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
  }

  return(ta)
}

int_tao1 <- function(x, i, est.tht, A, B, D, tao0 = 0, frailty = "LogN", power = NULL, mode = 0, fvar = 1) {
  
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
    if (mode == 3) {
      ta = (-0.5*log(2*pi*si2)-(log(x))^2/(2*si2))*1/(x*sqrt(2*pi*fvar))*exp(-(log(x))^2/(2*fvar))*x^D[i]*B[i]*exp(-A[i]*x)
    }
  }
  
  if (frailty == "InvGauss") {
    alp = est.tht
    if (mode == 0) {
      ta = 1/(sqrt(2*pi*alp)*x^(3/2))*exp(-(x-1)^2/(2*x*alp))*x^D[i]*B[i]*exp(-A[i]*x)
    }
    if (mode == 1) {
      ta = 1/(sqrt(2*pi*alp)*x^(3/2))*exp(-(x-1)^2/(2*x*alp))*x^D[i]*B[i]*exp(-A[i]*x)*x/tao0
    }
    if (mode == 2) {
      ta = ((x-1)^2/x)/(sqrt(2*pi*alp)*x^(3/2))*exp(-(x-1)^2/(2*x*alp))*x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
  }
  
  if (frailty == "PVF") {
    if (mode == 0) {
      ta = dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)
    }
    if (mode == 1) {
      ta = dtweedie(x, mu = 1, phi = est.tht, power = power)*x^D[i]*B[i]*exp(-A[i]*x)*x/tao0
    }
    if (mode == 2) {
      ta = ldTweedie(x, mu = 1, phi = est.tht, p = power)[2]*exp(ldTweedie(x, mu = 1, phi = fvar, p = power)[1])*
        x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
    if (mode == 3) {
      ta = ldTweedie(x, mu = 1, phi = est.tht, p = power)[3]*exp(ldTweedie(x, mu = 1, phi = fvar, p = power)[1])*
        x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
    if (mode == 4) {
      ta = ldTweedie(x, mu = 1, phi = est.tht, p = power)[1]*exp(ldTweedie(x, mu = 1, phi = fvar, p = power)[1])*
        x^D[i]*B[i]*exp(-A[i]*x)/tao0
    }
  }
  
  return(ta)
}