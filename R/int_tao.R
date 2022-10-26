int_tao <- function(w, i, est.tht, A, D, frailty = "LogN", mode = 0) {
  
  ta = 0
  if (frailty == "LogN") {
    si2 = est.tht
    if (mode == 0) {
      ta = w^(D[i]-1)*exp(5 - (log(w))^2/(2*si2) - A[i]*w)
    }
    if (mode == 1) {
      ta = w^(D[i]-1)*(log(w))^2*exp(5 - (log(w))^2/(2*si2) - A[i]*w)
    }
    if (mode == 2) {
      ta = w^(D[i])*exp(5 - (log(w))^2/(2*si2) - A[i]*w)
    }
  }

  if (frailty == "InvGauss") {
    mu = est.tht
    if (mode == 0) {
      ta = w^(D[i]-3/2)*exp(-mu^2/(2*w) - w*(1/2 + A[i]))
    }
    if (mode == 1) {
      ta = w^(D[i]-5/2)*exp(-mu^2/(2*w) - w*(1/2 + A[i]))
    }
    if (mode == 2) {
      ta = w^(D[i]-1/2)*exp(-mu^2/(2*w) - w*( 1/2 + A[i]))
    }
  }

  return(ta)
}