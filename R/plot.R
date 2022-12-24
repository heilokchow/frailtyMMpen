plot.fpen <- function(x, ...) {
  xaxis = log(x$tune)
  yaxis = t(x$coef)

  n = ncol(yaxis)
  matplot(xaxis, yaxis, lty=1, xlab="log(tune)",ylab="coefficients",type="l")
  abline(v = log(x$tune.min), col="red", lwd=3, lty=2)
}