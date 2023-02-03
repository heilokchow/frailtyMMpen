#' Plot the regularization path
#' 
#' @description Plot the whole regularization path run by frailtyMMpen
#' @usage
#' ##S3 method for class "fpen"
#' @param x Object with class "fpen"
#' @param ... Further arguments pass to or from other methods
#' @export
#' 
plot.fpen <- function(x, ...) {
  xaxis = log(x$tune)
  yaxis = t(x$coef)

  n = ncol(yaxis)
  matplot(xaxis, yaxis, lty=1, xlab="log(tune)",ylab="coefficients",type="l")
  abline(v = log(x$tune.min), col="red", lwd=3, lty=2)
}