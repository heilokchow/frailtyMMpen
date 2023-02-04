#' @export

print.fmm_summary <- function(x, ...) {
  
  cat("Call:\n")
  dput(attr(x, "call"))
  cat("\n")
  
  cat(paste("  ", x$frailty, "frailty model with", x$datatype, "data, estimated using MM algorithm\n"))
 
  cat("\n")
  
  coef_summary = data.frame(coef = x$coef,
                            "exp(coef)" = exp(x$coef),
                            "se(coef)" = x$coef_sd,
                            zval = x$zcoef,
                            pval = x$pval_coef)
  row.names(coef_summary) = x$coefname
  
  printCoefmat(coef_summary, digits = 3, signif.stars = TRUE)
  
  cat("\n")
  cat(paste("   Frailty parameter, theta:", round(x$est.tht, 3), "(",
            round(x$th_sd, 3), ")", "with p-value = ", sprintf("%.3e", x$pval_th), "\n"))
  cat("\n")
  cat(paste("   Observed likelihood:", round(x$likelihood, 3), "\n"))
  cat("\n")
  cat("   Convergence criteria:")
  cat(paste(" thershold:", sprintf("%.3e", x$convergence), ", number of iterations:", x$iter))

}

#'@export

print.fmm <- function(x, ...) {
  
  cat("Call:\n")
  dput(attr(x, "call"))
  cat("\n")
  
  cat(paste("  ", x$frailty, "frailty model with", x$datatype, "data, estimated using MM algorithm\n"))
  
  cat("\n")
  
  coef_summary = data.frame(coef = x$coef, "exp(coef)" = exp(x$coef))
  row.names(coef_summary) = x$coefname
  
  printCoefmat(coef_summary, digits = 3, signif.stars = TRUE)
  
  cat("\n")
  cat(paste("   Frailty parameter, theta:", round(x$est.tht, 3)))
  cat("\n")
  cat(paste("   Observed likelihood:", round(x$likelihood, 3), "\n"))
  cat("\n")
  cat("   Convergence criteria:")
  cat(paste(" thershold:", sprintf("%.3e", x$convergence), ", number of iterations:", x$iter))
  
}

#'@export

print.fpen <- function(x, ...) {
  
  cat("Call:\n")
  dput(attr(x, "call"))
  cat("\n")
  
  tune = x$tune
  df = colSums(abs(x$coef) > 1e-5)
  
  tune_summary = data.frame(tune = round(tune, 5),
                            Df = df,
                            Lik = round(x$likelihood, 3))
  
  print(tune_summary)
  
}