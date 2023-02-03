#' Provide the summary for the model fitting 
#' 
#' @importFrom numDeriv hessian
#' @usage
#' ##S3 method for class "fmm"
#' @param model Object with class "fmm", generated from \code{frailtyMM}.
#' @param ... Ignored
#' 
#' @details the summary for the model of frailtyMM.
#' The standard error and p-value of estimated parameters are based on Fisher Information martix.

#' @export
summary.fmm <- function(model, ...) {
  
  Fisher = -numDeriv::hessian(loglik,
                              x = c(model$est.tht, model$coef),
                              method="Richardson",
                              data = model$input,
                              lambda = model$lambda,
                              lambda2 = model$lambda2,
                              frailty = model$frailty,
                              power = model$power,
                              datatype = model$datatype)
  
  sd = diag(solve(Fisher))
  
  th_sd = sqrt(sd[1])
  coef_sd = sqrt(sd[-1])
  
  zth = model$est.tht / th_sd
  zcoef = model$coef / coef_sd
  
  pval_th = 1 - pnorm(zth)
  pval_coef = 1 - pnorm(zcoef)
  
  ret = list(est.tht = model$est.tht,
             coef = model$coef,
             th_sd = th_sd,
             coef_sd = coef_sd,
             zth = zth,
             zcoef = zcoef,
             pval_th = pval_th,
             pval_coef = pval_coef,
             frailty = model$frailty,
             power = model$power,
             datatype = model$datatype,
             iter = model$iter,
             convergence = model$convergence)
  
  class(ret) = "fmm_summary"
  return(ret)
}