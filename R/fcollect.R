MMroutine <- function(y, X, d, coef, lambda, tht, frailty, penalty, gam, tune, id, N, a, p, power, type) {
  
  TEST = NULL
  
  if (type == 1) {
    TEST = MMCL(y, X, d, coef, lambda, tht, frailty, penalty, gam, tune, id, N, a, p, power)
  }
  
  if (type == 2) {
    TEST = MMME(y, X, d, coef, lambda, tht, frailty, penalty, gam, tune, N, a, p, power)
  }
  
  if (type == 3) {
    # TEST = MMRE(y, X, d, coef, lambda, tht, frailty, penalty, gam, tune, id, N, a, p, power)
    TEST = MMRELS(y, X, d, coef, lambda, tht, frailty, penalty, gam, tune, id, N, a, p, power)
  }
  
  return(TEST)
}

#'@export
LogHessianCL <- function(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power) {
  
  neworder = order(id, decreasing = FALSE)
  
  y = y[neworder]
  X = X[neworder, ]
  d = d[neworder]
  newid = id[neworder]
  lambda = lambda[neworder]
  
  return(LogLikHessianCL(y, X, d, coef, lambda, est.tht, frailtyc, newid - 1, N, a, p, power))
  
}

logLikcal <- function(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type) {
  
  s = 0
  
  if (type == 1) {
    s = LogLikCL(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power)
  }
  
  if (type == 2) {
    s = LogLikME(y, X, d, coef, lambda, est.tht, frailtyc, N, a, p, power)
  }
  
  if (type == 3) {
    s = LogLikRE(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power)
  }
  
  return(s)
  
}

#'@export
logLik <- function(x1, data, lambda, frailtyc, id, N, a, p, power, type) {
  
  est.tht = x1[1]
  coef = x1[-1]
  y = data$y
  X = data$X
  d = data$d
  
  return(logLikcal(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type))
  
}

#' #' cluster function
#' #'
#' #' @param x name from original dataframe which specifies the cluster or object id.
#' #' @description {Specify cluster id for clustered data or object id for recurrent data in the input \code{formula}.}
#' #' @return No return value, called to construct \code{formula}.
#' #' @export
#' cluster <- function(x) {
#'   x
#' }

#' event function
#' 
#' @param x name from original dataframe which specifies the event id.
#' @description {Specify event id for multi-event data in the input \code{formula}.}
#' @return No return value, called to construct \code{formula}.
#' @export
event <- function(x) {
  x
}