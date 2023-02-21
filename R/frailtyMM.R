#' Fitting frailty models with clustered, multi-event and recurrent data using MM algorithm
#' 
#' @importFrom purrr pmap
#' @param formula Formula where the left hand side is an object of the type \code{Surv}
#' and the right hand side contains the variables and additional specifications. 
#' \code{+cluster()} function specify the group id for clustered data or individual id for recurrent data.
#' \code{+event()} function specify the event id for multi-event data (only two events are allowed).
#' @param data The \code{data.frame} where the formula argument can be evaluated.
#' @param frailty The frailty used for model fitting. The default is "LogN", other choices are
#' "InvGauss", "Gamma" and "PVF". (Note that the computation time for PVF family will be slow 
#' due to the non-explicit expression of likelihood function)
#' @param power The power used if PVF frailty is applied.
#' @param tol The tolerance level for convergence.
#' @param maixt Maximum iterations for MM algorithm.
#' @export

frailtyMM <- function(formula, data, frailty = "LogN", power = NULL, tol = 1e-5, maxit = 200, ...) {
  
  Call <- match.call()
  
  if(!inherits(formula, "formula")) {
    stop("please provide formula object for formula")
  }
  
  if(!inherits(data, "data.frame")) {
    stop("please provide data.frame type for data")
  }
  
  m <- model.frame(formula, data)
  mx <- model.matrix(formula, data)
  
  if (ncol(m[[1]]) == 2) {
    
    cluster_id <- grep("cluster", names(m))
    if (length(cluster_id) == 1) {
      type = "Cluster"
      pb = unlist(gregexpr('\\(', names(m)[cluster_id])) + 1
      pe = unlist(gregexpr('\\)', names(m)[cluster_id])) - 1
      clsname = substr(names(m)[cluster_id], pb, pe)
      remove_cluster_id <- grep(clsname, names(m))
      mx1 = mx[, -c(1, remove_cluster_id)]
      mxid = mx[, cluster_id]
      coef_name = colnames(mx1)
      
      nord = order(mxid)
      mxid = mxid[nord]
      
      N = length(mxid)
      newid = rep(0, N)
      
      if (N <= 2) {
        stop("Please check the sample size of data")
      }
      
      for (i in 2:N) {
        if (mxid[i] > mxid[i-1]) {
          newid[i:N] = newid[i:N] + 1
        }
      }
      
      y = m[[1]][nord, 1]
      X = mx1[nord, ]
      d = m[[1]][nord, 2]
      a = max(newid) + 1
      
      initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
      
      output = frailtyMMcal(y, X, d, N, a, newid,
                            coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                            frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
      
      ret = list(coef = output$coef,
                 est.tht = output$est.tht,
                 lambda = output$lambda,
                 likelihood = output$likelihood,
                 input = output$input,
                 frailty = frailty,
                 power = power,
                 iter = output$iter,
                 convergence = output$convergence,
                 formula = formula,
                 coefname = coef_name,
                 datatype = "Cluster")
    }
    
    event_id <- grep("event", names(m))
    if (length(event_id) == 1) {
      type = "Multiple"
      pb = unlist(gregexpr('\\(', names(m)[event_id])) + 1
      pe = unlist(gregexpr('\\)', names(m)[event_id])) - 1
      evsname = substr(names(m)[event_id], pb, pe)
      remove_event_id <- grep(evsname, names(m))
      mx1 = mx[, -c(1, remove_event_id)]
      mxid = mx[, event_id]
      coef_name = colnames(mx1)
      
      mxid_info = table(mxid)
      n = length(mxid_info)
      b = min(mxid_info)
      p = ncol(mx1)
      if (b != max(mxid_info)) {
        stop("every subject should have same number of events")
      }
      
      nord = order(mxid)
      N = length(nord)
      mx1 = mx1[nord, ]
      X = mx1[nord, ]
      y = m[[1]][nord, 1]
      d = m[[1]][nord, 2]
      
      initGam = frailtyMMcal(y, X, d, N, b, NULL, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 2)
      
      output = frailtyMMcal(y, X, d, N, b, NULL,
                            coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                            frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 2)
      
      ret = list(coef = output$coef,
                 est.tht = output$est.tht,
                 lambda = output$lambda,
                 likelihood = output$likelihood,
                 input = output$input,
                 frailty = frailty,
                 power = power,
                 iter = output$iter,
                 convergence = output$convergence,
                 formula = formula,
                 coefname = coef_name,
                 datatype = "Multi-event")
    }
  }
  
  if (ncol(m[[1]]) == 3) {
    type = "Recurrent"
    cluster_id <- grep("cluster", names(m))
    pb = unlist(gregexpr('\\(', names(m)[cluster_id])) + 1
    pe = unlist(gregexpr('\\)', names(m)[cluster_id])) - 1
    clsname = substr(names(m)[cluster_id], pb, pe)
    remove_cluster_id <- grep(clsname, names(m))
    mx1 = mx[, -c(1, remove_cluster_id)]
    mxid = mx[, cluster_id]
    coef_name = colnames(mx1)
    
    mxid_info = table(mxid)
    p = ncol(mx1)

    nord = order(mxid)
    mxid = mxid[nord]
    
    N = length(mxid)
    newid = rep(0, N)
    
    if (N <= 2) {
      stop("Please check the sample size of data")
    }
    
    for (i in 2:N) {
      if (mxid[i] > mxid[i-1]) {
        newid[i:N] = newid[i:N] + 1
      }
    }
    
    y = m[[1]][nord, 2]
    X = mx1[nord, ]
    d = m[[1]][nord, 3]
    a = max(newid) + 1
    
    initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 3)
    
    output = frailtyMMcal(y, X, d, N, a, newid,
                          coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                          frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 3)
    
    ret = list(coef = output$coef,
               est.tht = output$est.tht,
               lambda = output$lambda,
               likelihood = output$likelihood,
               input = output$input,
               frailty = frailty,
               power = power,
               iter = output$iter,
               convergence = output$convergence,
               formula = formula,
               coefname = coef_name,
               datatype = "Recurrent")
  }
  
  attr(ret, "call") <-  Call
  class(ret) <- "fmm"
  return(ret)
}