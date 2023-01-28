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

frailtyMM <- function(formula, data, frailty = "LogN", power = NULL, tol = 1e-6, maxit = 200, ...) {
  
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
      
      mxid_info = table(mxid)
      a = length(mxid_info)
      b = min(mxid_info)
      p = ncol(mx1)
      if (b != max(mxid_info)) {
        stop("each cluster should have same amount of objects")
      }
      
      nord = order(mxid)
      mx1 = mx1[nord, ]
      X = array(c(mx1), c(b, a, p))
      X = aperm(X, c(2, 1, 3))
      y = matrix(m[[1]][nord, 1], c(a, b), byrow = TRUE)
      d = matrix(m[[1]][nord, 2], c(a, b), byrow = TRUE)
      
      output = frailtyMM_CL(y, X, d, frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol)
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
      
      mxid_info = table(mxid)
      n = length(mxid_info)
      b = min(mxid_info)
      p = ncol(mx1)
      if (b != max(mxid_info) || b != 2) {
        stop("every subject should have exactly two events")
      }
      
      nord = order(mxid)
      mx1 = mx1[nord, ]
      X = array(c(mx1), c(n, 2, p))
      X = aperm(X, c(2, 1, 3))
      y = matrix(m[[1]][nord, 1], c(2, n), byrow = TRUE)
      d = matrix(m[[1]][nord, 2], c(2, n), byrow = TRUE)
      
      output = frailtyMM_ME(y, X, d, frailty = frailty, penalty = NULL, maxit = maxit, threshold = tol)
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
    
    mxid_info = table(mxid)
    n = length(mxid_info)
    b = min(mxid_info)
    p = ncol(mx1)
    
    nord = order(mxid)
    mx1 = mx1[nord, ]
    y1 = m[[1]][nord, 2]
    d1 = m[[1]][nord, 3]
    
    dfX = as.data.frame(cbind(mxid, mx1))
    X = split(dfX, f = dfX$mxid)
    X = lapply(X, function(x) {x = as.matrix(x[, -1])})
    
    y = split(y1, f = mxid)
    d = split(d1, f = mxid)
    
    output = frailtyMM_RE(y, X, d, frailty = frailty, penalty = NULL, maxit = maxit, threshold = tol)
  }
  
  return(output)
}