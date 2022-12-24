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
    
    output = frailtyMM_ME(y, X, d, frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol)
  }
  
  return(output)
}