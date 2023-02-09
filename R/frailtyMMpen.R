#' Fitting penalized frailty models with clustered, multi-event and recurrent data using MM algorithm
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
#' @param penalty The penalty used for regulization, the default is "LASSO", other choices are "MCP" and "SCAD".
#' @param tune The sequence of tuning parameters provided by user. If not provided, the default grid will be applied.
#' @param tol The tolerance level for convergence.
#' @param maixt Maximum iterations for MM algorithm.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib frailtyMMpen, .registration = TRUE
#' 
frailtyMMpen <- function(formula, data, frailty = "LogN", power = NULL, penalty = "LASSO", tune = NULL, tol = 1e-5, maxit = 200) {
  
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
  }
  
  threshold = tol
  
  if (is.null(tune)) {
    tuneseq = exp(seq(-5.5, 1, 0.25))
  } else {
    tuneseq = tune
  }
  
  if (type == "Cluster") {
    
    p = dim(X)[3]
    a = nrow(y)
    b = ncol(y)
    N = a*b
    
    ini = frailtyMM_CL(y, X, d, frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = threshold)
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMM_CL(y, X, d, 
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0,
                         frailty = frailty, power = power, penalty = penalty, tune = tuneseq[z], maxit = maxit, threshold = threshold)
      
      coef0 = cur$coef
      est.tht0 = cur$est.tht
      lambda0 = cur$lambda
      likelihood0 = cur$likelihood
      
      coef_all[[z]] = coef0
      est.tht_all[[z]] = est.tht0
      lambda_all[[z]] = lambda0
      likelihood_all[[z]] = likelihood0
      BIC_all[[z]] = -2*likelihood0 + max(1, log(log(p + 1)))*(sum(abs(coef0) > threshold) + 1)*log(N)
      
      if (sum(abs(coef0)) < threshold) {
        cat(sum(abs(coef0)), "????\n")
        break
      }
      
      cat(z, "---------\n")
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(coef_all)
    likelihood_all = unlist(likelihood_all)
    BIC_all = unlist(BIC_all)
    
    output = list(coef = coef_all,
                  est.tht = est.tht_all,
                  lambda = lambda_all,
                  likelihood = likelihood_all,
                  BIC = BIC_all,
                  tune = tuneseq[seq_len(z)],
                  tune.min = tuneseq[which.min(BIC_all)],
                  y = y,
                  X = X,
                  d = d)
  }
  
  if (type == "Recurrent") {
    # NO PVF for RE
    
    p = dim(X[[1]])[2]
    n = length(y)
    
    ini = frailtyMM_RE(y, X, d, frailty = frailty, penalty = NULL, maxit = maxit, threshold = threshold)
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMM_RE(y, X, d, 
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0,
                         frailty = frailty, penalty = penalty, tune = tuneseq[z], maxit = maxit, threshold = threshold)
      
      coef0 = cur$coef
      est.tht0 = cur$est.tht
      lambda0 = cur$lambda
      likelihood0 = cur$likelihood
      
      coef_all[[z]] = coef0
      est.tht_all[[z]] = est.tht0
      lambda_all[[z]] = lambda0
      likelihood_all[[z]] = likelihood0
      BIC_all[[z]] = -2*likelihood0 + max(1, log(log(p + 1)))*(sum(abs(coef0) > 1e-6) + 1)*log(n)
      
      if (sum(abs(coef0)) < 1e-6) {
        cat(sum(abs(coef0)), "????\n")
        break
      }
      
      cat(z, "---------\n")
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(coef_all)
    likelihood_all = unlist(likelihood_all)
    BIC_all = unlist(BIC_all)
    
    output = list(coef = coef_all,
                  est.tht = est.tht_all,
                  lambda = lambda_all,
                  likelihood = likelihood_all,
                  BIC = BIC_all,
                  tune = tuneseq[seq_len(z)],
                  tune.min = tuneseq[which.min(BIC_all)],
                  y = y,
                  X = X,
                  d = d)
  } 
  
  if (type == "Multiple") {
    # NO PVF for ME
    
    p = dim(X)[3]
    n = ncol(y)
    
    ini = frailtyMM_ME(y, X, d, frailty = frailty, penalty = NULL, maxit = maxit, threshold = threshold)
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda01 = ini$lambda1
    lambda02 = ini$lambda2
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda1_all = list()
    lambda2_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMM_ME(y, X, d, 
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda1.ini = lambda01, lambda2.ini = lambda02,
                         frailty = frailty, penalty = penalty, tune = tuneseq[z], maxit = maxit, threshold = threshold)
      
      coef0 = cur$coef
      est.tht0 = cur$est.tht
      lambda01 = cur$lambda1
      lambda02 = cur$lambda2
      likelihood0 = cur$likelihood
      
      coef_all[[z]] = coef0
      est.tht_all[[z]] = est.tht0
      lambda1_all[[z]] = lambda01
      lambda2_all[[z]] = lambda02
      likelihood_all[[z]] = likelihood0
      BIC_all[[z]] = -2*likelihood0 + max(1, log(log(p + 1)))*(sum(abs(coef0) > 1e-6) + 1)*log(n)
      
      if (sum(abs(coef0)) < 1e-6) {
        cat(sum(abs(coef0)), "????\n")
        break
      }
      
      cat(z, "---------\n")
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(coef_all)
    likelihood_all = unlist(likelihood_all)
    BIC_all = unlist(BIC_all)
    
    output = list(coef = coef_all,
                  est.tht = est.tht_all,
                  lambda1 = lambda1_all,
                  lambda2 = lambda2_all,
                  likelihood = likelihood_all,
                  BIC = BIC_all,
                  tune = tuneseq[seq_len(z)],
                  tune.min = tuneseq[which.min(BIC_all)],
                  y = y,
                  X = X,
                  d = d)
  } 
  
  attr(output, "call") <-  Call
  class(output) = "fpen"
  output
}