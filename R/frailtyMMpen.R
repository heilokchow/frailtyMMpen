frailtyMMpen <- function(y, X, d, frailty = "LogN", type = "Cluster", power = NULL, penalty = "LASSO", tune = NULL, maxit = 200, threshold = 1e-6) {

  tuneseq = exp(seq(-5.5, 1, 0.25))
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
                  tune.min = tuneseq[which.min(BIC_all)])
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
                  tune.min = tuneseq[which.min(BIC_all)])
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
                  tune.min = tuneseq[which.min(BIC_all)])
  } 
  
  class(output) = "fpen"
  output
}