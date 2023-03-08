#' @importFrom numDeriv jacobian
#' @export
#' 
predict.fmm <- function(object, newdata = NULL, surv = FALSE, ...) {
  
  coef = object$coef
  time_eval = sort(object$input$y)
  
  frailtyc = switch(object$frailty, "Gamma" = 0, "LogN" = 1, "InvGauss" = 2, "PVF" = 3)
  datatype = switch(object$datatype, "Cluster" = 1, "Multi-event" = 2, "Recurrent" = 3)
  p = length(object$coef)
  
  if (is.null(object$power)) {
    power = 0.0
  } else {
    power = object$power
  }
  
  Fisher = -numDeriv::hessian(logLik,
                              x = c(object$est.tht, object$coef),
                              method="Richardson",
                              data = object$input,
                              lambda = object$lambda,
                              frailtyc = frailtyc,
                              id = object$id,
                              N = object$N,
                              a = object$a,
                              p = p,
                              power = power,
                              type = datatype)
  
  coef_cv = solve(Fisher)
  
  
  new_covariate = c()
  
  if (!is.null(newdata)) {
    for (i in object$coefname) {
      pos = which(colnames(newdata) == i)
      
      if (length(pos) == 0) {
        stop("The new data dose not have the same format as original dataset.")
      }
      
      if (length(newdata[pos][[1]]) > 1) {
        stop("Please provide one vector of new data.")
      }
      new_covariate = c(new_covariate, newdata[pos][[1]])
    }
  } else {
    new_covariate = rep(0, length(coef))
  }
  
  
  if (object$datatype == "Multi-event") {
    
    m <- model.frame(object$formula, newdata)
    cluster_id <- grep("event", names(m))
    pb = unlist(gregexpr('\\(', names(m)[cluster_id])) + 1
    pe = unlist(gregexpr('\\)', names(m)[cluster_id])) - 1
    clsname = substr(names(m)[cluster_id], pb, pe)
    
    pos = which(colnames(newdata) == clsname)
    
    if (length(pos) == 0) {
      stop("The new data dose not have the information about the event's id.\n
               Please ckeck whether the newdata is of the same format as original data.")
    }
    new_event = newdata[pos][[1]]
    
    estimate = hazardME(c(1, coef), object = object, new_covariate = new_covariate, new_event = new_event, Surv = surv)
    
    gradT = jacobian(hazardME, c(1, coef), object = object, new_covariate = new_covariate, new_event = new_event, Surv = surv)
    
    estimateCOV = gradT %*% coef_cv %*% t(gradT)
    estimateSD = sqrt(diag(estimateCOV))
    
  }
  
  if (object$datatype == "Cluster" || object$datatype == "Recurrent") {
    
    if (object$datatype == "Cluster") {
      
      estimate = hazardCL(c(1, coef), object = object, new_covariate = new_covariate, Surv = surv)
      
      gradT = jacobian(hazardCL, c(1, coef), object = object, new_covariate = new_covariate, Surv = surv)
      
      estimateCOV = gradT %*% coef_cv %*% t(gradT)
      estimateSD = sqrt(diag(estimateCOV))
      
    }
    
    if (object$datatype == "Recurrent") {
      
      estimate = hazardRE(c(1, coef), object = object, new_covariate = new_covariate, Surv = surv)
      
      gradT = jacobian(hazardRE, c(1, coef), object = object, new_covariate = new_covariate, Surv = surv)
      
      estimateCOV = gradT %*% coef_cv %*% t(gradT)
      estimateSD = sqrt(diag(estimateCOV))
      
    }
  }
  
  if (surv == TRUE) {
    output = data.frame(estimate = estimate,
                        estimateSD = estimateSD,
                        estimatelb = pmax(estimate - 1.96 * estimateSD, 0),
                        estimateub = pmin(estimate + 1.96 * estimateSD, 1))
  } else {
    output = data.frame(estimate = estimate,
                        estimateSD = estimateSD,
                        estimatelb = pmax(estimate - 1.96 * estimateSD, 0),
                        estimateub = estimate + 1.96 * estimateSD)
  }
  
  return(output)
}

hazardCL <- function(coefall, object, new_covariate, Surv = FALSE) {
  
  vy = object$input$y
  vd = object$input$d
  N = length(vy)
  wi = coefall[1]
  coef = coefall[-1]
  int2 = object$Ar
  YpreExp = exp(object$input$X %*% coef)
  
  ME = object$id + 1
  for (i in 1:N) {
    ME[i] = int2[ME[i]]
  }
  
  E_0 = as.vector(ME*YpreExp)
  SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
  SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
  
  lambda = vd/SUM_0
  
  La = (cumsum(lambda[order(vy)]))
  
  Ha = La*wi*exp(new_covariate %*% coef)[[1]]
  
  if (Surv == FALSE) {
    return(Ha)
  }
  
  return(exp(-Ha))
  
}


hazardME <- function(coefall, object, new_covariate, new_event, Surv = FALSE) {
  
  vy = object$input$y
  vd = object$input$d
  
  a = object$a
  y = vy[((new_event-1)*a + 1):(new_event*a)]
  d = vd[((new_event-1)*a + 1):(new_event*a)]
  X = object$input$X[((new_event-1)*a + 1):(new_event*a),]
  
  wi = coefall[1]
  coef = coefall[-1]
  int2 = object$Ar
  YpreExp = exp(X %*% coef)
  
  E_0 = int2*YpreExp
  SUM_0 = cumsum((E_0[order(y)])[seq(a,1,-1)])
  SUM_0 = (SUM_0[seq(a,1,-1)])[rank(y)] 
  
  lambda = d/SUM_0
  
  La = (cumsum(lambda[order(y)]))
  
  Ha = La*wi*exp(new_covariate %*% coef)[[1]]
  
  if (Surv == FALSE) {
    return(Ha)
  }
  
  return(exp(-Ha))
  
}

hazardRE <- function(coefall, object, new_covariate, Surv = FALSE) {
  
  y = object$input$y
  d = object$input$d
  N = length(y)
  wi = coefall[1]
  coef = coefall[-1]
  id = object$id + 1
  a = object$a
  X = object$input$X
  int2 = object$Ar
  XM = array(0, c(a, p, N))
  
  for (i in N:1) {
    tempID = id[i]
    ct = y[i]
    if (i == N || tempID < id[i+1]) {
      for (j in 1:N) {
        XM[tempID,,j] = X[i, ]
      }
    } else {
      for (j in 1:N) {
        if (ct >= y[j]) {
          XM[tempID,,j] = X[i, ]
        }
      }
    }
  }
  
  XMexp = apply(XM, 3, function(x, coef) {exp(x %*% coef)}, coef = coef)
  
  Yend = rep(0, a)
  for (i in 1:N) {
    tempID = id[i]
    if (i == N || tempID < id[i+1]) {
      Yend[tempID] = y[i]
    }
  }
  
  SUM0 = rep(0, N)
  for (i in 1:N) {
    for (j in 1:a) {
      SUM0[i] = SUM0[i] + (Yend[j] >= y[i]) * int2[j] * XMexp[j, i]
    }
  }
  
  lambda = d / SUM0
  
  La = (cumsum(lambda[order(y)]))
  
  Ha = La*wi*exp(new_covariate %*% coef)[[1]]
  
  if (Surv == FALSE) {
    return(Ha)
  }
  
  return(exp(-Ha))
  
}