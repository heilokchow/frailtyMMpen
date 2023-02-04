#'@export

coef.fmm <- function(x, ...) {
  
  ret = data.frame(coef = x$coef)
  row.names(ret) = x$coefname
  
  return(ret)
} 

#'@export

coef.fpen <- function(x, tune = NULL) {
  
  n = length(x$tune)
  
  if (is.null(tune)) {
    
    retx = x$coef[, which(x$tune == x$tune.min)]
  
    } else if (length(tune) > 1) {
    
    stop("Please provide one tuning parameter only.")
  
    } else {
    
    if (tune < min(x$tune) || tune > max(x$tune)) {
      stop("The provided tuning parameter is out of the range of tuning parameters used for construct regularization path.")
    }
    
    if (tune == min(x$tune)) {
      retx = x$coef[, 1]
    } else {
      it = which.min(tune >= x$tune)
      s1 = tune - x$tune[it-1]
      s2 = x$tune[it] - tune
      s = s1 + s2
      
      retx = s1/s *  x$coef[it] + s2/s *  x$coef[it-1]
    }
  }
  
  ret = data.frame(coef = unname(retx))
  return(ret)
  
}