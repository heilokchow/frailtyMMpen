library(survival)

a <- factor(rep(1:3,4), labels=c("low", "medium", "high"))
b <- factor(rep(1:4,3))
levels(strata(b))
levels(strata(a,b,shortlabel=TRUE))

f <- function(x) {
  x^2
}
coxph(Surv(futime, fustat) ~ age + strata(rx) + f(age), data=ovarian) 

test = Surv(futime, fustat) ~ age + strata(rx) + f(age)
test = model.frame(test, data = ovarian)

data = ovarian
special <- c("strata", "cluster", "slope")

Terms <- if (missing(data)) terms(test, special) else terms(test, special, data = ovarian)   
ord <- attr(Terms, "order")
#  if (length(ord) & any(ord != 1)) 
#     stop("Interaction terms are not valid for this function")

call
m <- match.call(expand.dots = FALSE)
m$correlation <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$... <- NULL
m$formula <- Terms
m[[1]] <- as.name("model.frame")
m <- eval(test)

model.extract(test, "response")