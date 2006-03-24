nltm <- function(formula=formula(data), data=parent.frame(), subset,
                 na.action, init, control,
                 model=c("PH","PHC","PO","PHPHC","PHPOC","GFM","PHPO"),
                 verbose=FALSE, ...)
{
  cat(gettextf("Authors: G. Garibotti, A. Tsodikov\n"))
  model <- match.arg(model)
  call <- match.call()
  m <- match.call(expand=FALSE)
  temp <- c("", "formula", "data", "subset", "na.action")
  m <- m[ match(temp, names(m), nomatch=0)]
  Terms <- if(missing(data)) terms(formula)
           else              terms(formula, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  if(NROW(m)==0)
    stop("No (non-missing) observations")
  
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  attr(Terms,"intercept")<- 1  # model always has \Lambda_0
  newTerms<-Terms
  
  X <- model.matrix(newTerms,m)
  assign <- lapply(attrassign(X,newTerms)[-1],function(x) x-1)
  X <- X[,-1,drop=FALSE]

# G: needs to be done, summary of optimization options

  # nice idea, check how it works
  controls <- nltm.control(...)
  if(!missing(control)) controls[names(control)] <- control

  if(verbose==TRUE) res <- .C("openDebug", model)
  fit <- nltm.fit(X, Y, model, init, controls, verbose)
  if(verbose==TRUE) res <- .C("closeDebug")

  if(min(diag(fit$var))<=0){
    cat(gettextf("\nWarning message:\n"))
    cat(gettextf("Problem with covariance matrix. Determinant:%s\nExpect problems in the estimation of standard errors.\n\n\n",
                 det(fit$var)))
  }

  class(fit) <- "coxph"
  fit$formula <- formula(Terms)
  fit$call <- call
  fit
}

