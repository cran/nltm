# Predictor 1: long term effect
# Predictor 2: short term effect

nltm <- function(formula1=formula(data), formula2=formula(data),
                 data=parent.frame(), subset, na.action, init, control,
                 model=c("PH","PHC","PO","PHPHC","PHPOC","GFM","PHPO"),
                 verbose=FALSE, ...)
{
  cat(gettextf("Authors: G. Garibotti, A. Tsodikov\n"))
  model <- match.arg(model)
  call <- match.call()
  m <- match.call(expand=FALSE)
  # this is necessary because otherwise eval(m, parent.frame()) doesn't work
  names(m)[names(m)=="formula1"] <- "formula"  
  temp <- c("","formula","data","subset","na.action")
  m <- m[match(temp, names(m), nomatch=0)]
  Terms1 <- if(missing(data)) terms(formula1)
            else              terms(formula1, data=data)
  m$formula <- Terms1
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  if(NROW(m)==0)
    stop("No (non-missing) observations")
  
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  attr(Terms1,"intercept") <- 1  
  newTerms <- Terms1 
  
  X1 <- model.matrix(newTerms,m)
  assign <- lapply(attrassign(X1,newTerms)[-1],function(x) x-1)
  X1 <- X1[,-1,drop=FALSE]

  npred <- nPredictor(model)
  if(!missing(formula2)){
    if(npred==1){
      cat(gettextf("\nWarning message:\n"))
      cat(gettextf(paste("Model", model, sep=" ")))
      cat(gettextf(" has only one predictor however there are two formulas,\nformula2 will not be used.\n\n"))
    }else{
      m <- match.call(expand=FALSE)
      names(m)[names(m)=="formula1"] <- "formula"  
      temp <- c("","formula","data","subset","na.action")
      m <- m[match(temp, names(m), nomatch=0)]
      Terms2 <- if(missing(data)) terms(formula2)
      else              terms(formula2, data=data)
      m$formula <- Terms2
      m[[1]] <- as.name("model.frame")
      m <- eval(m, parent.frame())
      if(NROW(m)==0)
        stop("No (non-missing) observations")

      attr(Terms2,"intercept")<- 1  
  
      X2 <- model.matrix(Terms2,m)
      assign <- lapply(attrassign(X2,Terms2)[-1],function(x) x-1)
      X2 <- X2[,-1,drop=FALSE]
    }
  }else{
    if(npred>1){
      X2 <- NULL
      Terms2 <- Terms1
    }
  }
  
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(nltm.control))  #legal arg names
    indx <- match(names(extraArgs), controlargs, nomatch=0)
    if (any(indx==0))
      stop("Argument ", names(extraArgs)[indx==0], "not matched")
  }
  controls <- nltm.control(...)
  if(!missing(control)) controls[names(control)] <- control
  
  if(verbose==TRUE) res <- .C("openDebug", model)
  fit <- nltm.fit(X1, X2, Y, model, init, controls, verbose)
  if(verbose==TRUE) res <- .C("closeDebug")

  class(fit) <- "nltm"
  fit$formula <- formula(Terms1)
  if(npred>1)
    fit$formula <- list(predictor1=formula(Terms1), predictor2=formula(Terms2))
  fit$call <- call
  na.action <- attr(m, "na.action")
  if (length(na.action))
    fit$na.action <- na.action
  fit
}

