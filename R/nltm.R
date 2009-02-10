# Predictor 1: long term effect
# Predictor 2: short term effect

nltm <- function(formula1=formula(data), formula2=formula(data),
                 data=parent.frame(), subset, na.action, init=NULL,
                 control, nlt.model=c("PH","PHC","PO","PHPHC","PHPOC",
                            "GFM","PHPO"),
                 model=FALSE, x=FALSE, y=FALSE, verbose=FALSE, ...)
{
  if(sys.parent()==0)
    cat(gettextf("Authors: G. Garibotti, A. Tsodikov\n"))
  if(!nlt.model %in% eval(formals()[["nlt.model"]]))
    stop(gettextf("nlt.model should be one of %s",
                  paste(dQuote(eval(formals()[["nlt.model"]])),collapse=", ")),
         domain=NA)

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

  attr(Terms1,"intercept") <- 1  # model always has \Lambda_0
  newTerms <- Terms1 
  
  X1 <- model.matrix(newTerms,m)
  assign <- lapply(attrassign(X1,newTerms)[-1],function(x) x-1)
  X1 <- X1[,-1,drop=FALSE]

  npred <- nPredictor(nlt.model)
  if(model & npred>1)
    m1 <- m
  
  if(!missing(formula2)){
    if(npred==1){
      cat(gettextf("\nWarning message:\n"))
      cat(gettextf(paste("Model", nlt.model, sep=" ")))
      cat(gettextf(" has only one predictor however there are two formulas,\nformula2 will not be used.\n\n"))
    }else{
      m <- match.call(expand=FALSE)
      names(m)[names(m)=="formula2"] <- "formula"  
      temp <- c("","formula","data","subset","na.action")
      m <- m[match(temp, names(m), nomatch=0)]
      Terms2 <- if(missing(data)) terms(formula2)
                else              terms(formula2, data=data)
      m$formula <- Terms2
      m[[1]] <- as.name("model.frame")
      m <- eval(m, parent.frame())
      if(NROW(m)==0)
        stop("No (non-missing) observations")

      attr(Terms2,"intercept")<- 1  # model always has \Lambda_0
  
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
  
  if(verbose!=FALSE){
    res <- .C("openDebug", verbose)
    verbose <- TRUE
  }
  fit <- nltm.fit(X1, X2, Y, nlt.model, init, controls, verbose)
  if(verbose==TRUE) res <- .C("closeDebug")

  if(npred==1){
    fit$formula <- formula(Terms1)
    fit$terms <- Terms1
    fit$nvar <- dim(X1)[2]
  }else{
    fit$formula$pred.long <- formula(Terms1)
    fit$formula$pred.short <- formula(Terms2)
    fit$terms$pred.long <- Terms1
    fit$terms$pred.short <- Terms2
    fit$nvar$pred.long <- dim(X1)[2]
    if(!is.null(X2))
      fit$nvar$pred.short <- dim(X2)[2]
    else
      fit$nvar$pred.short <- dim(X1)[2]
  }
  fit$call <- call
  if(is.null(fit$call$nlt.model))
    fit$call$nlt.model <- "PH"
     
  na.action <- attr(m, "na.action")
  if (length(na.action))
    fit$na.action <- na.action
  fit$npred <- npred
  if(x)
    if(npred==1){
      fit$x <- X1
    }else{
      fit$x$pred.long <- X1
      if(!is.null(X2))
        fit$x$pred.short <- X2
      else
        fit$x$pred.short <- X1
    }
  if(model)
    if(npred==1){
      fit$model <- m
    }else{
      fit$model$pred.long <- m1
      fit$model$pred.short <- m
    }
  if(y)
    fit$y <- Y
  
  class(fit) <- "nltm"
  
  fit
}

