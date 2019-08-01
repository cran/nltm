summary.nltm <- function(object,  coef = TRUE, conf.int = 0.95, 
                         digits = max(options()$digits - 4, 3),...)
{
  if(!is.null(object$coef)){
    beta <- object$coef
    nabeta <- !(is.na(beta))          # non-missing coefs
    beta2 <- beta[nabeta]
    if(is.null(beta) | is.null(object$var))
      stop("Input is not valid")
    
    negVar <- which(diag(object$var)<=0)
    if(length(negVar)==0)
      se <- sqrt(diag(object$var))
    else{
      se <- array(NA, dim=length(beta))
      se[-negVar] <- sqrt(diag(object$var)[-negVar])
    }
    rval <- list(call=object$call, convergence=object$convergence,
                 na.action=object$na.action, n=object$n)
    
    if(coef){
      rval$coef <- cbind(beta, exp(beta), se, beta/se,
                         signif(1 - pchisq((beta/se)^2, 1), digits -1))
      dimnames(rval$coef) <- list(names(beta),
                                  c("coef", "exp(coef)", "se(coef)", "z", "p"))
    }
    if(conf.int){
      z <- qnorm((1 + conf.int)/2, 0, 1)
      tmp <- cbind(exp(beta), exp(-beta), exp(beta-z*se), exp(beta+z*se))
      dimnames(tmp) <- list(names(beta),
                            c("exp(coef)", "exp(-coef)",
                              paste("lower .",round(100*conf.int,2),sep = ""),
                              paste("upper .",round(100*conf.int,2),sep = "")))
      rval$conf.int <- tmp
    }
    df <- length(beta2)
    logtest <- -2 * (object$loglik[1] - object$loglik[2])
    rval$logtest <- c(test=logtest, df=df, pvalue=1 - pchisq(logtest, df))
  }else{
    rval <- list(call=object$call, n=object$n, loglik=object$loglik[1])
  }
  rval$npred <- object$npred
  rval$nvar <- object$nvar
  class(rval) <- "summary.nltm"
  rval
}


print.summary.nltm <- function(x, digits = max(options()$digits - 4, 3), ...)
{
  ## message(gettextf("Authors: G. Garibotti, A. Tsodikov"))
  if(!is.null(x$call)) {
    message("Call:")
    message(deparse(x$call))
    message()
  }

  if(!is.null(x$convergence)){
    if(x$convergence==1){
      message("Iteration limit maxit: ", x$maxit, " has been reached.")
      message("Consider a new initial value or optimization parameters such",
          "as bscale.",sep=" ")
      message("See nltm.control help")
    }else{
      if(x$convergence==51){
        message("Warning from optimization method: ", x$message)
        message("Consider a new initial value or optimization parameters such",
            "as bscale.", sep=" ")
        message("See nltm.control help")
      }else{
        if(x$convergence==52){
          message("Error in optimization method: ", x$message)
          message("Consider a new initial value or optimization parameters such",
              "as bscale.", sep=" ")
          message("See nltm.control help")
        }
      }
    }

    negVar <- which(diag(x$var)<=0)
    if(length(negVar)>0){
      message(gettextf("\nWarning message:"))
      message(gettextf("Problem with covariance matrix of coefficients."))
      message(gettextf("Some diagonal terms are negative."))
      message(gettextf("p-values and confidence intervals are unreliable.\n\n"))
    }
  
    savedig <- options(digits = digits)
    on.exit(options(savedig))
  
    message("Non Linear Transformation Model: ", x$call$nlt.model,
        ", fit by maximum likelihood\n", sep="")

    if(x$npred==1){
      if(!is.null(x$coef)) {
        prmatrix(x$coef)
      }
      if(!is.null(x$conf.int)) {
        message()
        prmatrix(x$conf.int)
      }
    }else{
      ind1 <- numeric()
      if(x$nvar$pred.long>0)
        ind1 <- 1:x$nvar$pred.long
      if(!is.na(match("cure",dimnames(x$coef)[[1]])))
        ind1 <- c(ind1,dim(x$coef)[1])    
      if(length(ind1)>0){
        message("Long term predictor")
        if(!is.null(x$coef[ind1,])) 
          prmatrix(matrix(x$coef[ind1,],nrow=length(ind1)),
                   rowlab=rownames(x$coef)[ind1],collab=colnames(x$coef))
        if(!is.null(x$conf.int[ind1,])) {
          message()
          prmatrix(matrix(x$conf.int[ind1,],nrow=length(ind1)),
                   rowlab=rownames(x$conf.int)[ind1],
                   collab=colnames(x$conf.int))
        }
      }
      ind1 <- setdiff(1:dim(x$coef)[1],ind1)    
      if(length(ind1>0)){
        message("\nShort term predictor")
        if(!is.null(x$coef[ind1,]))
          prmatrix(matrix(x$coef[ind1,],nrow=length(ind1)),
                   rowlab=rownames(x$coef)[ind1],collab=colnames(x$coef))
        if(!is.null(x$conf.int[ind1,])) {
          message()
          prmatrix(matrix(x$conf.int[ind1,],nrow=length(ind1)),
                   rowlab=rownames(x$conf.int)[ind1],
                   collab=colnames(x$conf.int))
        }
      }
    }
    message()
    message("Likelihood ratio test=", format(round(x$logtest["test"], 2)),
        " on ", x$logtest["df"], " df, p=", format(x$logtest["pvalue"]),
        sep="")
    message()
    
    omit <- x$na.action
    if(length(omit))
      message("n=", x$n, " (", naprint(omit), ")", sep="")
    else message("n=", x$n, "\n",sep="")
  }else{
    message("Null model: ",x$call$nlt.model,sep="")
    message("Likelihood=", format(round(x$loglik[1], digits)), " n=", x$n,
        sep="")
  }
  invisible(x)
}
