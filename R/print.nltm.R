print.nltm <- function(x, digits=max(options()$digits - 4, 3), ...)
{
  ## message(gettextf("Authors: G. Garibotti, A. Tsodikov"))

  if(!is.null(cl<- x$call)) {
    message("Call:")
    message(deparse(cl))
    message()
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if(!is.null(x$coef)){
    coef <- x$coef
    negVar <- which(diag(x$var)<=0)
    if(length(negVar)==0)
      se <- sqrt(diag(x$var))
    else{
      se <- array(NA, dim=length(coef))
      se[-negVar] <- sqrt(diag(x$var)[-negVar])
    }
    
    if(is.null(coef) | is.null(se))
      stop("Input is not valid")

    if(x$convergence==1){
      message("Iteration limit maxit: ", x$maxit, " has been reached.")
      message("Consider a new initial value or optimization parameters such as",
          "bscale.", sep=" ")
      message("See nltm.control help.")
    }else{
      if(x$convergence==51){
        message("Warning from optimization method: ", x$message)
        message("Consider a new initial value or optimization parameters such as",
            "bscale.", sep=" ")
        message("See nltm.control help.")
      }else{
        if(x$convergence==52){
          message("Error in optimization method: ", x$message)
          message("Consider a new initial value or optimization parameters such",
              "as bscale.", sep=" ")
          message("See nltm.control help.")
        }
      }
    }

    if(length(negVar)>0){
      message(gettextf("\nWarning message:"))
      message(gettextf("Problem with covariance matrix of coefficients."))
      message(gettextf("Some diagonal terms are negative."))
      message(gettextf("p-values and confidence intervals are unreliable.\n\n"))
    }
    
    message("Non Linear Transformation Model: ", x$call$nlt.model,
        ", fit by maximum likelihood\n",sep="")
    tmp <- cbind(coef, exp(coef), se, coef/se,
                 signif(1 - pchisq((coef/ se)^2, 1), digits-1))
    dimnames(tmp) <- list(names(coef),
                          c("coef","exp(coef)","se(coef)","z","p"))
    
    if(x$npred==1){
      message()
      prmatrix(tmp)
    }else{
      ind1 <- numeric()
      if(x$nvar$pred.long>0)
        ind1 <- 1:x$nvar$pred.long
      if(!is.na(match("cure",names(coef))))
        ind1 <- c(ind1,length(coef))
      if(length(ind1)>0){
        message("Long term predictor")
        prmatrix(matrix(tmp[ind1,],nrow=length(ind1)),
                 rowlab=rownames(tmp)[ind1],collab=colnames(tmp))
      }
      ind1 <- setdiff(1:length(coef),ind1)
      if(length(ind1>0)){
        message("\nShort term predictor")
        prmatrix(matrix(tmp[ind1,],nrow=length(ind1)),
                 rowlab=rownames(tmp)[ind1],collab=colnames(tmp))
      }
    }
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) df <- sum(!is.na(coef))
    else  df <- round(sum(x$df),2)
    message()
    message("Likelihood ratio test=", format(round(logtest, digits)), " on ",
        df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
    omit <- x$na.action
    if(length(omit))
      message("n=", x$n, " (", naprint(omit), ")", sep="")
    else message("n=", x$n, sep="")
  }else{
    message("Null model: ",x$call$nlt.model,sep="")
    message("Likelihood=", format(round(x$loglik[1], digits)), " n=", x$n,
        sep="")
  }
  invisible()
}
