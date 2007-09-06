print.nltm <- function(x, digits=max(options()$digits - 4, 3), ...)
{
  cat(gettextf("Authors: G. Garibotti, A. Tsodikov\n"))

  if(!is.null(cl<- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
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
      cat("Iteration limit maxit: ", x$maxit, " has been reached.\n")
      cat("Consider a new initial value or optimization parameters such as",
          "bscale.\n", sep=" ")
      cat("See nltm.control help.\n")
    }else{
      if(x$convergence==51){
        cat("Warning from optimization method: ", x$message, "\n")
        cat("Consider a new initial value or optimization parameters such as",
            "bscale.\n", sep=" ")
        cat("See nltm.control help.\n")
      }else{
        if(x$convergence==52){
          cat("Error in optimization method: ", x$message, "\n")
          cat("Consider a new initial value or optimization parameters such",
              "as bscale.\n", sep=" ")
          cat("See nltm.control help.\n")
        }
      }
    }

    if(length(negVar)>0){
      cat(gettextf("\nWarning message:\n"))
      cat(gettextf("Problem with covariance matrix of coefficients.\n"))
      cat(gettextf("Some diagonal terms are negative.\n"))
      cat(gettextf("p-values and confidence intervals are unreliable.\n\n\n"))
    }
    
    cat("Non Linear Transformation Model: ", x$call$nlt.model,
        ", fit by maximum likelihood\n\n",sep="")
    tmp <- cbind(coef, exp(coef), se, coef/se,
                 signif(1 - pchisq((coef/ se)^2, 1), digits-1))
    dimnames(tmp) <- list(names(coef),
                          c("coef","exp(coef)","se(coef)","z","p"))
    
    if(x$npred==1){
      cat("\n")
      prmatrix(tmp)
    }else{
      ind1 <- numeric()
      if(x$nvar$pred.long>0)
        ind1 <- 1:x$nvar$pred.long
      if(!is.na(match("cure",names(coef))))
        ind1 <- c(ind1,length(coef))
      if(length(ind1)>0){
        cat("Long term predictor\n")
        prmatrix(matrix(tmp[ind1,],nrow=length(ind1)),
                 rowlab=rownames(tmp)[ind1],collab=colnames(tmp))
      }
      ind1 <- setdiff(1:length(coef),ind1)
      if(length(ind1>0)){
        cat("\nShort term predictor\n")
        prmatrix(matrix(tmp[ind1,],nrow=length(ind1)),
                 rowlab=rownames(tmp)[ind1],collab=colnames(tmp))
      }
    }
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) df <- sum(!is.na(coef))
    else  df <- round(sum(x$df),2)
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, digits)), " on ",
        df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
    cat("\n")
    omit <- x$na.action
    if(length(omit))
      cat("n=", x$n, " (", naprint(omit), ")\n", sep="")
    else cat("n=", x$n, "\n", sep="")
  }else{
    cat("Null model: ",x$call$nlt.model,"\n",sep="")
    cat("Likelihood=", format(round(x$loglik[1], digits)), " n=", x$n, "\n",
        sep="")
  }
  invisible()
}
