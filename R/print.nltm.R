print.nltm <-
 function(x, digits=max(options()$digits - 4, 3), ...)
{
  if (!is.null(cl<- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
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
    cat("Consider a new initial value or optimization parameters such as bscale.\n")
    cat("See nltm.control help\n")
  }else{
    if(x$convergence==51){
      cat("Warning from optimization method: ", x$message, "\n")
      cat("Consider a new initial value or optimization parameters such as bscale.\n")
      cat("See nltm.control help\n")
    }else{
      if(x$convergence==52){
        cat("Error in optimization method: ", x$message, "\n")
        cat("Consider a new initial value or optimization parameters such as bscale.\n")
        cat("See nltm.control help\n")
      }
    }
  }

  if(length(negVar)>0){
    cat(gettextf("\nWarning message:\n"))
    cat(gettextf("Problem with covariance matrix of coefficients.\n"))
    cat(gettextf("Some diagonal terms are negative.\n"))
    cat(gettextf("p-values and confidence intervals are unreliable.\n\n\n"))
    }
  
  cat("Non-linear transformation model fit by maximum likelihood\n")
  if (is.null(x$naive.var)) {
    tmp <- cbind(coef, exp(coef), se, coef/se,
                 signif(1 - pchisq((coef/ se)^2, 1), digits -1))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                         "se(coef)", "z", "p"))
  }
  else {
    nse <- sqrt(diag(x$naive.var))
    tmp <- cbind(coef, exp(coef), nse, se, coef/se,
                 signif(1 - pchisq((coef/se)^2, 1), digits -1))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                         "se(coef)", "robust se", "z", "p"))
  }
  cat("\n")
  prmatrix(tmp)
  
  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  if (is.null(x$df)) df <- sum(!is.na(coef))
  else  df <- round(sum(x$df),2)
  cat("\n")
  cat("Likelihood ratio test=", format(round(logtest, 2)), " on ",
      df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
  cat("\n")
  omit <- x$na.action
  if(length(omit))
    cat("n=", x$n, " (", naprint(omit), ")\n", sep="")
  else cat("n=", x$n, "\n")
  
  invisible()
}
