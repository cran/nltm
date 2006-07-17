summary.nltm <- function(object,  coef = TRUE, conf.int = 0.95, 
                         digits = max(options()$digits - 4, 3),...)
{
  beta <- object$coef
  nabeta <- !(is.na(beta))          #non-missing coefs
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
                            paste("lower .",round(100*conf.int, 2),sep = ""),
                            paste("upper .",round(100*conf.int, 2),sep = "")))
    rval$conf.int <- tmp
  }
  df <- length(beta2)
  logtest <- -2 * (object$loglik[1] - object$loglik[2])
  rval$logtest <- c(test=logtest, df=df, pvalue=1 - pchisq(logtest, df))
  class(rval) <- "summary.nltm"
  rval
}


print.summary.nltm <- function(x, digits = max(options()$digits - 4, 3), ...)
{
  if(!is.null(x$call)) {
    cat("Call:\n")
    dput(x$call)
    cat("\n")
  }

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

  negVar <- which(diag(x$var)<=0)
  if(length(negVar)>0){
    cat(gettextf("\nWarning message:\n"))
    cat(gettextf("Problem with covariance matrix of coefficients.\n"))
    cat(gettextf("Some diagonal terms are negative.\n"))
    cat(gettextf("p-values and confidence intervals are unreliable.\n\n\n"))
    }
  
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Non-linear transformation model fit by maximum likelihood\n\n")
  if (nrow(x$coef)==0) {   # Null model
    cat ("Null model\n")
    return()
  }
  
  if(!is.null(x$coef)) {
    prmatrix(x$coef)
  }
  if(!is.null(x$conf.int)) {
    cat("\n")
    prmatrix(x$conf.int)
  }
  cat("\n")
  cat("Likelihood ratio test=", format(round(x$logtest["test"], 2)),
      " on ", x$logtest["df"], " df, p=", format(x$logtest["pvalue"]),
      sep="")
  cat("\n")
  
  omit <- x$na.action
  if(length(omit))
    cat("n=", x$n, " (", naprint(omit), ")\n", sep="")
  else cat("n=", x$n, "\n")
  
  invisible()
}
