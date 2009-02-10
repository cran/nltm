update.nltm.formula <- function(fit,form,pred)
{
  if(pred=="long"){
    if(fit$npred==1){
      fit$call$formula1 <- update.formula(fit$formula,form)
    }else{
      fit$call$formula1 <- update.formula(fit$formula$pred.long,form)
      fit$call$formula2 <- update.formula(fit$formula$pred.short,~1)
    }
  }else{
    if(!is.na(match("formula2",names(fit$call))))
      fit$call$formula2 <- update.formula(fit$formula$pred.short,form)
    else
      fit$call$formula2 <- update.formula(fit$formula$pred.long,form)
  }
  fit
}

anova.nltm <- function(object, ..., test=FALSE)
{
  cat(gettextf("Authors: G. Garibotti, A. Tsodikov\n"))
  dotargs <- list(...)
  named <- if (is.null(names(dotargs))) 
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named)) 
    warning(paste("The following arguments to anova.nltm(..)", 
                  "are invalid and dropped:",
                  paste(deparse(dotargs[named]), collapse=", ")))
  dotargs <- dotargs[!named]
  is.nltm <- unlist(lapply(dotargs, function(x) inherits(x,"nltm")))
  dotargs <- dotargs[is.nltm]
  if(length(dotargs)>0)
    return(anova.nltmlist(c(list(object), dotargs), test=test))

  if(object$npred==1){
    if(object$nvar==0)
      stop(paste("\nModel with no covariantes,",
                 "it is not possible to build an anova table.", sep=" "))
  }else{
    if(object$nvar$pred.long==0 & object$nvar$pred.short==0)
      stop(paste("\nModel with no covariantes,",
                 "it is not possible to build an anova table.", sep=" "))
  }
  
  cure <- !is.na(match("cure",names(object$coef)))
  row.name <- character()
  resdev <- resdf <- NULL

  if(object$npred>1){
    if(object$nvar$pred.short>0){
      form <- "~."
      termlist <- attr(object$terms$pred.short,"term.label")
      if(object$nvar$pred.long==0 & !cure)
        termlist <- termlist[-1]
      nvars <- length(termlist)
      if(nvars>0){
        for(i in rev(termlist)){
          form <- paste(form,i,sep="-")
          fit <- update.nltm.formula(object,form,"short")
          fit <- update(fit)
          resdev <- c(resdev, -2*fit$loglik[2])
          resdf <- c(resdf,object$n-sum(!is.na(coef(fit))))
        }
      }
      row.name <- c(paste(attr(object$terms$pred.short,"term.labels"),
                           "(short term)",sep=" "))
    }
    termlist <- attr(object$terms$pred.long,"term.labels")
  }else{
    termlist <- attr(object$terms,"term.labels")
  }
  nvars <- length(termlist)
  form <- "~."
  if(nvars>1){
    for(i in rev(termlist[-1])){
      form <- paste(form,i,sep="-")
      fit <- update.nltm.formula(object,form,"long")
      fit <- update(fit)
      resdev <- c(resdev, -2*fit$loglik[2])
      resdf <- c(resdf,object$n-sum(!is.na(coef(fit))))
    }
    if(cure){
      form <- paste(form,termlist[1],sep="-")
      fit <- update.nltm.formula(object,form,"long")
      fit <- update(fit)
      resdev <- c(resdev, -2*fit$loglik[2])
      resdf <- c(resdf,object$n-sum(!is.na(coef(fit))))
    }
    row.name <- c(paste(termlist, "(long term)", sep=" "),row.name)
  }else{
    if(nvars==1){
      if(cure){
        form <- paste(form,termlist[1],sep="-")
        fit <- update.nltm.formula(object,form,"long")
        fit <- update(fit)
        resdev <- c(resdev, -2*fit$loglik[2])
        resdf <- c(resdf,object$n-sum(!is.na(coef(fit))))
      }
      row.name <- c(paste(termlist, "(long term)", sep=" "),row.name)
    }
  }
  if(cure)
    row.name <- c("cure (long term)",row.name)
  row.name <- c("NULL",row.name)

  resdf <- c(object$n,rev(resdf),object$n-sum(!is.na(coef(object))))
  resdev <- c(-2*object$loglik[1],rev(resdev),-2*object$loglik[2])
  table <- data.frame(c(NA,-diff(resdf)),c(NA,pmax(0,-diff(resdev))),resdf,
                      resdev)
  
  dimnames(table) <- list(row.name,c("Df","Deviance","Resid. Df",
                                      "Resid. Dev"))
  title <- paste("Analysis of Deviance Table\nNLT model: ",
                 object$call$nlt.model, "\nResponse is ",
                 ifelse(object$npred==1,deparse(object$terms[[2]]),
                        deparse(object$terms$pred.long[[2]])),
                 "\nTerms added sequentially (first to last)\n", sep = "")
  dispersion <- 1
  if (test)
    table <- stat.anova(table=table, test="Chisq", scale=dispersion)
  
  structure(table, heading=title, class=c("anova","data.frame"))
}


anova.nltmlist <- function(object, ..., test=FALSE)
{
  responses <- as.character(lapply(object,
                                   function(x) {x$call$formula1[[2]]}))
  same.resp <- responses==responses[1]
  if(!all(same.resp)){
    object <- object[same.resp]
    warning(paste("Models with response", deparse(responses[!same.resp]), 
                  "removed because response differs from model 1."))
  }
  same.n <- sapply(object, function(x) x$n)
  if(any(same.n!=same.n[1]))
    stop("Models were not all fitted to the same size of dataset.")
  nlt.models <- sapply(object, function(x) x$call$nlt.model)
  same.nltm <- nlt.models==nlt.models[1]
  if(!all(same.nltm)){
    object <- object[same.nltm]
    warning(paste("Models with nltm", deparse(responses[!same.nltm]), 
                  "removed because nltm differs from", "model 1"))
  }
  nmodels <- length(object)
  
  if(nmodels==1)
    stop(paste("Only one model left after removing models with different",
               "sample sizes and nltm", sep=" "))
  
  resdf <- as.numeric(lapply(object, function(x) x$n-sum(!is.na(coef(x)))))
  resdev <- as.numeric(lapply(object, function(x) -2*x$loglik[2]))
  table <- data.frame(resdf, resdev, c(NA,-diff(resdf)), c(NA,-diff(resdev)))
  variables1 <- lapply(object, function(x) x$call$formula1)
  dimnames(table) <- list(1:nmodels,
                          c("Resid. Df","Resid. Dev","Df","Deviance"))
  title <- "Analysis of Deviance Table\n"

  if(object[[1]]$npred==1)
    topnote <- paste("Model ", format(1:nmodels), ": ", variables1, 
                     sep="", collapse="\n")
  else{
    num <- character()
    num[seq(1,nmodels*2,2)] <- format(1:nmodels)
    num[seq(2,nmodels*2,2)] <- " "
    variables1[seq(1,nmodels*2,2)] <- variables1
    variables1[seq(2,nmodels*2,2)] <- lapply(object,
                                             function(x) x$call$formula2)
    topnote <- paste(c("Model ","      "), num,
                     c(": Long term:  ","  Short term: "), variables1, sep="",
                     collapse="\n")
  }
  if (test) {
    dispersion <- 1
    table <- stat.anova(table=table, test="Chisq", scale=dispersion)
  }
  structure(table, heading=c(title,topnote), class=c("anova","data.frame"))
}
