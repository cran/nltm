# y[,1]: times
# y[,2]: status
# x1: regression variables for first predictor (long term)
# x2: regression variables for second predictor (short term)

nltm.fit <- function(x1, x2, y, model, init, control, verbose)
{
  n <-  nrow(y)
  if(is.matrix(x1))
    nvar1 <- ncol(x1)   
  else
    if (length(x1)==0)
      nvar1 <- 0
    else
      nvar1 <- 1

  npred <- nPredictor(model)
  cure <- cureModel(model)

  if(npred>1){
    if(is.null(x2)){
      nbeta <- npred*nvar1+cure
      x2 <- x1
      nvar2 <- ncol(x2)
    }else{
      nvar2 <- ncol(x2)
      nbeta <- nvar1+nvar2+cure
    }
  }else{
    nbeta <- nvar1+cure
    x2 <- 0
    nvar2 <- 0
  }

  # Eliminate censored observations prior to first event time.
  # Individuals with censored time before the first event time are not
  # used in the estimation
  sorted <- order(y[,1])
  nc.small <- length(which(y[,1]<min(y[y[,2]==1,1])))+1
  # All observations are censored
  if(nc.small>n)
    stop("All observations are censored. The model cannot be fitted.\n")
  sorted <- sorted[nc.small:n]
  
  y <- y[sorted,]

  if(nvar1>0){
    names <- dimnames(x1)
    names[[1]] <- names[[1]][sorted]
    x1 <- matrix(x1[sorted,], nrow=length(sorted), dimnames=names)
  }
  if(npred>1 & nvar2>0){
    names <- dimnames(x2)
    names[[1]] <- names[[1]][sorted]
    x2 <- matrix(x2[sorted,], nrow=length(sorted), dimnames=names)
  }
  time <- eventTimes(y)  # time and status correspond to stime and sstat in 
  status <- y[,2]        # coxph.fit
    
  count <- counts(time, status)
  isurv <- initSurvival(count, cure)
  s0 <- isurv$s0
  
  if(!is.null(init)){
    if(length(init)!=nbeta)
      stop("Wrong length for initial values")
  }else{
    init <- rep(0,nbeta)
    if(cure)
      init[nbeta] <- log(-log(isurv$tailDefect))
  }

  if(nbeta>0){
    bound <- boundary(x1, x2, npred, cure, control$bscale)
    
    fit <- optim(par=init, fn=profileLikR, gr=NULL, method="L-BFGS-B",
                 lower=-bound, upper=bound,
                 control=control[-match(c("bscale","s0.tol"),names(control))],
                 hessian=FALSE, x1, x2, status, count, s0, model, cure,
                 control$s0.tol, nvar1, nvar2, as.integer(n-nc.small+1), npred,
                 verbose)
    
    coef <- fit$par
    
    if(npred==1)
      names(coef) <- dimnames(x1)[[2]]
    else
      names(coef) <- c(dimnames(x1)[[2]],dimnames(x2)[[2]])
    if(cure)
      names(coef)[nbeta] <- "cure"

    # Find survival function at betaMLE
    survMLE <- .C("profileLik", coef=as.double(coef), t(x1), t(x2),
                  as.integer(status), count$dd, count$rr, surv=as.double(s0),
                  model, as.integer(cure), control$s0.tol, nvar1, nvar2,
                  nrow(count), as.integer(n-nc.small+1), as.integer(npred),
                  as.integer(verbose), plik=double(1), PACKAGE="nltm")$surv
  
    # Find covariance matrix
    imat1 <- .C("informationMatrix", coef, t(x1), t(x2), as.integer(status),
                count$dd, count$rr, survMLE, model, as.integer(cure), nvar1,
                nvar2, nrow(count), as.integer(n-nc.small+1),
                as.integer(npred), as.integer(verbose),
                imat=double(nbeta*nbeta), PACKAGE="nltm")$imat

    # Invert information matrix
    cov <- solve(matrix(imat1, nbeta, nbeta))

    counts <- fit$count[1]
    names(counts) <- NULL
    res <- list(coefficients=coef, surv=survMLE, var=cov, n=n,
                maxit=control$maxit, counts=counts,
                convergence=fit$convergence, message=fit$message)
  }
  # Find profile likelihood of null model
  init <- rep(0,nbeta)
  lik0 <- .C("profileLik", coef=as.double(init), t(x1), t(x2),
             as.integer(status), count$dd, count$rr, surv=as.double(s0),
             model, as.integer(cure), control$s0.tol, nvar1, nvar2,
             nrow(count), as.integer(n-nc.small+1), as.integer(npred),
             as.integer(verbose), plik=double(1), PACKAGE="nltm")$plik

  if(nbeta>0)
    res$loglik <- c(lik0, fit$value)
  else
    res <- list(loglik=lik0, n=n)
  res
}


# The actual value of the censored times doesn't matter, people with
# censored time in [t_i, t_{i+1}) get t_i
eventTimes <- function(y)
{
  td <- unique(y[y[,2]==1,1])  # td is sorted because y[,1] is
  time <- array(NA, nrow(y)) 
  
  for(i in 1 : length(td))
    time[y[,1]>=td[i]] <- td[i]
  time
}

# count$dd: number of deaths at t_i
# count$rr: number of dead or censored people at t_i
counts <- function(time, status)
{
  count <- data.frame(table(time[status==1]),table(time))
  count <- count[,c(2,4)]
  colnames(count) <- c("dd","rr")
  count
}

reverseCumsum <- function(a)
{
  res <- cumsum(a[length(a):1])[length(a):1]
  res
}


# use the Nelson-Aalen etimator (Klein, Moeschberger, page 86) to
# initialize the survival jumps
# Note: cumprod(res) is the KM in [SP]
initSurvival <- function(count, cure)
{
  if(cure){
    nt <- nrow(count)
    aux <- cumsum(count$dd/reverseCumsum(count$rr))
    res <- aux[nt]-aux
    res <- list(s0=res/c(aux[nt],res[1:(nt-1)]), tailDefect=exp(-aux[nt]))
  }else
    res <- list(s0=exp(-count$dd/reverseCumsum(count$rr)))
  res
}

# compute boundaries for maximization
# Note: this function assumes that either length(bscale)=1 or
# bscale=(bscale_theta, bscale_eta, bscale_cure); that is, long-term
# predictor, short-term predictor, cure term
boundary <- function(x1, x2, npred, cure, bscale)
{
  r <- numeric()
  if(ncol(x1)>0)
    r <- apply(abs(x1), MARGIN=2, max)
  if(npred>1)
    r <- c(r, apply(abs(x2), MARGIN=2, max))
  r <- ifelse(r<1e-10, 1e-10, r)
  
  bound <- bscale/r
  
  if(cure)
    bound <- c(bound, ifelse(length(bscale)>1, bscale[length(bscale)], bscale))
  bound
}


# control parameters for nltm
nltm.control <- function(fnscale=-1, maxit=1000, reltol, factr=1e7, pgtol=0,
                         s0.tol=1e-5, bscale=5)
{
  # control parameters for optimization (see optim help)
  if(fnscale>0) warning("Minimization will take place")
  if(maxit<=0) stop("Invalid value for iterations")
  if(!missing(reltol)) if(reltol<0) stop("Invalid convergence tolerance")
  if(factr<0) stop("Invalid factor convergence tolerance")
  if(pgtol<0) stop("Invalid tolerance on the projected gradient")

  if(s0.tol<0) stop("Invalid convergence tolerance of hazard self-consistency")
  if(sum(bscale<=0)>0) stop("Invalid bound scale")
  
  ifelse(missing(reltol), 
         control <- list(fnscale=fnscale, maxit=maxit, factr=factr,
                         pgtol=pgtol, s0.tol=s0.tol, bscale=bscale),
         control <- list(fnscale=fnscale, maxit=maxit, reltol=reltol,
                         factr=factr, pgtol=pgtol, s0.tol=s0.tol,
                         bscale=bscale))
  control
}


profileLikR <- function(beta, x1, x2, status, count, s0, model, cure, tol,
                        nvar1, nvar2, nobs, npred, verbose)
{
  # need to transpose x otherwise when passed to C it is stored in a vector
  # by columns and I need it by rows because of dmat
  res <- .C("profileLik", coef=as.double(beta), t(x1), t(x2),
            as.integer(status), count$dd, count$rr, surv=as.double(s0), model,
            as.integer(cure), tol, nvar1, nvar2, nrow(count), nobs,
            as.integer(npred), as.integer(verbose), plik=double(1),
            PACKAGE="nltm")$plik 

  res
}

nPredictor <- function(model)
{
  switch(model, PH=1, PHC=1, PO=1, PHPHC=2, PHPOC=2, GFM=2, PHPO=2)
}

cureModel <- function(model)
{
  switch(model, PH=FALSE, PHC=TRUE, PO=FALSE, PHPHC=TRUE, PHPOC=TRUE,
         GFM=FALSE, PHPO=FALSE)
}
