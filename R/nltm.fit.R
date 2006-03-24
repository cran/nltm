# y[,1]: times
# y[,2]: status
# x: regression variables

nltm.fit <- function(x, y, model, init, control, verbose)
{
#  cat(gettextf("problems if x is not a matrix x <- x[sorted,] would not work"))

  n <-  nrow(y)
  if (is.matrix(x))
    nvar <- ncol(x)   
  else
    if (length(x)==0)
      nvar <- 0
    else
      nvar <- 1

  npred <- nPredictor(model)
  cure <- cureModel(model)
  nbeta <- npred*nvar+cure
  

### THIS NEEDS TO BE CHECKED
  if(nvar==0){
    # A special case: Null model.
    #  (This is why I need the rownames arg- can't use x' names)
    # Set things up for 0 iterations on a dummy variable
    x <- as.matrix(rep(1.0, n))
    nvar <- 1
    init <- 0
    maxiter <- 0
    nbeta <- npred
  }
###

  # Eliminate censored observations prior to first event time.
  # Individuals with censored time before the first event time are not
  # used in the estimation
  sorted <- order(y[,1])
  nc.small <- length(which(y[,1]<min(y[y[,2]==1,1])))+1
  ### See!!!! What happens if all observations are censored
  sorted <- sorted[nc.small:n]
  
  y <- y[sorted,]
  x <- x[sorted,]
  time <- eventTimes(y)  # time and status correspond to stime and sstat in 
  status <- y[,2]        # coxph.fit

  count <- counts(time, status)
  isurv <- initSurvival(count, cure)
  s0 <- isurv$s0
#  print("s0 init"); print(s0)

  if(nvar!=0){
    if (!missing(init) && !is.null(init)) {
      if (length(init) != nbeta) stop("Wrong length for inital values")
    }
    else{
      init <- rep(0,nbeta)
      if(cure)
        init[nbeta] <- log(-log(isurv$tailDefect))
    }
  }    
  bound <- boundary(x, npred, cure)
  
  fit <- optim(par=init, fn=profileLikR, gr=NULL, method="L-BFGS-B",
               lower=-bound, upper=bound, control=control,
               hessian=FALSE, x, status, count, s0, model, cure,
               control$s0.tol, npred, verbose)
#  print("s0 init"); print(s0)
#  print("opt value"); print(fit$value)

  coef <- fit$par
  names(coef) <- rep(dimnames(x)[[2]], npred)
  if(cure)
    names(coef)[nbeta] <- "cure"
  
  # Find profile likelihood at initial values
  lik0 <- .C("profileLik", coef=as.double(init), t(x), as.integer(status),
             count$dd, count$rr, surv=as.double(s0), model, as.integer(cure),
             control$s0.tol, ncol(x), nrow(count), nrow(x), as.integer(npred), 
             as.integer(verbose), plik=double(1), PACKAGE="nltm")$plik
  
  # Find survival function at betaMLE
  survMLE <- .C("profileLik", coef=as.double(coef), t(x), as.integer(status),
                count$dd, count$rr, surv=as.double(s0), model,
                as.integer(cure), control$s0.tol, ncol(x), nrow(count),
                nrow(x), as.integer(npred), as.integer(verbose),
                plik=double(1), PACKAGE="nltm")$surv
 
  # Find covariance matrix
  imat1 <- .C("informationMatrix", coef, t(x), as.integer(status), count$dd,
              count$rr, survMLE, model, as.integer(cure), ncol(x), nrow(count),
              nrow(x), as.integer(npred), as.integer(verbose),
              imat=double(nbeta*nbeta), PACKAGE="nltm")$imat

#  print("imat"); print(imat1)
#  print("det(imat)"); print(det(matrix(imat1, nbeta, nbeta)))
  # Invert information matrix
  cov <- solve(matrix(imat1, nbeta, nbeta))
#  print("cov mat"); print(cov)
#  print("1/det(cov)"); print(1/det(cov))
  
  list(coefficients=coef, loglik=c(lik0, fit$value), surv=survMLE,
       iter=fit$counts, var=cov, n=n)
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
boundary <- function(x, npred, cure)
{
  r <- apply(abs(x), MARGIN=2, max)
  r <- ifelse(r<1e-10, 1e-10, r)
  bound <- rep(5/r, npred)
  if(cure)
    bound <- c(bound, 5)
  bound
}


# control parameters for nltm
nltm.control <- function(fnscale=-1, maxit=1000, reltol, factr=1e7, pgtol=0,
                         s0.tol=1e-5)
{
  # control parameters for optimization (see optim help)
  # fnscale=-1 for maximization, default is minimization
  if(fnscale>0) warning("Minimization will take place")
  if(maxit<=0) stop("Invalid value for iterations")
  if(!missing(reltol)) if(reltol<0) stop("Invalid convergence tolerance")
  if(factr<0) stop("Invalid factor convergence tolerance")
  if(pgtol<0) stop("Invalid tolerance on the projected gradient")

  if(s0.tol<0) stop("Invalid convergence tolerance of hazard self-consistency")
  
  ifelse(missing(reltol), 
         control <- list(fnscale=fnscale, maxit=maxit, factr=factr,
                         pgtol=pgtol, s0.tol=s0.tol),
         control <- list(fnscale=fnscale, maxit=maxit, reltol=reltol,
                         factr=factr, pgtol=pgtol, s0.tol=s0.tol))
  control
}


profileLikR <- function(beta, x, status, count, s0, model, cure, tol, npred,
                        verbose)
{
#  print("profileLikR")
#  print(status)
  # need to transpose x otherwise when passed to C it is stored in a vector
  # by columns and I need it by rows because of dmat
  res <- .C("profileLik", coef=as.double(beta), t(x), as.integer(status),
            count$dd, count$rr, surv=as.double(s0), model, as.integer(cure),
            tol, ncol(x), nrow(count), nrow(x), as.integer(npred),
            as.integer(verbose), plik=double(1), PACKAGE="nltm")$plik 

#  print(c(beta, res))
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
