\name{nltm.control}
\alias{nltm.control}
\title{Package options for nltm}
\description{
  Sets default values for arguments related to calculation of the
  maximum profile likelihood estimator of the regression parameters,
  \eqn{\beta_\theta} (and \eqn{\beta_\eta}, and \eqn{\beta_c}), and
  the baseline hazard \eqn{S_0(t)} (see
  \code{\link{nltm}}). Optimization is performed using the "L-BFGS-B"
  method by Byrd et. al. (1995). See \code{optim}.
}
\usage{
nltm.control(fnscale=-1, maxit=1000, reltol, factr=1e7, pgtol=0,
             s0.tol=1e-5, bscale=5) 
}
\arguments{
  \item{fnscale}{
    An overall scaling to be applied to the profile likelihood function
    (profileLik) during optimization. If positive, turns the problem
    into a minimization problem. Optimization is performed on
    'profileLik(par)/fnscale'. Default is -1.
  }
  \item{maxit}{
    The maximum number of iterations. Default is 1000.
  }
  \item{reltol}{
    Relative convergence tolerance.  The algorithm stops if it is unable
    to reduce the value by a factor of 'reltol * (abs(val) + reltol)' at
    a step. Default is sqrt(.Machine\$double.eps), typically about 1e-8.
  }
  \item{factr}{
    Controls the convergence of the "L-BFGS-B" method. Convergence
    occurs when the reduction in the objective is within this factor of
    the machine tolerance. Default is 1e7, that is a tolerance of about
    1e-8. 
  }
  \item{pgtol}{
    Helps control the convergence of the "L-BFGS-B" method. It is a
    tolerance on the projected gradient in the current search
    direction. Default is 0.
  }
  \item{s0.tol}{
    Convergence tolerance of baseline hazard self-consistency
    equation. Default is 1e-5.
  }
  \item{bscale}{
    The maximum profile likelihood estimator is obtained by maximizing
    the profile likelihood over a region determined by the magnitude of
    the observed covariates. These constraints are imposed in order to
    avoid numerical problems in the calculation of the profile
    likelihood function.
    
    For a given regression parameter, corresponding to a covariate with
    observed values x, the upper bound is the bscale of the parameter
    divided by max(abs(x)) if max(abs(x))>1e-10, otherwise 1e-10. The
    lower bound is minus the upper bound. 

    Different values of bscale are allowed for different parameters. If
    different values of bscale are provided, the vector needs to have the
    scale for the regression parameters of the long-term predictor,
    \eqn{\beta_\theta}, followed by the scale for the short-term
    predictor, \eqn{\beta_\eta}, if available and the one for the cure
    parameter, \eqn{\beta_c}, last if a cure model. 

    Default is 5.
  }
}
\value{
  A list with the same elements as the input.
}
\references{
  Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995) A limited memory
  algorithm for bound constrained optimization. SIAM J. Scientific
  Computing, 16, 1190-1208. 
}

\seealso{
  \code{optim}, \code{\link{nltm}}.
}
\keyword{survival}
