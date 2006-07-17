\name{nltm.object}
\alias{nltm.object}
\alias{print.nltm}
\title{
  Non-Linear Transformation Model Object
}
\description{
  This object is returned by the \code{nltm} function. It represents a
  fitted non-linear transformation model. 

  Objects of this class have methods for the functions \code{print} and
  \code{summary}. 
}
\section{Components}{
  Components of a  \code{nltm} object.
  \describe{
    \item{coefficients}{
      The maximum profile likelihood estimators of the model regression
      parameters, \eqn{\beta_\theta} (and \eqn{\beta_\eta} and
      \eqn{\beta_c}). See \code{\link{nltm}}.
    }
    \item{loglik}{
      A vector of length 2 containing the log-likelihood with the initial
      values and with the final values of the coefficients.
    }
    \item{surv}{
      MLE of the baseline survival function at the profile maximum
      likelihood parameters. It is obtained from the hazard jumps that
      satisfy self-consistency equation (5) in Tsodikov A. and Garibotti
      G. (2006).
    }
    \item{var}{
      The variance matrix of the coefficients.
    }
    \item{n}{
      The number of observations used in the fit.
    }
    \item{maxit}{
      The maximum number of iterations of the optimization
      procedure. Default is 1000.
    }
    \item{counts}{
      Number of calls to the profile likelihood function during the
      optimization process. This excludes those calls needed to compute a
      finite-difference approximation to the gradient.
    }
    \item{convergence}{
      An integer code. '0' indicates successful convergence.
      Error codes are
      
      '1' indicates that the iteration limit 'maxit' had been reached.
      
      '51' indicates a warning from the optimization method; see component
      'message' for further details.
      
      '52' indicates an error from the "L-BFGS-B" method; see component
      'message' for further details.
    }
    \item{message}{
      A character string giving any additional information returned
      by the optimizer, or 'NULL'.
    }
    \item{formula}{
      If the model has only one predictor a single formula indicating the
      model for that predictor. If the model has two predictors, then
      formula is a list with terms predictor1 and predictor2 indicating
      the model for the long and short term predictor respectively. 
    }
    \item{call}{
      The call of the nltm model.
    }
    \item{na.action}{
      The \code{na.action} attribute, if any, that was returned by the
      \code{na.action} routine.
    }
  }
}
\references{
  Tsodikov A., Garibotti G. (2006) "Profile information matrix for
  nonlinear transformation models". to appear in Journal of Lifetime
  Data Analysis. 
}
\seealso{
  \code{\link{nltm}}, \code{\link{summary.nltm}}. 
}
\keyword{survival}


