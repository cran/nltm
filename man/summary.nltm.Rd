\name{summary.nltm}
\alias{summary.nltm}
\title{
  Summary of a nltm objects. 
}
\description{
  This function finds confidence intervals for relative
  risks and calculates the log-likelihood test. It is a method for the
  generic function \code{summary} of class \code{nltm}. It can be
  invoked by calling summary for an object of class nltm.
}
\usage{
\method{summary}{nltm}(object,  coef = TRUE, conf.int = 0.95, 
        digits = max(options()$digits - 4, 3),...)
}
\arguments{
  \item{object}{
    Fitted model object of class \code{nltm}.  This is assumed to be the
    result of some function that produces an object with the same named
    components as that returned by the \code{nltm} function. 
  }
  \item{coef}{
    If true it gives a table with coefficients, standard errors and
    p-values. Default is TRUE.
  }
  \item{conf.int}{
    Confidence level of confidence intervals for relative risks. If 0
    confidence intervals are not computed. Default is 0.95.
  }
  \item{digits}{
    Number of digits used for formatting output.
  }
  \item{\dots}{
    Arguments to be passed to or from other methods.
  }
}
\seealso{
  \code{\link{nltm}}, \code{\link{nltm.object}}, \code{\link{summary}}.
}
\keyword{survival}

