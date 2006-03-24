\name{nltm-internal}
\alias{nltm-internal}
\alias{nltm.fit}
\alias{profileLikR}
\alias{boundary}
\alias{counts}
\alias{eventTimes}
\alias{initSurvival}
\alias{nPredictor}
\alias{reverseCumsum}
\title{Internal non linear transformation model functions}
\description{Internal non linear transformation model functions}
\usage{
nltm.fit(x, y, model, init, control, verbose)
profileLikR <- function(beta, x, status, count, s0, model, cure, tol,
npred)
boundary <- function(x, npred, cure)
counts <- function(time, status)
eventTimes(y)
initSurvival <- function(count, cure)
nPredictor(model)
reverseCumsum(a)
}
\details{These are not to be called by the user}
\keyword{survival}
\keyword{internal}
