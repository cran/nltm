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
\alias{cureModel}
\title{Internal non-linear transformation model functions}
\description{Internal non-linear transformation model functions.}
\usage{
nltm.fit(x1, x2, y, model, init, control, verbose)
profileLikR(beta, x1, x2, status, count, s0, model, cure, tol, nvar1,
nvar2, nobs, npred, verbose)
boundary(x1, x2, npred, cure, bscale)
counts(time, status)
eventTimes(y)
initSurvival(count, cure)
nPredictor(model)
reverseCumsum(a)
cureModel(model)
}
\details{These are not to be called by the user.}
\keyword{internal}
