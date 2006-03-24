\name{nltm}
\alias{nltm}
\title{
  Fit Non Linear Transformation Model for Analyzing Survival Data 
}
\description{
  Fits a non-linear transformation (NLT) model for analyzing survival
  data, see Tsodikov (2003). The class of NLT models includes the
  following currently supported models Cox proportional hazard and
  proportional hazard cure models, proportional odds model, proportional
  hazard - proportional hazard cure model, proportional hazard -
  proportional odds model, Gamma frailty model, and proportional hazard
  - proportional odds model.
}
\usage{
nltm(formula=formula(data), data=parent.frame(), subset, na.action,
init, control, model=c("PH","PHC","PO","PHPHC","PHPOC","GFM","PHPO"),
verbose=FALSE, ...)
}
\author{
  Gilda Garibotti (garibott AT math.utah.edu), Alexander Tsodikov
}
\arguments{
  \item{formula}{A formula object, with the response on the left of a
    \code{~} operator, and the terms on the right.  The response must be
    a survival object as returned by the \code{Surv} function.}
  \item{data}{A data.frame in which to interpret the variables named in
    the \code{formula}, or in the \code{subset} argument.}
  \item{subset}{Expression saying that only a subset of the rows of the
    data should be used in the fit.}
  \item{na.action}{A missing-data filter function, applied to the
    model.frame, after any subset argument has been used. Default is
    \code{options()$na.action}.}
  \item{init}{Vector of initial values for the calculation of the
    maximum likelihood estimator of the regression parameters. Default
    initial value is zero.}
  \item{control}{Object of class \code{\link{coxph.control}} specifying
    iteration limit and other control options. Default is
    \code{nltm.control(...)}.}
  \item{model}{A character string specifying a non-linear transformation
    model. Default Proportional Hazards Model.

    The conditional survival function \eqn{S(t|z)} given the covariates
    \eqn{z} of each of the models currently supported are given
    below. Let \eqn{S_0(t)} be the non-parametric baseline survival
    function, and \eqn{\theta(z)} and \eqn{\eta(z)} predictors. We take
    \eqn{\theta(z)=\exp(\beta_\theta z)} and
    \eqn{\eta(z)=\exp(\beta_\eta z)}.
    \itemize{
      \item Proportional hazard model (PH):
      \deqn{S(t|z)=S_0(t)^\theta(z).}{S(t|z)=S_0(t)^theta(z).} 
      \item Proportional hazard cure model (PHC):
      \deqn{S(t|z)=\exp(-\theta(z)(1-S_0(t))).}
      \item Proportional odds model (PO):
      \deqn{S(t|z)=\frac{\theta(z)}{\theta(z)-\ln(S_0(t)).}}
      \item Proportional hazard - proportional hazard cure model (PHPHC):
      \deqn{S(t|z)=\exp(-\theta(z)(1-S_0^{\eta(z)}(t))).}
      \item Proportional hazard - proportional odds cure model (PHPOC):
      \deqn{S(t|z)=\exp(-\theta(z)(1-S_0(t))/(1-(1-\eta(z))S_0(t))).}
      \item Gamma frailty model (GFM):
      \deqn{S(t|z)=(\frac{\theta(z)^{\eta(z)}}
	{\theta(z)-\ln(S_0(t))})^\eta(z).}
      \item Proportional hazard - proportional odds model (PHPO):
      \deqn{S(t|z)=\frac {\theta S_0^{\eta(z)}(t)}
	{1-(1-\theta)S_0^{\eta(z)}(t)}.}
    }
  }
  \item{verbose}{If TRUE it stores information from maximization of
    likelihood and calculation of information matrix in a file. Default
    is FALSE.}
  \item{...}{Other arguments}
}
\value{
  an object of class \code{"coxph"}.
}
\references{
  
  Tsodikov AD (2003) "Semiparametric models: a generalized
  self-consistency approach". Journal of the Royal Statistical Society
  B, 65, Part 3, 759-774.

  Tsodikov AD, Ibrahim JG, Yakovlev AY (2003) "Estimating cure rates
  from survival data: an alternative to two-component mixture
  models". Journal of the American Statistical Association, Vol. 98,
  No. 464, 1063-1078.

  Tsodikov AD (2002) "Semi-parametric models of long- and short-term
  survival: an application to the analysis of breast cancer survival in
  Utah by age and stage". Statistics in Medicine, 21, 895-920.

  Wendland MMM, Tsodikov AD, Glenn MJ, Gaffney DK (2004) "Time interval
  to the development of breast carcinoma after treatment
  for Hodgkin disease". Cancer Vol. 101, No. 6, 1275-1282.
  
}
\seealso{
\code{coxph}, \code{nltm.control}
}
\examples{
\dontrun{
# Simple test data set 
test1 <- list(time=c(10,7,32,23,22,6,16,34,32,25,11,20,19,6,17,35,6,13,9,6,1),
              status=c(1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
              size=c(1.79,7.93,2.02,6.89,2.30,7.82,1.25,9.85,6.02,3.43,4.72,7.45,8.83,9.53,1.10,1.06,5.25,5.86,2.03,3.62,3.52),
              age=factor(c(65,65,65,65,99,45,65,99,99,99,65,45,65,55,45,45,55,55,55,99,65)))
nltm(Surv(time,status) ~ size + age, data=test1, model="PO")
}
}
\keyword{survival}
