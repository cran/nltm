\name{nltm}
\alias{nltm}
\title{
  Fit Non-Linear Transformation Model for analyzing survival data 
}
\description{
  Fits a non-linear transformation model (nltm) for analyzing survival
  data, see Tsodikov (2003). The class of nltm includes the following
  currently supported models: Cox proportional hazard, proportional
  hazard cure, proportional odds, proportional hazard - proportional
  hazard cure, proportional hazard - proportional odds cure, Gamma
  frailty, and proportional hazard - proportional odds. 
}
\usage{
nltm(formula1=formula(data), formula2=formula(data),
     data=parent.frame(), subset, na.action, init=NULL, control,
     nlt.model=c("PH","PHC","PO","PHPHC","PHPOC","GFM","PHPO"),
     model=FALSE, x=FALSE, y=FALSE, verbose=FALSE, ...)
}
\author{
  Gilda Garibotti (garibotti@crub.uncoma.edu.ar) and Alexander Tsodikov.
}
\arguments{
  \item{formula1}{
    A formula object with the response on the left of a \code{~}
    operator, and the terms on the right. The response must be a
    survival object as returned by the \code{Surv} function. In models
    with two predictors, this corresponds to the long term effect.
  }
  \item{formula2}{
    A formula corresponding to the short term effect. Will be ignored in
    models with only one predictor. If not present in models with two
    predictors, then formula1 will be used both for the long and short
    term effect.
  }
  \item{data}{
    A data.frame in which to interpret the variables named in
    \code{formula1} and \code{formula2}, or in the \code{subset}
    argument.
  }
  \item{subset}{
    Expression saying that only a subset of the rows of the data should
    be used in the fit.
  }
  \item{na.action}{
    A missing-data filter function, applied to the model.frame, after
    any subset argument has been used. Default is
    \code{options()$na.action}.
  }
  \item{init}{
    Vector of initial values for the calculation of the maximum
    likelihood estimator of the regression parameters. Default is zero.
  }
  \item{control}{
    Object specifying iteration limit and other control options. Default
    is \code{nltm.control(...)}.
  }
  \item{nlt.model}{
    A character string specifying a non-linear transformation
    model. Default is Proportional Hazards Model.

    The conditional survival function \eqn{S(t|z)} given the covariates
    \eqn{z} of each of the models currently supported are given
    below. Let \eqn{S_0(t)} be the non-parametric baseline survival
    function, and \eqn{\theta(z)} and \eqn{\eta(z)} predictors. We take
    \eqn{\theta(z)=\exp(\beta_\theta z)} and
    \eqn{\eta(z)=\exp(\beta_\eta z)}, where \eqn{\beta_\theta} and
    \eqn{\beta_\eta} are the regresssion coefficients. In cure models,
    there is an additional regression parameter \eqn{\beta_c} and
    \eqn{\theta(z)=\exp(\beta_\theta z+\beta_c)}.
    
    \itemize{
      \item Proportional hazard model (PH):
      \deqn{S(t|z)=S_0(t)^{\theta(z)}.} 
      \item Proportional hazard cure model (PHC):
      \deqn{S(t|z)=\exp\bigl(-\theta(z)(1-S_0(t))\bigr).}{S(t|z)=\exp(-\theta(z)(1-S_0(t))).} 
      \item Proportional odds model (PO):
      \deqn{S(t|z)=\frac{\theta(z)}{\theta(z)-\ln(S_0(t))}.}{S(t|z)=\theta(z)/(\theta(z)-\ln(S_0(t))).} 
      \item Proportional hazard - proportional hazard cure model (PHPHC):
      \deqn{S(t|z)=\exp\bigl(-\theta(z)(1-S_0^{\eta(z)}(t))\bigr).}{S(t|z)=\exp(-\theta(z)(1-S_0^{\eta(z)}(t))).} 
      \item Proportional hazard - proportional odds cure model (PHPOC):
      \deqn{S(t|z)=\exp\biggl(-\frac{\theta(z)(1-S_0(t))}
	{1-(1-\eta(z))S_0(t)}\biggr).}{S(t|z)=\exp(-\theta(z)(1-S_0(t))/(1-(1-\eta(z))S_0(t))).}  
      \item Gamma frailty model (GFM):
      \deqn{S(t|z)=\biggl(\frac{\theta(z)^{\eta(z)}}
	{\theta(z)-\ln(S_0(t))}\biggr)^{\eta(z)}.}{S(t|z)=(\theta(z)^{\eta(z)}/(\theta(z)-\ln(S_0(t))))^{\eta(z)}.} 
      \item Proportional hazard - proportional odds model (PHPO):
      \deqn{S(t|z)=\frac {\theta(z) S_0^{\eta(z)}(t)}
	{1-(1-\theta(z))S_0^{\eta(z)}(t)}.}{S(t|z)=\theta(z) S_0^{\eta(z)}(t)/(1-(1-\theta(z))S_0^{\eta(z)}(t)).}
    }
  }
  \item{model}{If TRUE the model frame is stored. Default is FALSE.} 
  \item{x}{If TRUE the model matrix is stored. Default is FALSE.} 
  \item{y}{If TRUE the response matrix is stored. Default is FALSE.}
  \item{verbose}{If a file name is given, it stores information from
    maximization of likelihood and calculation of information matrix in
    a file. Otherwise, verbose=FALSE. Default is FALSE.}
  \item{\dots}{Other arguments.}
}
\value{
  An object of class \code{"nltm"}. See \code{\link{nltm.object}} for
  details.
}
\references{  
  Tsodikov A. (2003) "Semiparametric models: a generalized
  self-consistency approach". Journal of the Royal Statistical Society
  B, 65, Part 3, 759-774.

  Tsodikov A. (2002) "Semi-parametric models of long- and short-term
  survival: an application to the analysis of breast cancer survival in
  Utah by age and stage". Statistics in Medicine, 21, 895-920.

  Tsodikov A., Garibotti G. (2006) "Profile information matrix for
  nonlinear transformation models". to appear in Journal of Lifetime
  Data Analysis. 

  Tsodikov A., Ibrahim J., Yakovlev A. (2003) "Estimating cure rates
  from survival data: an alternative to two-component mixture
  models". Journal of the American Statistical Association, Vol. 98,
  No. 464, 1063-1078.

  Wendland M., Tsodikov A., Glenn M., Gaffney D. (2004) "Time interval
  to the development of breast carcinoma after treatment
  for Hodgkin disease". Cancer Vol. 101, No. 6, 1275-1282.
}
\seealso{
  \code{\link{nltm.object}}, \code{\link{summary.nltm}},
  \code{\link{nltm.control}}
}
\examples{
# fit a Proportional Odds Model 
data(melanoma)
fit <- nltm(Surv(time,status) ~ size + age, data=melanoma, nlt.model="PO")
summary(fit)
}
\keyword{survival}

