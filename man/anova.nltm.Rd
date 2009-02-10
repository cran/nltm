\name{anova.nltm}
\alias{anova.nltm}
\alias{anova.nltmlist}
\alias{update.nltm.formula}
\title{Analysis of deviance table for NLT models}
\description{
  If a single model is specified, gives a sequential analysis of
  deviance table for that fit.

  Each row of the table represents a model, the first corresponds to the
  NULL model and the following rows correspond to models in which
  variables have been added sequentially. If a model has both, long term
  and short term covariates, long term covariates are added first. The
  last row corresponds to the input model (object). 

  Each row of the table gives the likelihood ratio test for comparing
  the model corresponding to that row and the one without the indicated
  covariate (Deviance) and its degrees of freedom (Df). It also gives -2
  times the loglikelihood for the model (Resid. Dev) and  its degrees of
  freedom (Resid. Df). Optionally it gives the p-value for the
  likelihood ratio test. 

  If more than one object is specified, the table has a row for each
  model. It specifies -2 times the loglikelihood for the model
  (Resid. Dev) and its degrees of freedom (Resid. Df). It also gives the
  likelihood ratio test for comparing the model corresponding to the
  given row with that of the previous row. (This only makes statistical
  sense if the models are nested.) Optionally it gives the p-value for
  the likelihood ratio test. 
}
\usage{
\method{anova}{nltm}(object, \dots, test=FALSE)
}

\arguments{
  \item{object}{An object of class \code{nltm}.}
  \item{\dots}{Further \code{nltm} objects.}
  \item{test}{If TRUE the p-value for the likelihood ratio test is
    computed. Default is FALSE.}
}
\value{
  An object of class \code{"anova"} inheriting from class
  \code{"data.frame"} with the information given in the description.
}
\section{Warning}{
  The comparison between two or more models by \code{anova} will only be
  valid if they are fitted to the same dataset. This may be a problem if
  there are missing values.}
\seealso{
  \code{\link{nltm}}, \code{anova}.
}
\examples{
data(melanoma)
fit1 <- nltm(formula1=Surv(time,status) ~ size + stage, 
formula2=Surv(time,status) ~ size, data=melanoma, nlt.model="PHPHC")
anova(fit1)
anova(fit1,test=TRUE)
fit2 <- nltm(formula1=Surv(time,status) ~ size + stage + age, 
formula2=Surv(time,status) ~ size + age, data=melanoma, nlt.model="PHPHC")
anova(fit1,fit2,test=TRUE)
 }
\keyword{survival}
