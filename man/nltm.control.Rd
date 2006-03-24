\name{nltm.control}
\alias{nltm.control}
\title{Package options for nltm}
\description{
  Sets default values for parametes related to optimization and update
  of baseline hazard
}
\usage{
nltm.control(fnscale=-1, maxit=1000, reltol, factr=1e7, pgtol=0,
s0.tol=1e-5) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fnscale}{An overall scaling to be applied to the value of 'fn'
    and 'gr' during optimization. If negative, turns the problem into a
    maximization problem. Optimization is performed on 'fn(par)/fnscale'.}
  \item{maxit}{The maximum number of iterations. Defaults to 1000.}
  \item{reltol}{Relative convergence tolerance.  The algorithm stops if
    it is unable to reduce the value by a factor of 'reltol * (abs(val)
    + reltol)' at a step.  Defaults to sqrt(.Machine\$double.eps),
    typically about 1e-8.}
  \item{factr}{Controls the convergence of the '"L-BFGS-B"'
    method. Convergence occurs when the reduction in the objective is
    within this factor of the machine tolerance. Default is 1e7, that is
    a tolerance of about 1e-8.}
  \item{pgtol}{Helps controls the convergence of the '"L-BFGS-B"'
    method. It is a tolerance on the projected gradient in the current
    search direction. This defaults to zero, when the check is
    suppressed.}
  \item{s0.tol}{Convergence tolerance of hazard self-consistency
    equation. Default 1e-5.}
}
\value{
  A list with the same elements as the input
}

\seealso{\code{\link{optim}} }

\examples{

}
\keyword{survival}
