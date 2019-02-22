\name{boot.mple}
\alias{boot.mple}
\alias{summary.boot.mple}
\title{Bootstrap resampling for MPLE}
\description{
  Carry out bootstrap replicates of MPLE on simulated data. 
}
\usage{
boot.mple(mple.out, n = 100, conf.level = 0.95, se = TRUE)

\method{summary}{boot.mple}(object, \dots)
}
\arguments{
  \item{mple.out}{an object of class "mple", usually the result of a call to 
    \code{\link{mple.cppm}}.}
  \item{n}{number of bootstrap replicates performed.}
  \item{conf.level}{the confidence level required.}
  \item{se}{logical. If \code{TRUE} standard errors are returned.}
  \item{object}{an object of class \code{"boot.mple"}.}
  \item{\dots}{ignored.}
}
\value{
  \code{boot.mple} returns an object of class \code{"boot.mple"} containing the
    following components:
  \item{mple}{the values of each parameter (= \code{mple.out$mple}) passed as
    arguments to function \code{\link{sim.cppm}}.}
  \item{boot.mples}{a matrix of bootstrap replicates of MPLEs.}
  \item{confint}{confidence intervals for MPLEs.}
}
\examples{
### Thomas Model
# simulation
pars <- c(mu = 50.0, nu = 30.0, sigma = 0.03)
t.sim <- sim.cppm("Thomas", pars, seed = 117)

\dontrun{ # estimation (need long CPU time)
init.pars <- c(mu = 40.0, nu = 40.0, sigma = 0.05)
t.mple <- mple.cppm("Thomas", t.sim$offspring$xy, init.pars)
t.boot <- boot.mple(t.mple)
summary(t.boot)
}
}
\keyword{spatial}
