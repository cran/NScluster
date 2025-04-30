\name{boot.mple}
\alias{boot.mple}
\alias{summary.boot.mple}
\title{Bootstrap resampling for MPLE}
\description{
  Carry out bootstrap replicates of MPLE on simulated data. 
}
\usage{
boot.mple(mple.out, n = 100, conf.level = 0.95, se = TRUE, trace = FALSE)

\method{summary}{boot.mple}(object, \dots)
}
\arguments{
  \item{mple.out}{an object of class "mple", usually the result of a call to 
    \code{\link{mple.cppm}}.}
  \item{n}{number of bootstrap replicates performed.}
  \item{conf.level}{the confidence level required.}
  \item{se}{logical. If \code{TRUE} standard errors are returned.}
  \item{trace}{logical: if \code{TRUE}, a progress bar is shown.}
  \item{object}{an object of class \code{"boot.mple"}.}
  \item{\dots}{ignored.}
}
\value{
  \code{boot.mple} returns an object of class \code{"boot.mple"} containing the
    following components:
  \item{boot.mples}{a matrix of \code{n} rows each of which is a bootstrap
    replicate of the result of calling \code{\link{mple.cppm}}.}
  \item{confint}{confidence intervals for MPLEs.}
  \item{mple}{MPLE of \code{mple.out} passed as \code{'pars'} argument to
    \code{\link{sim.cppm}}.}
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
