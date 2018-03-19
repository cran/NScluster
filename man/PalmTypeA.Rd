\name{PalmTypeA}
\alias{PalmTypeA}
\title{Non-Parametric and Parametric Estimate of the Palm Intensity of Type A
  Model}
\description{
  Calculate the non-parametric and parametric Palm intensity function of Type A
  model estimated directly from a set of point pattern data.}
\usage{
  PalmTypeA(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001, uplimit = 0.3,
            plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points
    in a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars1}{a named vector of the true parameters (\code{mu}, \code{nu},
    \code{a}, \code{sigma1}, \code{sigma2}), where \code{mu} is an intensity of
    parents, \code{nu} is an expected number of descendants for each parent,
    \code{a} is a mixture parameter, \code{sigma1} and \code{sigma2} are
    parameters of the dispersal kernel for each component.}
  \item{pars2}{a named vector of MPLEs (the maximum Palm likelihood estimates)
    (\code{mu}, \code{nu}, \code{a}, \code{sigma1}, \code{sigma2}).}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{uplimit}{upper limit value in place of \eqn{\infty}.}
  \item{plot}{logical. If \code{TRUE} (default), the non-parametric estimate and
    the curves of true parameters and MPLEs are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}, where
    \eqn{j=1,2,\dots,[R/\Delta]}{j=1,2,...,[R/\Delta]},
    where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program
    for the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity
    function of \eqn{r}, which is normalized by the total intensity estimate of
    the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions
    \eqn{\lambda_{\bm{o}}(r)/\hat{\lambda}}{\lambda_o(r)/\lambda^}
    calculated from the given sets of parameter values.
    See \url{../doc/NScluster-guide.pdf}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
## simulation
pars <- c(mu = 50.0, nu = 30.0, a = 0.3, sigma1 = 0.005, sigma2 = 0.1)
z <- SimulateTypeA(pars, seed=575)

## estimation
## need very long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu=60.0, nu=40.0, a=0.5, sigma1=0.01, sigma2=0.1)
z1 <- EstimateTypeA(z$offspring$xy, init.pars, skip=100)
# Parameter              mu           nu            a       sigma1       sigma2
# Initial value     60.0000      40.0000       0.5000       0.0100       0.1000
# MPLE              51.2441      25.1439       0.3431       0.0054       0.0824
}

## Palm intensity
par1 <- c(50.0, 30.0, 0.3, 0.005, 0.1)  # pars
par2 <- c(51.2441, 25.1439, 0.3431, 0.0054, 0.0824)  # z1$mple
PalmTypeA(z$offspring$xy, par1, par2)
}

\keyword{spatial}
