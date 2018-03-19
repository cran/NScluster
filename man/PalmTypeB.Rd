\name{PalmTypeB}
\alias{PalmTypeB}
\title{Non-Parametric and Parametric Estimate of the Palm Intensity of Type B
  Model}
\description{
  Calculate the non-parametric and parametric Palm intensity function of Type B
  model estimated directly from a set of point pattern data.}
\usage{
  PalmTypeB(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001, plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points
    in a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars1}{a named vector of the true parameters (\code{mu1}, \code{mu2},
    \code{nu}, \code{sigma1}, \code{sigma2}), where (\code{mu}\eqn{i}, \code{nu},
    \code{sigma}\eqn{i}) is an intensity of parents, an expected number of
    descendants, a parameter of the dispersal kernel for superposed component
    \eqn{i} (\eqn{i = 1,2}), respectively.}
  \item{pars2}{a named vector of MPLEs (the maximum Palm likelihood estimates)
    (\code{mu1}, \code{mu2}, \code{nu}, \code{sigma1}, \code{sigma2}).}
  \item{delta}{a width for the non-parametric Palm intensity function.}
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
    See 'Details' in \code{\link{EstimateTypeB}}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
## simulation
pars <- c(mu1 = 10.0, mu2 = 40.0, nu = 30.0, sigma1 = 0.01, sigma2 = 0.03)
z <- SimulateTypeB(pars, seed = 257)

## estimation
## need very long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu1 = 20.0, mu2 = 30.0, nu = 30.0, sigma1 = 0.02, sigma2 = 0.02)
z1 <- EstimateTypeB(z$offspring$xy, init.pars)
# Parameter             mu1          mu2          nu        sigma1       sigma2
# Initial value     20.0000      30.0000      30.0000       0.0200       0.0200
# MPLE              16.1778      44.3974      28.3942       0.0101       0.0312
}

## Palm intensity
par1 <- c(10.0, 40.0, 30.0, 0.01, 0.03) # pars
par2 <- c(16.1778, 44.3974, 28.3942, 0.0101, 0.0312) # z1$mple
PalmTypeB(z$offspring$xy, par1, par2)
}

\keyword{spatial}
