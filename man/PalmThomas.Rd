\name{PalmThomas}
\alias{PalmThomas}
\title{Non-Parametric and Parametric Estimate of the Palm Intensity of Thomas
  Model}
\description{
  Calculate the non-parametric and parametric Palm intensity function of the
  Thomas model estimated directly from a set of point pattern data.}
\usage{
  PalmThomas(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001, plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points
    in a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars1}{a named vector of the true parameters (\code{mu}, \code{nu},
    \code{sigma}), where \code{mu} is an intensity of parents, \code{nu} is
    an expected number of descendants for each parent and \code{sigma} is
    a parameter of the dispersal kernel.}
  \item{pars2}{a named vector of MPLEs (the maximum Palm likelihood estimates)
    (\code{mu}, \code{nu}, \code{sigma}).}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{plot}{logical. If \code{TRUE} (default), the non-parametric estimate and
    the curves of true parameters and MPLEs are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}, where
    \eqn{j=1,2,\dots,[R/\Delta]}{j=1,2,...,[R/\Delta]},
    where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program for
    the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity
    function of \eqn{r}, which is normalized by the total intensity estimate of
    the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions
    \eqn{\lambda_{\bm{o}}(r)/\hat{\lambda}}{\lambda_o(r)/\lambda^}
    calculated from the given sets of parameter values.
    See 'Details' in \code{\link{EstimateThomas}}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
## simulation
pars <- c(mu = 50.0, nu = 30.0, sigma = 0.03)
z <- SimulateThomas(pars, seed = 117)

## estimation
## need long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu = 40.0, nu = 40.0, sigma = 0.05)
z1 <- EstimateThomas(z$offspring$xy, init.pars)
# Parameter               mu            nu         sigma
# Initial value      40.0000       40.0000        0.0500
# MPLE               46.7657       27.2079        0.0279
}


## Palm intensity
par1 <- c(50.0, 30.0, 0.03)  # pars
par2 <- c(46.7657, 27.2079, 0.0279)  # z1$mple
PalmThomas(z$offspring$xy, par1, par2)
}

\keyword{spatial}
