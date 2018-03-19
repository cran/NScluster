\name{PalmIP}
\alias{PalmIP}
\title{Non-Parametric and Parametric Estimate of the Palm Intensity of
  Inverse-Power Type Model}
\description{
  Calculate the non-parametric and parametric Palm intensity function of
  the Inverse-power type model estimated directly from a set of point pattern
  data.}
\usage{
  PalmIP(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001, uplimit = 0.3,
         plot=TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points
    in a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars1}{a named vector of true parameters (\code{mu}, \code{nu}, \code{p},
    \code{c}), where \code{mu} is an intensity of parents, \code{nu} is
    an expected number of descendants for each parent, \code{p} is the decay
    order and \code{c} is the scaling parameter.}
  \item{pars2}{a named vector of MPLEs (the maximum Palm likelihood estimates)
    (\code{mu}, \code{nu}, \code{p}, \code{c}).}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{uplimit}{upper limit value in place of \eqn{\infty}.}
  \item{plot}{logical. If \code{TRUE} (default), the non-parametric estimate and
    the curves of true parameters and MPLEs are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta},
    where \eqn{j=1,2,\dots,[R/\Delta]}{j=1,2,...,[R/\Delta]}, where [ ] is
    the Gauss' symbol and \eqn{R=1/2} is given in the program for the normalized
    rectangular region for the point pattern.}
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
pars <- c(mu = 50.0, nu = 30.0, p = 1.5, c = 0.005)
z <- SimulateIP(pars, seed = 353)

## estimation
## need very long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu = 55.0, nu = 35.0, p = 1.0, c = 0.01)
z1 <- EstimateIP(z$offspring$xy, init.pars, skip = 100)
# Parameter               mu            nu             p             c
# Initial value      55.0000       35.0000        1.0000        0.0100
# MPLE               48.3143       30.5047        1.6467        0.0079
}

## Palm intensity
par1 <- c(50.0, 30.0, 1.5, 0.005)  # pars
par2 <- c(48.3143, 30.5047, 1.6467, 0.0079) # z1$mple
PalmIP(z$offspring$xy, par1, par2)
}

\keyword{spatial}
