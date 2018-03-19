\name{PalmTypeC}
\alias{PalmTypeC}
\title{Non-Parametric and Parametric Estimate of the Palm Intensity of Type C
  Model}
\description{
  Calculate the non-parametric and parametric Palm intensity function of Type C
  model estimated directly from a set of point pattern data.}
\usage{
  PalmTypeC(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001, plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points
    in a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars1}{a named vector of the true parameters (\code{mu1}, \code{mu2},
    \code{nu1}, \code{nu2}, \code{sigma1}, \code{sigma2}), where
    (\code{mu}\eqn{i}, \code{nu}\eqn{i}, \code{sigma}\eqn{i}) is an intensity of
    parents, an expected number of descendants, a parameter of the dispersal
    kernel for superposed component \eqn{i} (\eqn{i = 1,2}), respectively.}
  \item{pars2}{a named vector of MPLEs (the maximum Palm likelihood estimates)
    (\code{mu1}, \code{mu2}, \code{nu1}, \code{nu2}, \code{sigma1},
    \code{sigma2}).}
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
    \eqn{\lambda_{\bm{o}}(r)/\hat{\lambda}}{\lambda_o(r)/\lambda^} calculated from
    the given sets of parameter values.
    See 'Details' in \code{\link{EstimateTypeC}}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
## simulation
pars <- c(mu1 = 5.0, mu2 = 9.0, nu1 = 30.0, nu2 = 150.0,
          sigma1 = 0.01, sigma2 = 0.05)
z <- SimulateTypeC(pars, seed = 555)

## estimation
## need long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu1 = 10.0, mu2 = 10.0, nu1 = 30.0, nu2 = 120.0,
               sigma1 = 0.03, sigma2 = 0.03)
z1 <- EstimateTypeC(z$offspring$xy, init.pars)
## Parameter          mu1        mu2        nu1        nu2     sigma1     sigma2
## Initial value  10.0000    10.0000    30.0000   120.0200     0.0300     0.0300
## MPLE            4.6481     5.1690    28.1250   296.7022     0.0088     0.0477
}

## Palm intensity
par1 <- c(5.0, 9.0, 30.0, 150.0, 0.01, 0.05)  # pars
par2 <- c(4.6481, 5.1690, 28.1250, 296.7022, 0.0088, 0.0477)  # z1$mple
PalmTypeC(z$offspring$xy, par1, par2)
}

\keyword{spatial}
