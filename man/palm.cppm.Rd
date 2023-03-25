\name{palm.cppm}
\alias{palm.cppm}
\alias{print.Palm}
\title{Non-parametric and Parametric Estimation for Palm Intensity}
\description{
  Compute the non-parametric and the parametric Palm intensity function of the
  Neyman-Scott cluster point process models and their extensions.
}
\usage{
palm.cppm(mple, pars = NULL, delta = 0.001, uplimit = 0.3)

\method{print}{Palm}(x, \dots)
}
\arguments{
  \item{mple}{an object of class "mple".}
  \item{pars}{a named vector of the true parameters, if any.}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{uplimit}{upper limit in place of \eqn{\infty} of the integral in the
    probability distribution function relative to the random distance between
    two descendant points within the same cluster. The \code{uplimit} is valid 
    for \code{"IP"} and \code{"TypeA"}.}
  \item{x}{an object of class \code{"Palm"}.}
  \item{\dots}{ignored.}
}
\value{
  An object of class \code{"Palm"} containing the following components:
  \item{r}{the distance \eqn{r=j\Delta},
    where \eqn{j=1,2,\dots,[R/\Delta]}, \eqn{[\cdot]}{[.]} is
    the Gauss' symbol and \eqn{R=1/2} in the program.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity
    function, which is normalized by the total intensity estimate (the mean
    number of points in \eqn{W}) of a given point pattern data.}
  \item{norm.palm}{the corresponding values of the normalized Palm intensity
    function, i.e., \eqn{\lambda_{\bm{o}}(r)/\hat{\lambda}}{%
    \lambda_o(r)/\lambda^}, where \eqn{\lambda_{\bm{o}}(r)}{\lambda_o(r)} is the
    Palm intensity and \eqn{\lambda} is an intensity of a cluster point process
    model. See 'Details' in \code{\link{mple.cppm}.}}
  There is another method \code{\link{plot.Palm}} for this class.
}
\references{
Tanaka, U., Ogata, Y. and Katsura, K. (2008)
 Simulation and estimation of the Neyman-Scott type spatial cluster models.
 \emph{Computer Science Monographs} \bold{34}, 1-44.
 The Institute of Statistical Mathematics, Tokyo.
 \url{https://www.ism.ac.jp/editsec/csm/}.
}
\seealso{
See \code{\link{sim.cppm}} and \code{\link{mple.cppm}} to simulate the
 Neyman-Scott cluster point process models and their extensions and to compute
 the MPLEs, respectively.
}
\examples{
\dontrun{
# The computation of MPLEs takes a long CPU time in the minimization procedure,
# especially for the Inverse-power type and the Type A models.

### Thomas Model
 #simulation
 pars <- c(mu = 50.0, nu = 30.0, sigma = 0.03)
 t.sim <- sim.cppm("Thomas", pars, seed = 117)
 ## estimation => Palm intensity
 init.pars <- c(mu = 40.0, nu = 40.0, sigma = 0.05)
 t.mple <- mple.cppm("Thomas", t.sim$offspring$xy, init.pars)
 t.palm <- palm.cppm(t.mple, pars)
 plot(t.palm)

### Inverse-Power Type Model
 # simulation
 pars <- c(mu = 50.0, nu = 30.0, p = 1.5, c = 0.005)
 ip.sim <- sim.cppm("IP", pars, seed = 353)
 ## estimation => Palm intensity
 init.pars <- c(mu = 55.0, nu = 35.0, p = 1.0, c = 0.01)
 ip.mple <- mple.cppm("IP", ip.sim$offspring$xy, init.pars, skip = 100)
 ip.palm <- palm.cppm(ip.mple, pars)
 plot(ip.palm)

### Type A Model
# simulation
 pars <- c(mu = 50.0, nu = 30.0, a = 0.3, sigma1 = 0.005, sigma2 = 0.1)
 a.sim <- sim.cppm("TypeA", pars, seed=575)
 ## estimation => Palm intensity
 init.pars <- c(mu=60.0, nu=40.0, a=0.5, sigma1=0.01, sigma2=0.1)
 a.mple <- mple.cppm("TypeA", a.sim$offspring$xy, init.pars, skip=100)
 a.palm <- palm.cppm(a.mple, pars)
 plot(a.palm)

### Type B Model
 # simulation
 pars <- c(mu1 = 10.0, mu2 = 40.0, nu = 30.0, sigma1 = 0.01, sigma2 = 0.03)
 b.sim <- sim.cppm("TypeB", pars, seed = 257)
 ## estimation => Palm intensity
 init.pars <- c(mu1 = 20.0, mu2 = 30.0, nu = 30.0, sigma1 = 0.02, sigma2 = 0.02)
 b.mple <- mple.cppm("TypeB", b.sim$offspring$xy, init.pars)
 b.palm <- palm.cppm(b.mple, pars)
 plot(b.palm)


### Type C Model
 # simulation
 pars <- c(mu1 = 5.0, mu2 = 9.0, nu1 = 30.0, nu2 = 150.0,
           sigma1 = 0.01, sigma2 = 0.05)
 c.sim <- sim.cppm("TypeC", pars, seed = 555)
 ## estimation => Palm intensity
 init.pars <- c(mu1 = 10.0, mu2 = 10.0, nu1 = 30.0, nu2 = 120.0,
                sigma1 = 0.03, sigma2 = 0.03)
 c.mple <- mple.cppm("TypeC", c.sim$offspring$xy, init.pars)
 c.palm <- palm.cppm(c.mple, pars)
 plot(c.palm)
}
}

\keyword{spatial}
