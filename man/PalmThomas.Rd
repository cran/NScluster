\name{PalmThomas}
\alias{PalmThomas}
\title{Non-Parametric Estimate of The Palm Intensity of The Thomas Model}
\description{
  Calculate the non-parametric Palm intensity function of the Thomas model estimated directly from a set of point pattern data.}
\usage{
  PalmThomas(offspring, pa, delta, ty, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the parameters \eqn{(\mu_i, \nu_i, \sigma_i)}{(mu_i, nu_i, sigma_i)}, \eqn{i=1,2,...,I \le 8}{i=1,2,...,I <= 8}.}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{ty}{the variable \eqn{T_y}{Ty} for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
  \item{plot}{logical. If TRUE (default) the non-parametric estimate and the curves of the true and MPLE (the maximum Palm likelihood estimator) parameters are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}{r=j*Delta}, where \eqn{j=1,2,...,[R/\Delta]}{j=1,2,...,[R/Delta]},
           where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program for the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity function of \eqn{r}, which is normalized
              by the total intensity estimate of the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions \eqn{\lambda_0(r) / \hat{\lambda}}{lambda_0(r)/hat{lambda}} calculated from the given sets of parameter values
                     \eqn{(\mu_i, \nu_i, \sigma_i)}{(mu_i, nu_i, sigma_i)}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  pa1 <- c(50, 30, 0.03)
  ty <- 1
  z <- SimulateThomas(seeds, pa1, ty)

  ## estimation
  delta <- 0.001
  pa2 <- c(0.51311e+02, 0.26104e+02, 0.29071e-01)
  pa <- matrix(c(pa1,pa2), 2, 3, byrow=TRUE) 
  PalmThomas(z$offspring, pa, delta, ty)
}

\keyword{spatial}
