\name{PalmTypeB}
\alias{PalmTypeB}
\title{Non-Parametric Estimate of The Palm Intensity of Type B Model}
\description{
  Calculate the non-parametric Palm intensity function of Type B Model estimated directly from a set of point pattern data.}
\usage{
  PalmTypeB(offspring, pa, delta, Ty, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the parameters \eqn{(\mu_i, \nu_i, a_i, \sigma_{1,i}, \sigma_{2,i})}{(mu_i, nu_i, a_i, sigma_{1,i}, sigma_{2,i})}, \eqn{i=1,2,...,I \le 8}{i=1,2,...,I <= 8}.}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{Ty}{the variable for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
  \item{plot}{logical. If TRUE (default) the non-parametric estimate and the curves of the true and MPLE (the maximum Palm likelihood estimator) parameters are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}{r=j*Delta}, where \eqn{j=1,2,...,[R/\Delta]}{j=1,2,...,[R/Delta]},
           where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program for the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity function of \eqn{r}, which is normalized
              by the total intensity estimate of the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions \eqn{\lambda_0(r) / \hat{\lambda}}{lambda_0(r)/hat{lambda}} calculated from the given sets of parameter values
                     \eqn{(\mu_i, \nu_i, a_i, \sigma_{1,i}, \sigma_{2,i})}{(mu_i, nu_i, a_i, sigma_{1,i}, sigma_{2,i})}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  mu1 <- 10; mu2 <- 40; nu <- 30; sig1 <- 0.01; sig2 <- 0.03
  Ty <- 1
  z <- SimulateTypeB( seeds, c(mu1,mu2,nu,sig1,sig2 ), Ty, 100, 150 )

  ## estimation
  delta <- 0.001
  pa1 <- c(50, 30, 0.2, 0.01, 0.03)
  pa2 <- c(0.44100e+02, 0.30691e+02, 0.15119, 0.94494e-02, 0.26822e-01)
  pa <- matrix(c(pa1,pa2), 2, 5, byrow=TRUE) 
  PalmTypeB( z$offspring, pa, delta, Ty )
}

\keyword{spatial}
