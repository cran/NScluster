\name{PalmTypeA}
\alias{PalmTypeA}
\title{Non-Parametric Estimate of The Palm Intensity of Type A Model}
\description{
  Calculate the non-parametric Palm intensity function of Type A Model estimated directly from a set of point pattern data.}
\usage{
  PalmTypeA(offspring, pa, delta, ty, x2, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the parameters \eqn{(\mu_i, \nu_i, a_i, \sigma_{1,i}, \sigma_{2,i})}{(mu_i, nu_i, a_i, sigma_{1,i}, sigma_{2,i})}, \eqn{i=1,2,...,I \le 7}{i=1,2,...,I <= 7}.}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{ty}{the variable \eqn{T_y}{Ty} for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
  \item{x2}{upper limit value in place of \eqn{\infty}.}
  \item{plot}{logical. If TRUE (default) the non-parametric estimate and the curves of the true and MPLE (the maximum Palm likelihood estimator) parameters are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}{r=j*Delta}, where \eqn{j=1,2,...,[R/\Delta]}{j=1,2,...,[R/Delta]},
           where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program for the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity function of \eqn{r}, which is normalized
              by the total intensity estimate of the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions \eqn{\lambda_0(r) / \hat{\lambda}}{lambda_0(r)/hat{lambda}} calculated
 from the given sets of parameter values \eqn{(\mu_i, \nu_i, a_i, \sigma_{1,i}, \sigma_{2,i})}{(mu_i, nu_i, a_i, sigma_{1,i}, sigma_{2,i})}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  pa1 <- c(50, 30, 0.3, 0.005, 0.1)
  ty <- 1
  z <- SimulateTypeA(seeds, pa1, ty, 100, 150)

  ## estimation
  delta <- 0.001
  x2 <- 0.3
  pa2 <- c(0.5639245e+02, 0.2356894e+02, 0.3549922, 0.5203741e-02, 0.1070969)
  pa <- matrix(c(pa1,pa2), 2, 5, byrow=TRUE) 
  PalmTypeA(z$offspring, pa, delta, ty, x2)
}

\keyword{spatial}
