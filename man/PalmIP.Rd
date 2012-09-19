\name{PalmIP}
\alias{PalmIP}
\title{Non-Parametric Estimate of The Palm Intensity of The Inverse-Power Type}
\description{
  Calculate the non-parametric Palm intensity function of the inverse-power type estimated directly from a set of point pattern data.}
\usage{
  PalmIP(offspring, pa, delta, Ty, x2, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the parameters \eqn{(\mu_i, \nu_i, p_i, c_i)}{(mu_i, nu_i, p_i, c_i)}, \eqn{i=1,2,...,I \le 7}{i=1,2,...,I <= 7}.}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{Ty}{the variable for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
  \item{x2}{upper limit value in place of \eqn{\infty}{infty}.}
  \item{plot}{logical. If TRUE (default) the non-parametric estimate and the curves of the true and MPLE (the maximum Palm likelihood estimator) parameters are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}{r=j*Delta}, where \eqn{j=1,2,...,[R/\Delta]}{j=1,2,...,[R/Delta]},
           where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program for the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity function of \eqn{r}, which is normalized
              by the total intensity estimate of the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions \eqn{\lambda_0(r) / \hat{\lambda}}{lambda_0(r)/hat{lambda}} calculated
 from the given sets of parameter values \eqn{(\mu_i, \nu_i, p_i, c_i)}{(mu_i, nu_i, p_i, c_i)}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  pa1 <- c(50, 30, 1.5, 0.005)
  Ty <- 1
  z <- SimulateIP( seeds, pa1, Ty ) 

  ## estimation
  delta <- 0.001
  x2 <- 0.3
  pa2 <- c(0.5651715e+02, 0.2392856e+02, 0.1600531e+01, 0.5565606e-02)
  pa <- matrix(c(pa1,pa2), 2, 4, byrow=TRUE) 
  PalmIP( z$offspring, pa, delta, Ty, x2 )
}

\keyword{spatial}
