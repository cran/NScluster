\name{PalmTypeC}
\alias{PalmTypeC}
\title{Non-Parametric Estimate of The Palm Intensity of Type C Model}
\description{
  Calculate the non-parametric Palm intensity function of Type C Model estimated directly from a set of point pattern data.}
\usage{
  PalmTypeC(offspring, pa, delta, ty, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the parameters \eqn{(\lambda_i, \nu_{1,i}, a_i, \sigma_{1,i}, \sigma_{2,i})}{(lambda_i, nu_{1,i}, a_i, sigma_{1,i}, sigma_{2,i})}, \eqn{i=1,2,...,I \le 8}{i=1,2,...,I<=8}.}
  \item{delta}{a width for the non-parametric Palm intensity function.}
  \item{ty}{the variable \eqn{T_y}{Ty} for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
  \item{plot}{logical. If TRUE (default) the non-parametric estimate and the curves of the true and MPLE (the maximum Palm likelihood estimator) parameters are shown.}
}
\value{
  \item{r}{the distance \eqn{r=j\Delta}{r=j*Delta}, where \eqn{j=1,2,...,[R/\Delta]}{j=1,2,...,[R/Delta]},
           where [ ] is the Gauss' symbol and \eqn{R=1/2} is given in the program for the normalized rectangular region for the point pattern.}
  \item{np.palm}{the corresponding values of the non-parametric Palm intensity function of \eqn{r}, which is normalized
              by the total intensity estimate of the point pattern data.}
  \item{palm.normal}{the normalized Palm intensity functions \eqn{\lambda_0(r) / \hat{\lambda}}{lambda_0(r)/hat{lambda}}
 calculated from the given sets of parameter values \eqn{(\lambda_i, \nu_{1,i}, a_i, \sigma_{1,i}, \sigma_{2,i})}{(lambda_i, nu_{1,i}, a_i, sigma_{1,i}, sigma_{2,i})}.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  mu1 <- 5;  nu1 <- 30;  sig1 <- 0.01 
  mu2 <- 9;  nu2 <- 150;  sig2 <- 0.05
  ty <- 1
  z <- SimulateTypeC(seeds, c(mu1,nu1,sig1), c(mu2,nu2,sig2), ty, 200, 300)

  ## estimation
  delta <- 0.001
  pa1 <- c(1500, 30, 0.1, 0.01, 0.05)
  pa2 <- c(0.1281e+04, 0.43134e+02, 0.30049, 0.14926e-01, 0.51719e-01)
  pa <- matrix(c(pa1,pa2), 2, 5, byrow=TRUE) 
  PalmTypeC(z$offspring, pa, delta, ty)
}

\keyword{spatial}
