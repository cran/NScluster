\name{SimplexTypeA}
\alias{SimplexTypeA}
\title{Parameter estimation of Type A model}
\description{
  Parameter estimation of type A model via numerical calculation of the Ripley's K-function.
}
\usage{
  SimplexTypeA(offspring, pa, Ty=1, x2, skip=1, eps=0.1e-2, process=0, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the initial guesses of the parameters \eqn{(\mu, \nu, a, \sigma_1, \sigma_2)}{(mu, nu, a, sigma_1, sigma_2)}.}
  \item{Ty}{the variable for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]}
            without loss of generality except for the scaling.}
  \item{x2}{upper limit value in place of \eqn{\infty}{infty} in the integral in the distribution function.}
  \item{skip}{the variable for the fast likelihood but rough approximation of the initial estimates.
              The \code{skip} calculate the Palm intensity function in the log-likelihood unction for every \code{skip}-th \eqn{r_{ij}}{r_ij} 
              in the ordered distances of the pairs \eqn{i} and \eqn{j}.}
  \item{eps}{the optimization procedure is iterated at most 1000 times until \code{stderr} becomes smaller than \code{eps}.}
  \item{process}{repot the process of minimizing. Allowed values are
    \tabular{rl}{
      0 : \tab no report.\cr
      1 : \tab output the process of minimizing the negative Palm log-likelihood function until the values converge to\cr
          \tab the MPLE values for given data.\cr
      2 : \tab output the process of optimizing by the simplex with the normalized parameters depending on \code{pa}.\cr
          \tab The actual estimates are obtained by the indicated x-values times \code{pa}.\cr
      3 : \tab output the both processes.
    }}
  \item{plot}{plot the process of optimizing by the simplex with the normalized parameters depending on \code{pa}.}
}
\value{
  \item{logL.p}{the minimized -log L in the process of minimizing the negative Palm log-likelihood function.}
  \item{mple}{the MPLE (maximum Palm likelihood estimator) values corresponding to the above \code{logL.p}.}
  \item{logL.s}{the minimized -log L by the simplex method.}
  \item{stderr}{the standard deviations.}
  \item{pa.normal}{the normalized variables corresponding to the MPLE values relative to the initial estimates.}
}
\note{
  \code{SimplexIP} and \code{SimplexTypeA} have to use numerical integration and difference to compute the functions, which need very long c.p.u. time in the minimization procedure.
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  mu <- 50;  nu <- 30;  a <- 0.3;  sig1 <- 0.005;  sig2 <- 0.1
  z <- SimulateTypeA( seeds, c(mu,nu,a,sig1,sig2), pmax=100, omax=150 )

  \dontrun{
    ## estimation
    ## need very long c.p.u time in the minimization procedure
    SimplexTypeA( z$offspring, c(mu,nu,a,sig1,sig2), x2=0.3, skip=1000, process=3 )
  }
}

\keyword{spatial}
