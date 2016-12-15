\name{SimplexIP}
\alias{SimplexIP}
\title{Parameter estimation of The Inverse-Power Type}
\description{
  Parameter estimation of the inverse-power type via numerical calculation of the Ripley's K-function.
}
\usage{
  SimplexIP(offspring, pa, ty=1, x2, skip=1, eps=0.1e-2, process=0, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the initial guesses of the parameters \eqn{(\mu, \nu, p, c)}{(mu, nu, p, c)}.}
  \item{ty}{the variable \eqn{T_y}{Ty} for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]}
            without loss of generality except for the scaling.}
  \item{x2}{upper limit value in place of \eqn{\infty}{infty} in the integral in distribution function below.}
  \item{skip}{the variable for the fast likelihood but rough approximation of the initial estimates.
              The \code{skip} calculate the Palm intensity function in the log-likelihood function for every \code{skip}-th \eqn{r_{ij}} 
              in the ordered distances of the pairs \eqn{i} and \eqn{j}.}
  \item{eps}{the optimization procedure is iterated  at most 1000 times until \code{stderr} becomes smaller than \code{eps}.}
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
  \item{logl.palm}{the minimized -log L in the process of minimizing the negative Palm log-likelihood function.}
  \item{mple}{the MPLE (maximum Palm likelihood estimator) values corresponding to the above \code{logl.palm}.}
  \item{logl.simplex}{the minimized -log L by the simplex method.}
  \item{stderr}{the standard deviations among function values of the vertices \eqn{x(i)} of the simplex.}
  \item{pa.normal}{the values of the corresponding normalized parameters \eqn{x(i)}.
                       the last values show the normalized variables corresponding to the MPLE (the maximum Palm likelihood estimator) values relative to the initial estimates.}
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
  mu <- 50;  nu <- 30;  p <- 1.5;  c <- 0.005
  z <- SimulateIP(seeds, c(mu,nu,p,c))

  \dontrun{
    ## estimation
    ## need very long c.p.u time in the minimization procedure
    SimplexIP(z$offspring, c(mu,nu,p,c), x2=0.3, skip=100, process=3)
  }
}

\keyword{spatial}
