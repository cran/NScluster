\name{SimplexThomas}
\alias{SimplexThomas}
\title{Parameter estimation of The Thomas Model}
\description{
  Parameter estimation of the Thomas model by using the Palm log-likelihood function.}
\usage{
  SimplexThomas(offspring, pa, Ty=1, eps=0.1e-2, process=0, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the initial guesses of the parameters \eqn{(\mu, \nu, \sigma)}{(mu, nu, sigma)}.}
  \item{Ty}{the variable for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
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

\details{
  The Palm intensity function of the Thomas model is calculated as follows:

  For any \eqn{r \ge 0}{r >= 0},

  \deqn{\lambda_0(r) = \mu\nu + \frac{\nu}{4\pi \sigma^2} \exp \left( -\frac{r^2}{4 \sigma^2} \right).}{%
        lambda_0(r) = mu*nu + nu/(4 \pi sigma^2) exp(-r^2/(4sigma^2)).}

  The Palm log-likelihood function of the Thomas model is analytically calculated as follows:

  \deqn{\log L(\mu,\nu,\sigma) = \sum_{\{i,j; i \ne j, r_{ij}<R\}} \log \nu \left\{ \mu + \frac{1}{4 \pi \sigma^2} \exp \left( -\frac{r_{ij}^2}{4 \sigma^2} \right) \right\}}{%
        log L(mu,nu,sigma) = sum_{i,j; i<>j, r_{ij}<R} log nu{mu + 1/(4 \pi sigma^2) exp(-r_{ij}^2/(4sigma^2))}}

  \deqn{- N(W)\nu \left\{ \pi \mu R^2 + 1 - \exp \left( -\frac{R^2}{4 \sigma^2} \right) \right\},}{%
        - N(W)nu{\pi mu R^2 + 1 - exp(-R^2/(4sigma^2))},}

  with \eqn{R = 1/2} which means the half of the \eqn{t_x} (TX) in side length of the normalized rectangular. 
}
\value{
  \item{logL.p}{the minimized -log L in the process of minimizing the negative Palm log-likelihood function.}
  \item{mple}{the MPLE (maximum Palm likelihood estimator) values corresponding to the above \code{logL.p}.}
  \item{logL.s}{the minimized -log L by the simplex method.}
  \item{stderr}{the standard deviation.}
  \item{pa.normal}{the normalized variables corresponding to the MPLE values relative to the initial estimates.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  mu <- 50;  nu <- 30;  sig <- 0.03
  z <- SimulateThomas( seeds, c(mu,nu,sig) )

  \dontrun{
    ## estimation
    ## need long c.p.u time in the minimization procedure
    SimplexThomas( z$offspring, c(mu,nu,sig), process=3 )
  }
}

\keyword{spatial}
