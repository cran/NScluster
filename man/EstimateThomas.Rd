\name{EstimateThomas}
\alias{EstimateThomas}
\title{Parameter Estimation of the Thomas Model}
\description{
  Parameter estimation of the Thomas model by using the Palm log-likelihood
 function.}
\usage{
  EstimateThomas(xy.points, pars, eps = 0.001, process.report = 0, plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points in
    a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars}{a named vector of containing the initial guess of the model
    parameters (\code{mu}, \code{nu}, \code{sigma}), where \code{mu} is
    an intensity of parents, \code{nu} is an expected number of descendants for
    each parent and \code{sigma} is a parameter of the dispersal kernel.}
  \item{eps}{the optimization procedure is iterated at most 1000 times until
    \code{process2$stderr} becomes smaller than \code{eps}.}
  \item{process.report}{the level of reporting the process of minimizing.
     Allowed values are as follows:
     \describe{
       \item{0}{no report (default).}
       \item{1}{output the process of minimizing the negative Palm
                log-likelihood function until the values converge to MPLEs.
                (\code{process1})}
       \item{2}{output the process of optimizing by the simplex with
                the normalized parameters. (\code{process2})}
       \item{3}{output both processes.}
     }}
  \item{plot}{logical. If \code{TRUE} (default), the process of optimizing by
    the simplex with the normalized parameters is plotted.}
}

\details{
  The Palm intensity function of the Thomas model is calculated as follows:

  For all \eqn{r \ge 0}{r >= 0},

  \deqn{\lambda_{\bm{o}}(r) = \mu\nu + \frac{\nu}{4\pi \sigma^2}
    \exp \left( -\frac{r^2}{4 \sigma^2} \right).}{%
    \lambda_o(r) = \mu\nu + \nu/(4\pi\sigma^2) * exp(-r^2/(4\sigma^2)).}

  The Palm log-likelihood function of the Thomas model on \eqn{W} is given by

  \deqn{\log L(\mu,\nu,\sigma) = \sum_{\{i,j; i<j, r_{ij} \le 1/2\}}
    \log \nu \left\{ \mu + \frac{1}{4 \pi \sigma^2}
    \exp \left( -\frac{{r_{ij}}^2}{4 \sigma^2} \right) \right\}}{%
    log L(\mu,\nu,\sigma) = \sum_{i, j; i < j, r(i, j) <= 1/2}
    log\nu { \mu + 1/(4\pi\sigma^2) * exp(-r(i, j)^2/(4\sigma^2)) }}

  \deqn{- N(W)\nu \left\{ \frac{\pi \mu}{4} + 1 -
    \exp \left( -\frac{1}{16 \sigma^2} \right) \right\}.}{%
    - N(W) \nu { \pi\mu/4 + 1 - exp(-1/(16\sigma^2)) }.} 
}
\value{
  \item{mple}{MPLE (maximum Palm likelihood estimate).}
  \item{process1}{a list with following components. (Only returned if
    \code{process.report} = 1 or 3.)
    \describe{
      \item{cflg}{1 (="update") or -1 (="testfn"), where "update" indicates that
        -log L value has attained the minimum so far, otherwise not.} 
      \item{logl.palm}{the minimized -log L in the process of minimizing the
        negative Palm log-likelihood function.}
      \item{mples}{corresponding MPLEs.}
    }
  }
  \item{process2}{a list with following components. (Only returned if
    \code{process.report} = 2 or 3.)
    \describe{
      \item{logl.simplex}{the minimized -log L by the simplex method.}
      \item{stderr}{the standard deviations.}
      \item{pa.normal}{the normalized variables corresponding the initial
        estimates.}
    }
  }
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
## simulation
pars <- c(mu = 50.0, nu = 30.0, sigma = 0.03)
z <- SimulateThomas(pars, seed = 117)

## estimation
## need long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu = 40.0, nu = 40.0, sigma = 0.05)
EstimateThomas(z$offspring$xy, init.pars)
}
}

\keyword{spatial}
