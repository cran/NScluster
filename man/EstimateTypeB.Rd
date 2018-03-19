\name{EstimateTypeB}
\alias{EstimateTypeB}
\title{Parameter Estimation of the Type B Model}
\description{
  Parameter estimation of the Type B model by using the Palm Log-Likelihood
  Function.
}
\usage{
  EstimateTypeB(xy.points, pars, eps = 0.001, process.report = 0, plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points in
    a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars}{a named vector of containing the initial guess of the model
    parameters (\code{mu1}, \code{mu2}, \code{nu}, \code{sigma1}, \code{sigma2}),
    where (\code{mu}\eqn{i}, \code{nu}, \code{sigma}\eqn{i}) is an intensity of
    parents, an expected number of descendants, a parameter of the dispersal
    kernel for superposed component \eqn{i} (\eqn{i = 1,2}), respectively.}
  \item{eps}{the optimization procedure is iterated at most 1000 times until
    \code{process2$stderr} becomes smaller than \code{eps}.}
  \item{process.report}{the level of reporting the process of minimizing.
    Allowed values are as follows:
    \describe{
      \item{0}{no report (default).}
      \item{1}{output the process of minimizing the negative Palm log-likelihood
        function until the values converge to MPLEs. (\code{process1})}
      \item{2}{output the process of optimizing by the simplex with
        the normalized parameters. (\code{process2})}
      \item{3}{output both processes.}
    }
  }
  \item{plot}{logical. If \code{TRUE} (default), the process of optimizing by
    the simplex with the normalized parameters is plotted.}
}
\details{
  The Palm intensity function of the Type B model is calculated as follows:

  For all \eqn{r \ge 0}{r >= 0},

  \deqn{\lambda_{\bm{o}}(r) = \lambda + \frac{\nu}{4 \pi} 
    \left\{ \frac{a}{{\sigma_1}^2} \exp \left( -\frac{r^2}{4{\sigma_1}^2}
    \right)+ \frac{(1-a)}{{\sigma_2}^2} \exp \left( -\frac{r^2}{4{\sigma_2}^2}
    \right) \right\},}{%
    \lambda_o(r) = \lambda + \nu/(4\pi) * { a/\sigma1^2 * exp(-r^2/(4\sigma1^2))
    + (1-a)/\sigma2^2 * exp(-r^2/(4\sigma2^2)) },}

  where \eqn{\lambda = \nu(\mu_1+\mu_2)}{\lambda = \nu(\mu1+\mu2)} is the total
  population size and \eqn{a = \mu_1/(\mu_1+\mu_2)}{a = \mu1/(\mu1+\mu2)} is
  the ratio of the parent points of the smaller sized cluster to the total ones.

  The Palm log-likelihood function of the Type B model on \eqn{W} is given by

  \eqn{\log L(\lambda, \alpha, \beta, \sigma_1, \sigma_2)}{log L(\lambda, \alpha,
 \beta, \sigma1, \sigma2)}
  \deqn{=\sum_{\{i,j; i<j, r_{ij} \le 1/2\}} \log
    \left[ \lambda + \frac{1}{4 \pi} \left\{ \frac{\alpha}{{\sigma_1}^2} \exp
    \left( -\frac{{r_{ij}}^2}{4{\sigma_1}^2} \right) + \frac{\beta}{{\sigma_2}^2}
    \exp \left( -\frac{{r_{ij}}^2}{4{\sigma_2}^2} \right) \right\} \right]}{%
    =\sum_{i, j; i < j, r(i, j) <= 1/2} log[ \lambda + 1/(4\pi) { \alpha/\sigma1^2 *
    exp(-r(i, j)^2/(4\sigma1^2)) + \beta/\sigma2^2 * exp(-r(i, j)^2/(4\sigma2^2))
    } ]}
  \deqn{- N(W) \left[ \frac{\pi \lambda}{4} + \alpha \left\{ 1 - \exp
    \left( -\frac{1}{16{\sigma_1}^2} \right) \right\} + \beta \left\{ 1- \exp
    \left( -\frac{1}{16{\sigma_2}^2} \right) \right\} \right],}{%
    - N(W) [ \pi\lambda/4 + \alpha{ 1-exp(-1/(16\sigma1^2)) } +
    \beta{ 1-exp(-1/(16\sigma2^2)) } ],}

  where \eqn{\alpha = a\nu} and \eqn{\beta = (1-a)\nu}. 
}
\value{
  \item{mple}{MPLE (maximum Palm likelihood estimate).}
  \item{process1}{a list with following components.
    (Only returned if \code{process.report} = 1 or 3.)
    \describe{
      \item{cflg}{1 (="update") or -1 (="testfn"), where "update" indicates that
        -log L value has attained the minimum so far, otherwise not.} 
      \item{logl.palm}{the minimized -log L in the process of minimizing
        the negative Palm log-likelihood function.}
      \item{mples}{corresponding MPLEs (\code{mu}, \code{nu}, \code{a},
        \code{sigma1}, \code{sigma2}), where \code{mu} = \code{mu1}+\code{mu2}
        and \code{a} = \code{mu1}/(\code{mu1}+\code{mu2}).}
    }
  }
  \item{process2}{a list with following components.
    (Only returned if \code{process.report} = 2 or 3.)
    \describe{
      \item{logl.simplex}{the minimized -log L by the simplex method.}
      \item{stderr}{the standard deviations.}
      \item{pa.normal}{the normalized variables (\code{mu}, \code{nu}, \code{a},
        \code{sigma1}, \code{sigma2}) as described above.}
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
pars <- c(mu1 = 10.0, mu2 = 40.0, nu = 30.0, sigma1 = 0.01, sigma2 = 0.03)
z <- SimulateTypeB(pars, seed = 257)

## estimation
## need very long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu1 = 20.0, mu2 = 30.0, nu = 30.0, sigma1 = 0.02, sigma2 = 0.02)
EstimateTypeB(z$offspring$xy, init.pars)
}
}

\keyword{spatial}
