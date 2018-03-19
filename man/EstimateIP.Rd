\name{EstimateIP}
\alias{EstimateIP}
\title{Parameter Estimation of the Inverse-Power Type Model}
\description{
  Parameter estimation of the Inverse-power type model via numerical calculation
  of the Ripley's \eqn{K}-function.
}
\usage{
  EstimateIP(xy.points, pars, eps = 0.001, uplimit = 0.3, skip = 1,
             process.report = 0, plot = TRUE)
}
\arguments{
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points in
    a unit square: \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars}{a named vector of containing the initial guess of the model
 parameters (\code{mu}, \code{nu}, \code{p}, \code{c}), where \code{mu} is an
   intensity of parents, \code{nu} is an expected number of descendants for each
   parent, \code{p} is the decay order and \code{c} is the scaling parameter.}
  \item{eps}{the optimization procedure is iterated  at most 1000 times until
    \code{process2$stderr} becomes smaller than \code{eps}.}
  \item{uplimit}{upper limit value in place of \eqn{\infty} in the integral in
    distribution function.}
  \item{skip}{the variable for the fast likelihood but rough approximation of
    the initial estimates. The \code{skip} calculate the Palm intensity function
    in the log-likelihood function for every \code{skip}-th \eqn{r_{ij}}{r(i,j)}
    in the ordered distances of the pairs \eqn{i} and \eqn{j}.}
  \item{process.report}{the level of reporting the process of minimizing.
    Allowed values are as follows:
    \describe{
      \item{0}{no report (default).}
      \item{1}{output the process of minimizing the negative Palm
        log-likelihood function until the values converge to MPLEs.
        (\code{process1})}
        \item{2}{output the process of optimizing by the simplex with the
          normalized parameters. (\code{process2})}
        \item{3}{output both processes.}
    }
  }
  \item{plot}{logical. If \code{TRUE} (default), the process of optimizing by
    the simplex with the normalized parameters is plotted.}
}
\value{
  \item{mple}{MPLE (maximum Palm likelihood estimate).}
  \item{process1}{a list with following components.
    (Only returned if \code{process.report} = 1 or 3.)
    \describe{
      \item{cflg}{1 (="update") or -1 (="testfn"), where "update" indicates that
        -log L value has attained the minimum so far, otherwise not.} 
      \item{logl.palm}{the minimized -log L in the process of minimizing the
        negative Palm log-likelihood function.}
      \item{mples}{corresponding MPLEs.}
    }
  }
  \item{process2}{a list with following components.
    (Only returned if \code{process.report} = 2 or 3.)
    \describe{
      \item{logl.simplex}{the minimized -log L by the simplex method.}
      \item{stderr}{the standard deviations.}
      \item{pa.normal}{the normalized variables corresponding the initial
        estimates.}
    }
  }
}
\note{
  \code{EstimateIP} and \code{EstimateTypeA} have to use numerical integration
  and difference to compute the functions, which need very long c.p.u. time
  in the minimization procedure.
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
## simulation
pars <- c(mu = 50.0, nu = 30.0, p = 1.5, c = 0.005)
z <- SimulateIP(pars, seed = 353)

## estimation
## need very long c.p.u time in the minimization procedure
\dontrun{
init.pars <- c(mu = 55.0, nu = 35.0, p = 1.0, c = 0.01)
EstimateIP(z$offspring$xy, init.pars, skip = 100)
}
}

\keyword{spatial}
