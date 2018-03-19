\name{SimulateThomas}
\alias{SimulateThomas}
\title{Simulation of the Thomas Model}
\description{
  Simulation of the Thomas model.
}
\usage{
  SimulateThomas(pars, seed = NULL, plot = TRUE)
}
\arguments{
  \item{pars}{a named vector of containing the values of the model parameters
    (\code{mu}, \code{nu}, \code{sigma}), where \code{mu} is an intensity of
    parents, \code{nu} is an expected number of descendants for each parent and
    \code{sigma} is a parameter of the dispersal kernel.}
  \item{seed}{a positive integer, which is the seed for a sequence of uniform
    random numbers. The default seed is based on the current time.}
  \item{plot}{logical. If \code{TRUE} (default), simulated parent points and
    offspring points are plotted.}
}
\details{
  Let random variable \eqn{U} be independently and uniformly distributed in [0,1].
  We put

  \deqn{U = \int_0^r q_\sigma(t)dt = 1 -
                                  \exp \left( -\frac{r^2}{2\sigma^2} \right).}{%
        U = integral_0^r q_\sigma(t)dt = 1 - exp(-r^2/(2\sigma^2)).}

  Then we have

  \deqn{r = \sigma \sqrt{-2 \log(1-U)}.}{%
        r = \sigma \sqrt{ -2log(1-U) }.}

  Each of the offspring coordinates \eqn{(x_j^i, y_j^i)} is given like that of
  \code{\link{SimulateIP}}.

}
\value{
  \item{parents}{a list containing two components named "\code{n}" and
    "\code{xy}" giving the number and the matrix of \code{(x,y)} coordinates of
    simulated parents points respectively.}
  \item{offspring}{a list containing two components named "\code{n}" and
    "\code{xy}" giving the number and the matrix of \code{(x,y)} coordinates of
    simulated offspring points respectively.}
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models,
  \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
pars <- c(mu = 50.0, nu = 30.0, sigma = 0.03)
SimulateThomas(pars, seed = 117)
}

\keyword{spatial}
