\name{SimulateIP}
\alias{SimulateIP}
\title{Simulation of the Inverse-Power Type Model}
\description{
  Simulation of the Inverse-Power Type Model.
}
\usage{
  SimulateIP(seeds, pa, ty=1, pmax=100, omax=3000, plot=TRUE)
}
\arguments{
  \item{seeds}{the three positive integer variables, which are initial seeds for a sequence of uniform random numbers.}
  \item{pa}{the parameter values \eqn{(\mu, \nu, p, c)}{(mu, nu, p, c)}.}
  \item{ty}{the variable \eqn{T_y}{Ty} for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]} without loss of generality except for the scaling.}
  \item{pmax}{maximum number of parent points.}
  \item{omax}{maximum number of offspring points.}
  \item{plot}{logical. If TRUE (default) simulated parent points and offspring points are plotted.}
}
\value{
  \item{n.parents}{the number of simulated parent points.}
  \item{parents}{the coordinates of simulated parent points.}
  \item{n.offspring}{the total number of simulated offspring points.}
  \item{offspring}{the coordinates of simulated offspring points.}
}
\details{
  Let random variable \eqn{U} be independently and uniformly distributed in [0,1]. For any \eqn{r \ge 0}{r>=0},

  \deqn{Q_{p,c}(r) := \int_0^r q_{p,c}(t)dt}{%
        Q_{p,c}(r) := int_0^r q_{p,c}(t)dt}
  \deqn{= c^{p-1}(p-1) \frac{(r+c)^{1-p} - c^{1-p}}{1-p}}{%
        = c^{p-1}(p-1) {(r+c)^{1-p} - c^{1-p}}/(1-p)}
  \deqn{= 1 - c^{p-1} (r+c)^{1-p}.}{%
        = 1 - c^{p-1}*(r+c)^{1-p}.}

  Here, we put \eqn{Q_{p,c}(r) = U}. From this, we have

  \deqn{r = c\{(1-U)^{1/(1-p)} - 1\}.}{%
        r = c{(1-U)^{1/(1-p)} - 1}.}

  Similarly, coordinate of the offspring points \eqn{(x_j^i, y_j^i), j=1,2,...,Poisson(\nu)}{(x_j^i, y_j^i), j=1,2,...,Poisson(nu)} with
  the inverse-power type is given for each \eqn{i=1,2,...,Poisson(\mu)}{i=1,2,...,Poisson(mu)},

  \deqn{x_j^i = x_i^p + r \cos(2 \pi U),}{%
        x_j^i = x_i^p + r cos(2 \pi U),}
  \deqn{y_j^i = y_i^p + r \sin(2 \pi U),}{%
        y_j^i = y_i^p + r sin(2 \pi U),}

  using series of different uniform random number \eqn{\{U\}}{{U}} for different \eqn{i} and \eqn{j}.
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  seeds <- c(822, 913, 905)
  mu <- 50;  nu <- 30;  p <- 1.5;  c <- 0.005
  SimulateIP(seeds, c(mu,nu,p,c))
}

\keyword{spatial}
