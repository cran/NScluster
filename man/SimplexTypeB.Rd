\name{SimplexTypeB}
\alias{SimplexTypeB}
\title{Parameter estimation of Type B model}
\description{
  Parameter estimation of Type B Model by using The Palm Log-Likelihood Function.
}
\usage{
  SimplexTypeB(offspring, pa, ty=1, eps=0.1e-2, process=0, plot=TRUE)
}
\arguments{
  \item{offspring}{the list of coordinates \eqn{(x,y)} of offspring points.}
  \item{pa}{the initial guesses of the parameters \eqn{(\mu_1, \mu_2, \nu, \sigma_1, \sigma_2)}{(mu_1, mu_2, nu, sigma_1, sigma_2)}.}
  \item{ty}{the variable \eqn{T_y}{Ty} for the standardized coordinates of points in the rectangular region \eqn{[0,1] \times [0,T_y]}{[0,1]*[0,Ty]}
            without loss of generality except for the scaling.}
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
  \item{logl.palm}{the minimized -log L in the process of minimizing the negative Palm log-likelihood function.}
  \item{mple}{the MPLE (maximum Palm likelihood estimator) values corresponding to the above \code{logl.palm}.}
  \item{logl.simplex}{the minimized -log L by the simplex method.}
  \item{stderr}{the standard deviations.}
  \item{pa.normal}{the normalized variables corresponding to the MPLE values relative to the initial estimates.}
}
\details{
  The Palm intensity function of the Type B model is calculated as follows:
  For any \eqn{r \ge 0}{r >= 0},

  \deqn{\lambda_0(r) = \lambda + \frac{\nu}{4 \pi} \left\{ \frac{a}{\sigma_1^2} \exp \left( -\frac{r^2}{4 \sigma_1^2} \right)+ \frac{(1-a)}{\sigma_2^2} \exp \left( -\frac{r^2}{4 \sigma_2^2} \right) \right\},}{%
         lambda_0(r) = lambda + nu/(4 \pi){a/sigma_1^2 exp(-r^2/(4sigma_1^2))+ (1-a)/sigma_2^2 exp(-r^2/(4sigma_2^2))},}

  where \eqn{\lambda = \nu(\mu_1+\mu_2)}{lambda = nu(mu_1+mu_2)} is the total population size and \eqn{a = \mu_1/(\mu_1+\mu_2)}{a = mu_1/(mu_1+mu_2)} is
the ratio of the parent points of the smaller sized cluster to the total ones. 

  The Palm log-likelihood function of the Type B model is analytically calculated as follows: \cr

  \eqn{\log L(\lambda, \alpha, \beta, \sigma_1, \sigma_2)}{log L(lambda, alpha, beta, sigma_1, sigma_2)}
  \deqn{=\sum_{\{i,j; i \ne j, r_{ij}<R\}} \log \left[ \lambda + \frac{1}{4 \pi} \left\{ \frac{\alpha}{\sigma_1^2} \exp \left( -\frac{r_{ij}^2}{4 \sigma_1^2} \right) + \frac{\beta}{\sigma_2^2} \exp \left( -\frac{r_{ij}^2}{4 \sigma_2^2} \right) \right\} \right]}{%
       = sum_{i,j; i<>j, r_{ij}<R} log[lambda + 1/(4 \pi){alpha/sigma_1^2 exp(-r_{ij}^2/(4sigma_1^2))+ beta/sigma_2^2 exp(-r_{ij}^2/(4sigma_2^2))}]}
  \deqn{- N(W) \left[ \pi R^2\lambda + \alpha \left\{ 1 - \exp \left( -\frac{R^2}{4 \sigma_1^2} \right) \right\} + \beta \left\{ 1- \exp \left( -\frac{R^2}{4 \sigma_2^2} \right) \right\} \right],}{%
       - N(W)[\pi R^2 lambda + alpha{1-exp(-R^2/(4sigma_1^2))} + beta{1-exp(-R^2/(4sigma_2^2))}],}

  with \eqn{R = 1/2}, where \eqn{\alpha = a\nu}{alpha = a*nu} and \eqn{\beta = (1-a)\nu}{beta = (1-a)nu}. 
}
\references{
  U. Tanaka, Y. Ogata and K. Katsura,
  Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44.
  The Institute of Statistical Mathematics.
}
\examples{
  ## simulation
  seeds <- c(822, 913, 905)
  mu1 <- 10;  mu2 <- 40; nu <- 30;  sig1 <- 0.01;  sig2 <- 0.03
  z <- SimulateTypeB(seeds, c(mu1,mu2,nu,sig1,sig2), pmax=100, omax=150)

  \dontrun{
    ## estimation
    ## need very long c.p.u time in the minimization procedure
    SimplexTypeB(z$offspring, c(mu1,mu2,nu,sig1,sig2), process=3)
  }
}

\keyword{spatial}
