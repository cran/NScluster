\name{mple.cppm}
\alias{mple.cppm}
\alias{coef.mple}
\alias{summary.mple}
\title{MPLE of Neyman-Scott Cluster Point Process Models and Their Extensions}
\description{
MPLE of the five cluster point process models.
}
\usage{
mple.cppm(model = "Thomas", xy.points, pars = NULL, eps = 0.001, uplimit = 0.3,
          skip = 1)

\method{coef}{mple}(object, \dots)
\method{summary}{mple}(object, \dots)
}
\arguments{
  \item{model}{a character string indicating each cluster point process model:
    \code{"Thomas"}, \code{"IP"}, \code{"TypeA"}, \code{"TypeB"}, and
    \code{"TypeC"}.}
  \item{xy.points}{a matrix containing the coordinates \code{(x,y)} of points
    in \eqn{W=[0,1]\times[0,1]}{W=[0,1]*[0,1]}.}
  \item{pars}{a named vector containing a given initial guess of each
    parameter. If \code{NULL}, any suitable parameters are used. See Details' in
    \code{\link{sim.cppm}} for the parameters of each model.}
  \item{eps}{the sufficiently small number to implement the optimization
    procedure for the log-Palm likelihood function. The procedure is iterated
    at most 1000 times until the \code{process2$stderr} becomes smaller than the
    \code{eps}.}
  \item{uplimit}{upper limit in place of \eqn{\infty} of the integral in the
    probability distribution function relative to the random distance between
    two descendant points within the same cluster. The \code{uplimit} is valid
    for \code{"IP"} and \code{"TypeA"}.}
  \item{skip}{the variable enables one to obtain speedily the initial MPLEs, but
    rough approximation. The \code{skip} calculates the Palm intensity function 
    of the log-Palm likelihood function for every \code{skip}-th
    \eqn{r_{ij}}{r(i,j)} in the ordered distances of the pairs \eqn{i} and 
    \eqn{j}. The \code{skip} is valid for \code{"IP"} and \code{"TypeA"}.}
  \item{object}{an object of class \code{"mple"}.}
  \item{\dots}{ignored.}
}

\details{
  \itemize{
    \item{\code{"Thomas"}}{
      (Thomas model)

      The Palm intensity function is given as follows:

      For all \eqn{r \ge 0}{r >= 0},

      \deqn{\lambda_{\bm{o}}(r) = \mu\nu + \frac{\nu}{4\pi \sigma^2}
        \exp \left( -\frac{r^2}{4 \sigma^2} \right).}{%
        \lambda_o(r) = \mu\nu + \nu/(4\pi\sigma^2) * exp(-r^2/(4\sigma^2)).}

      The log-Palm likelihood function is given by

      \deqn{\log L(\mu,\nu,\sigma) = \sum_{\{i,j; i<j, r_{ij} \le 1/2\}}
        \log \nu \left\{ \mu + \frac{1}{4 \pi \sigma^2}
        \exp \left( -\frac{{r_{ij}}^2}{4 \sigma^2} \right) \right\}}{%
        log L(\mu,\nu,\sigma) = \sum_{i, j; i < j, r(i, j) <= 1/2}
        log\nu { \mu + 1/(4\pi\sigma^2) * exp(-r(i, j)^2/(4\sigma^2)) }}

      \deqn{- N(W)\nu \left\{ \frac{\pi \mu}{4} + 1 -
        \exp \left( -\frac{1}{16 \sigma^2} \right) \right\}.}{%
        - N(W) \nu { \pi\mu/4 + 1 - exp(-1/(16\sigma^2)) }.}
	}
  }
  \itemize{
    \item{\code{"TypeB"}}{
      (Type B model)

      The Palm intensity function is given as follows:

      For all \eqn{r \ge 0}{r >= 0},

      \deqn{\lambda_{\bm{o}}(r) = \lambda + \frac{\nu}{4 \pi} 
        \left\{ \frac{a}{{\sigma_1}^2} \exp \left( -\frac{r^2}{4{\sigma_1}^2}
        \right)+ \frac{(1-a)}{{\sigma_2}^2} \exp
		\left( -\frac{r^2}{4{\sigma_2}^2} \right) \right\},}{%
        \lambda_o(r) = \lambda + \nu/(4\pi) * 
		{ a/\sigma1^2 * exp(-r^2/(4\sigma1^2)) +
		(1-a)/\sigma2^2 * exp(-r^2/(4\sigma2^2)) },}

      where \eqn{\lambda = \nu(\mu_1+\mu_2)}{\lambda = \nu(\mu1+\mu2)} and
      \eqn{a = \mu_1/(\mu_1+\mu_2)}{a = \mu1/(\mu1+\mu2)} are the total
      intensity and the ratio of the intensity of the parent points of the
      smaller cluster to the total one, respectively.

      The log-Palm likelihood function is given by

      \eqn{\log L(\lambda, \alpha, \beta, \sigma_1, \sigma_2)}{log L(\lambda,
      \alpha, \beta, \sigma1, \sigma2)} 

      \deqn{=\sum_{\{i,j; i<j, r_{ij} \le 1/2\}} \log
        \left[ \lambda + \frac{1}{4 \pi} \left\{ \frac{\alpha}{{\sigma_1}^2}
		\exp \left( -\frac{{r_{ij}}^2}{4{\sigma_1}^2} \right) + 
        \frac{\beta}{{\sigma_2}^2} \exp \left( -\frac{{r_{ij}}^2}{4{\sigma_2}^2}
	    \right) \right\} \right]}{%
        =\sum_{i, j; i < j, r(i, j) <= 1/2} log[ \lambda + 1/(4\pi)
	    { \alpha/\sigma1^2 * exp(-r(i, j)^2/(4\sigma1^2)) + \beta/\sigma2^2 * 
	    exp(-r(i, j)^2/(4\sigma2^2))}]}
      \deqn{- N(W) \left[ \frac{\pi \lambda}{4} + \alpha \left\{ 1 - \exp
        \left( -\frac{1}{16{\sigma_1}^2} \right) \right\} + \beta
		\left\{ 1- \exp \left( -\frac{1}{16{\sigma_2}^2} \right) \right\}
		\right],}{%
        - N(W) [ \pi\lambda/4 + \alpha{ 1-exp(-1/(16\sigma1^2)) } +
        \beta{ 1-exp(-1/(16\sigma2^2)) } ],}

      where \eqn{\alpha = a\nu} and \eqn{\beta = (1-a)\nu}.
    }
  }
  \itemize{
    \item{\code{"TypeC"}}{
      (Type C model)

      The Palm intensity function is  given as follows:

      For all \eqn{r \ge 0}{r >= 0},

      \deqn{\lambda_{\bm{o}}(r) = \lambda + \frac{1}{4 \pi} \left\{
        \frac{a\nu_1}{{\sigma_1}^2} \exp
		\left( -\frac{r^2}{4{\sigma_1}^2} \right) +
        \frac{(1-a)\nu_2}{{\sigma_2}^2} \exp \left( -\frac{r^2}{4{\sigma_2}^2}
        \right) \right\},}{%
        \lambda_o(r) = \lambda + 1/(4\pi) * { a\nu1/\sigma1^2 *
        exp(-r^2/(4\sigma1^2)) + (1-a)\nu2/\sigma2^2 * exp(-r^2/(4\sigma2^2))},}

      where \eqn{\lambda = \mu_1\nu_1 + \mu_2\nu_2}{%
	  \lambda = \mu1\nu1 + \mu2\nu2} and \eqn{a = \mu_1\nu_1/\lambda}{%
      a = \mu1\nu1/\lambda} are the total intensity and the
      ratio of the intensity of the smaller cluster to the total one,
      respectively.

      The log-Palm likelihood function is given by

      \eqn{\log L(\lambda, \alpha, \beta, \sigma_1, \sigma_2)}{%
         log L(\lambda, \alpha, \beta, \sigma1, \sigma2)}
      \deqn{= \sum_{\{i,j; i<j, r_{ij} \le 1/2\}} \log \left[ \lambda + \frac{1}
        {4 \pi} \left\{ \frac{\alpha}{{\sigma_1}^2} \exp
        \left( -\frac{{r_{ij}}^2}{4{\sigma_1}^2} \right) +
        \frac{\beta}{{\sigma_2}^2} \exp
		\left( -\frac{{r_{ij}}^2}{4{\sigma_2}^2} \right) \right\} \right]}{% 
        =\sum_{i, j; i < j, r(i, j) <= 1/2} log[ \lambda + 1/(4\pi) *
        { \alpha/\sigma1^2 * exp(-r(i, j)^2/(4\sigma1^2)) + \beta/\sigma2^2 *
        exp(-r(i, j)^2/(4\sigma2^2)) } ]}

      \deqn{-N(W) \left[ \frac{\pi\lambda}{4} + \alpha \left\{ 1 - \exp
        \left( -\frac{1}{16{\sigma_1}^2} \right) \right\} + \beta \left\{ 1-
        \exp \left( -\frac{1}{16{\sigma_2}^2} \right) \right\} \right],}{%
        -N(W) [ \pi\lambda/4 + \alpha{ 1-exp(-1/(16\sigma1^2)) } +
        \beta{ 1-exp(-1/(16\sigma2^2)) } ],}

      where \eqn{\alpha = a\nu_1}{\alpha = a\nu1} and
      \eqn{\beta = (1-a)\nu_2}{\beta = (1-a)\nu2}.
    }
  }
  For the inverse-power model and the Type A models, we need to take the 
  alternative form without explicit representation of the Palm intensity 
  function. See the second reference below for details.
}
\value{
  \code{mple.cppm} returns an object of class \code{"mple"} containing the
    following main components:
  \item{mple}{MPLE (maximum Palm likelihood estimate).}
  \item{log.mpl}{the log maximum Palm likelihood.}
  \item{aic}{AIC.}
  \item{process1}{a list with following components.
    \describe{
      \item{cflg}{1 (="update") or -1 (="testfn"), where "update" indicates that
        -log L value has attained the minimum so far, otherwise not.} 
      \item{logl}{the minimized -log L in the process to minimize the
        negative log-Palm likelihood function.}
      \item{mples}{corresponding MPLEs.}
    }
  }
  \item{process2}{a list with following components.
    \describe{
      \item{logl.simplex}{the minimized -log L by the simplex method.}
      \item{stderr}{standard error.}
      \item{pa.normal}{the normalized variables corresponding to the MPLEs
        relative to the initial estimates.}
    }
  }

  There are other methods \code{\link{plot.mple}} and \code{\link{print.mple}}
  for this class.
}
\references{
Tanaka, U., Ogata, Y. and Katsura, K. (2008)
 Simulation and estimation of the Neyman-Scott type spatial cluster models.
 \emph{Computer Science Monographs} \bold{34}, 1-44.
 The Institute of Statistical Mathematics, Tokyo.
 \url{https://www.ism.ac.jp/editsec/csm/}.

Tanaka, U., Ogata, Y. and Stoyan, D. (2008)
 Parameter estimation and model selection for Neyman-Scott point processes.
 \emph{Biometrical Journal} \bold{50}, 43-57.
}
\examples{
\dontrun{
# The computation of MPLEs takes a long CPU time in the minimization procedure,
# especially for the Inverse-power type and the Type A models.

### Thomas Model
 # simulation
 pars <- c(mu = 50.0, nu = 30.0, sigma = 0.03)
 t.sim <- sim.cppm("Thomas", pars, seed = 117)
 ## estimation
 init.pars <- c(mu = 40.0, nu = 40.0, sigma = 0.05)
 t.mple <- mple.cppm("Thomas", t.sim$offspring$xy, init.pars)
 coef(t.mple)

### Inverse-Power Type Model
 # simulation
 pars <- c(mu = 50.0, nu = 30.0, p = 1.5, c = 0.005)
 ip.sim <- sim.cppm("IP", pars, seed = 353)
 ## estimation
 init.pars <- c(mu = 55.0, nu = 35.0, p = 1.0, c = 0.01)
 ip.mple <- mple.cppm("IP", ip.sim$offspring$xy, init.pars, skip = 100)
 coef(ip.mple)

### Type A Model
 # simulation
 pars <- c(mu = 50.0, nu = 30.0, a = 0.3, sigma1 = 0.005, sigma2 = 0.1)
 a.sim <- sim.cppm("TypeA", pars, seed = 575)
 ## estimation
 init.pars <- c(mu = 60.0, nu = 40.0, a = 0.5, sigma1 = 0.01, sigma2 = 0.1)
 a.mple <- mple.cppm("TypeA", a.sim$offspring$xy, init.pars, skip = 100)
 coef(a.mple)

### Type B Model
 # simulation
 pars <- c(mu1 = 10.0, mu2 = 40.0, nu = 30.0, sigma1 = 0.01, sigma2 = 0.03)
 b.sim <- sim.cppm("TypeB", pars, seed = 257)
 ## estimation
 init.pars <- c(mu1 = 20.0, mu2 = 30.0, nu = 30.0, sigma1 = 0.02, sigma2 = 0.02)
 b.mple <- mple.cppm("TypeB", b.sim$offspring$xy, init.pars)
 coef(b.mple)

### Type C Model
 # simulation
 pars <- c(mu1 = 5.0, mu2 = 9.0, nu1 = 30.0, nu2 = 150.0,
           sigma1 = 0.01, sigma2 = 0.05)
 c.sim <- sim.cppm("TypeC", pars, seed = 555)
 ## estimation
 init.pars <- c(mu1 = 10.0, mu2 = 10.0, nu1 = 30.0, nu2 = 120.0,
                sigma1 = 0.03, sigma2 = 0.03)
 c.mple <- mple.cppm("TypeC", c.sim$offspring$xy, init.pars)
 coef(c.mple)
}
}

\keyword{spatial}
