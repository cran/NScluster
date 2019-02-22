\name{print.mple}
\alias{print.mple}
\title{Print Process for Maximizing Log-Palm Likelihood Function}
\description{
  Print the process for minimizing the negative log-Palm likelihood function 
  and/or the process for optimizing the normalized parameters depending on a
  given initial guess of each parameter by the simplex method.
}
\usage{
\method{print}{mple}(x, print.level = 0, \dots)
}
\arguments{
 \item{x}{an object of class \code{"mple"} returned by \code{\link{mple.cppm}}.}
 \item{print.level}{We have the following processes:
    \describe{
      \item{0}{output initial values and MPLE.}
      \item{1}{output the process for minimizing the negative log-Palm 
               likelihood function, in addition. (\code{x$process1})}
      \item{2}{output the process for optimizing the normalized parameters
               depending on a given initial guess of each parameter by the 
               simplex method, in addition. (\code{x$process2})}
      \item{3}{output both processes.}
    }}
 \item{\dots}{ignored.}
}
\keyword{spatial}
