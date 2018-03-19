\name{NScluster-package}
\alias{NScluster-package}
\alias{NScluster}
\docType{package}

\title{Simulation and Estimation of the Neyman-Scott Type Spatial Cluster Models}

\description{
  This package provides functions for simulation and estimation of spatial cluster point pattern of Neyman-Scott models and their extensions. 
We adopt the Simplex estimation to maximize the Palm likelihood function.
}

\details{
  The documentation 'NScluster: An R Package for Simulation and Estimation of the Neyman-Scott type spatial cluster models' is available in \url{../doc/NScluster-guide.pdf}.

  Simulation:

 \code{\link{SimulateThomas}}, \code{\link{SimulateIP}}, \code{\link{SimulateTypeA}}, \code{\link{SimulateTypeB}}
 and \code{\link{SimulateTypeC}} simulate spatial cluster point pattern of Neyman-Scott models
 and their extensions. We describe overview of those models briefly in the NScluster documentation \url{../doc/NScluster-guide.pdf}.

 Simulation method of each model is described under the corresponding topic.

  Parameter estimation:

 We adopt the Simplex estimation to maximize the Palm likelihood function (or minimize the negative Palm likelihood function).
 The maximum Palm likelihood estimate is called 'MPLE', for short.
 The Palm intensity function and the analytical form of the Palm log-likelihood of the Tomas model,
 Type B model and Type C model are described under the topic \code{\link{EstimateThomas}}, \code{\link{EstimateTypeB}} and \code{\link{EstimateTypeC}}, respectively.
 On the other hand, for \code{\link{EstimateIP}} and \code{\link{EstimateTypeA}},
 we need to take the alternative form without explicit representation of the Palm intensity function, which need very long c.p.u. time in the minimization procedure.
 We parallelize the minimization procedure with OpenMP.

 \code{\link{PalmThomas}}, \code{\link{PalmIP}}, \code{\link{PalmTypeA}}, \code{\link{PalmTypeB}} and \code{\link{PalmTypeC}}
 calculate the non-parametric Palm intensity function estimated directory from a set of point pattern data.
}

\references{
U. Tanaka, Y. Ogata and K. Katsura, Simulation and estimation of the Neyman-Scott type spatial cluster models, \emph{Computer Science Monographs} \bold{No.34}, 2008, 1-44, 
The Institute of Statistical Mathematics, Tokyo. \url{http://www.ism.ac.jp/editsec/csm/}

U.Tanaka, Y. Ogata and D. Stoyan, Parameter estimation and model selection for Neyman-Scott point processes, \emph{Biometrical Journal}, \bold{50}, 2008, 43-57.
}

\keyword{ package }

