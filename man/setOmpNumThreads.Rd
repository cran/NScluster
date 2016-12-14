\name{setOmpNumThreads}
\alias{setOmpNumThreads}
\title{Set the Number of OpenMP Threads}
\description{
  Set the number of threads used for parallel computation based on OpenMP.
}
\usage{
setOmpNumThreads(nthreads=NULL)
}

\arguments{
  \item{nthreads}{the number of OpenMP threads. If \code{nthreads} is NULL or greater than the number of processors, it is set to the number as omp_get_num_procs().}
}
\value{
  \code{setOmpNumThreads} returns the maximum number of threads available.
}
\keyword{ts}
