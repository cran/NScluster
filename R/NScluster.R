SimulateIP <- function(pars, seed = NULL, plot = TRUE) {

# (mu, nu) for the random variable Poisson
# p : a decay order with respect to the distance
# c : scaling factor with respect to the distance
  pname <- c("mu", "nu", "p", "c")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  p <- pa[3]
  c <- pa[4]
# seed : initial seeds for a sequence of uniform random numbers
  ix <- seed
  if (is.null(ix))
    ix <- -1
# ty : the variable for the standardized coordinates of points
#  in the rectangular region
  ty <- 1
# maximum number of parent points / maximum number of offspring points
  pmax <- 100
  omax <- 100


  z <- .Call("simIP",
             as.integer(ix),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(p),
             as.double(c),
             as.integer(pmax),
             as.integer(omax))

  ier <- z[[7L]]
  if (ier == -1) {
    stop(paste("too many parents (the default maximum number is ", pmax, ")"),
         call. = FALSE)
  } else if (ier == -2) {
    stop(paste("too many offspring (the default maximum number is ", omax, ")"),
         call. = FALSE)
  } else {

    npts <- z[[1L]]
    ncl <- z[[2L]]
    parents.x <- z[[3L]][1:npts]
    parents.y <- z[[4L]][1:npts]
    xcl <- array(z[[5L]], dim = c(pmax, omax))
    ycl <- array(z[[6L]], dim = c(pmax, omax))

    parents.xy <- array(c(parents.x, parents.y), dim = c(npts, 2))

    offspring.x <- NULL
    offspring.y <- NULL
    for (i in 1:npts) {
      offspring.x <- c(offspring.x, xcl[i, 1:ncl[i]])
      offspring.y <- c(offspring.y, ycl[i, 1:ncl[i]])
    }
    np <- sum(ncl[1:npts])
    offspring.xy <- array(c(offspring.x, offspring.y), dim = c(np, 2))

    if (plot == TRUE)
      plot.Simulate("IP", parents.xy, offspring.xy, c(mu, nu, p, c))

    SimulateIP.out <- list(parents = list(n = npts, xy = parents.xy),
                           offspring = list(n = np, xy = offspring.xy))

    class(SimulateIP.out) <- "Simulate"
    return(SimulateIP.out)
  }
}


SimulateThomas <- function(pars, seed = NULL, plot = TRUE)  {

# the parameter values
# (mu, nu, sigma) for the random variable Poisson
  pname <- c("mu", "nu", "sigma")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  sigma <- pa[3]
# seeds : initial seeds for a sequence of uniform random numbers
  ix <- seed
  if (is.null(ix))
    ix <- -1
# ty : the variable for the standardized coordinates of points
#  in the rectangular region
  ty <- 1
# maximum number of parent points / maximum number of offspring points
  pmax <- 100
  omax <- 100

  z <- .Call("simThom",
             as.integer(ix),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(sigma),
             as.integer(pmax),
             as.integer(omax))

  ier <- z[[7L]]
  if (ier == -1) {
    stop(paste("too many parents (the default maximum number is ", pmax, ")"),
         call. = FALSE)
  } else if (ier == -2) {
    stop(paste("too many offspring (the default maximum number is ", omax, ")"),
         call. = FALSE)
  } else {

    npts <- z[[1L]]
    ncl <- z[[2L]]
    parents.x <- z[[3L]][1:npts]
    parents.y <- z[[4L]][1:npts]
    xcl <- array(z[[5L]], dim = c(pmax, omax))
    ycl <- array(z[[6L]], dim = c(pmax, omax))

    parents.xy <- array(c(parents.x, parents.y), dim = c(npts, 2))

    offspring.x <- NULL
    offspring.y <- NULL
    for (i in 1:npts) {
      offspring.x <- c(offspring.x, xcl[i, 1:ncl[i]])
      offspring.y <- c(offspring.y, ycl[i, 1:ncl[i]])
    }
    np <- sum(ncl[1:npts])
    offspring.xy <- array(c(offspring.x, offspring.y), dim = c(np, 2))

    if (plot == TRUE)
      plot.Simulate("T", parents.xy, offspring.xy, c(mu, nu, sigma))

    SimulateThomas.out <- list(parents = list(n = npts, xy = parents.xy),
                               offspring = list(n = np, xy = offspring.xy))

    class(SimulateThomas.out) <- "Simulate"
    return(SimulateThomas.out)
  }
}


SimulateTypeA <- function(pars, seed = NULL, plot = TRUE) {

# the parameter values
# (mu, nu, a, sigma1, sigma2) for the random variable Poisson
  pname <- c("mu", "nu", "a", "sigma1", "sigma2")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  a <- pa[3]
  sigma1 <- pa[4]
  sigma2 <- pa[5]
# seeds : initial seeds for a sequence of uniform random numbers
  ix <- seed
  if (is.null(ix))
    ix <- -1
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# maximum number of parent points / maximum number of offspring points
  pmax <- 100
  omax <- 100

  z <- .Call("simA",
             as.integer(ix),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(a),
             as.double(sigma1),
             as.double(sigma2),
             as.integer(pmax),
             as.integer(omax))

  ier <- z[[7L]]
  if (ier == -1) {
    stop(paste("too many parents (the default maximum number is ", pmax, ")"),
         call. = FALSE)
  } else if (ier == -2) {
    stop(paste("too many offspring (the default maximum number is ", omax, ")"),
         call. = FALSE)
  } else {

    npts <- z[[1L]]
    ncl <- z[[2L]]
    parents.x <- z[[3L]][1:npts]
    parents.y <- z[[4L]][1:npts]
    xcl <- array(z[[5L]], dim = c(pmax, omax))
    ycl <- array(z[[6L]], dim = c(pmax, omax))

    parents.xy <- array(c(parents.x, parents.y), dim = c(npts, 2))

    offspring.x <- NULL
    offspring.y <- NULL
    for (i in 1:npts) {
      offspring.x <- c(offspring.x, xcl[i, 1:ncl[i]])
      offspring.y <- c(offspring.y, ycl[i, 1:ncl[i]])
    }
    np <- sum(ncl[1:npts])
    offspring.xy <- array(c(offspring.x, offspring.y), dim = c(np, 2))

    if (plot == TRUE)
      plot.Simulate("A", parents.xy, offspring.xy, c(mu, nu, a, sigma1, sigma2))

    SimulateTypeA.out <- list(parents = list(n = npts, xy = parents.xy),
                              offspring = list(n = np, xy = offspring.xy))

    class(SimulateTypeA.out) <- "Simulate"
    return(SimulateTypeA.out)
  }
}


SimulateTypeB <- function(pars, seed = NULL, parents.distinct = FALSE,
                 plot = TRUE) {

# the parameter values
# (mu1, mu2, nu, sigma1, sigma2)  for the random variable Poisson
  pname <- c("mu1", "mu2", "nu", "sigma1", "sigma2")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu1 <- pa[1]
  mu2 <- pa[2]
  nu <- pa[3]
  sigma1 <- pa[4]
  sigma2 <- pa[5]
# seeds : initial seeds for a sequence of uniform random numbers
  ix <- seed
  if (is.null(ix))
    ix <- -1
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# maximum number of parent points / maximum number of offspring points
  pmax <- 100
  omax <- 100

  z <- .Call("simB",
             as.integer(ix),
             as.double(ty),
             as.double(mu1),
             as.double(mu2),
             as.double(nu),
             as.double(sigma1),
             as.double(sigma2),
             as.integer(pmax),
             as.integer(omax))

  if (z[[13L]] == -1 || z[[13L]] == -2) {
    stop(paste("too many parents (the default maximum number is ", pmax, ")"),
         call. = FALSE)
  } else if (z[[13L]] == -11 || z[[13L]] == -22) {
    stop(paste("too many offspring (the default maximum number is ", omax, ")"),
         call. = FALSE)
  } else {

    m1 <- z[[1L]]
    m2 <- z[[7L]]
    parents1.x <- z[[3L]][1:m1]
    parents1.y <- z[[4L]][1:m1]
    parents2.x <- z[[9L]][1:m2]
    parents2.y <- z[[10L]][1:m2]
    parents.xy <- array(0, dim = c(m1 + m2, 2))
    parents.xy[, 1] <- c(parents1.x, parents2.x)
    parents.xy[, 2] <- c(parents1.y, parents2.y)

    np1 <- sum(z[[2L]][1:m1])
    offspring1.x <- z[[5L]][1:np1]
    offspring1.y <- z[[6L]][1:np1]
    np2 <- sum(z[[8L]][1:m2])
    offspring2.x <- z[[11L]][1:np2]
    offspring2.y <- z[[12L]][1:np2]
    offspring.xy <- array(0, dim = c(np1 + np2, 2))
    offspring.xy[, 1] <- c(offspring1.x, offspring2.x)
    offspring.xy[, 2] <- c(offspring1.y, offspring2.y)

    if (plot == TRUE) {
      if (parents.distinct == FALSE) {
        plot.Simulate("B", parents.xy, offspring.xy,
                      c(mu1, mu2, nu, sigma1, sigma2))
      } else {
        parents1.xy <- array(c(parents1.x, parents1.y), dim = c(m1, 2))
        parents2.xy <- array(c(parents2.x, parents2.y), dim = c(m2, 2))
        offspring1.xy <- array(c(offspring1.x, offspring1.y), dim = c(np1, 2))
        offspring2.xy <- array(c(offspring2.x, offspring2.y), dim = c(np2, 2))
        plot.Simulate2("B", parents=list(xy1 = parents1.xy, xy2 = parents2.xy),
                       offspring=list(xy1 = offspring1.xy, xy2 = offspring2.xy),
                       c(mu1, mu2, nu, sigma1, sigma2))
      }
    }

    SimulateTypeB.out <- list(parents = list(n = c(m1, m2), xy = parents.xy),
                         offspring = list(n = c(np1, np2), xy = offspring.xy))
    class(SimulateTypeB.out) <- "Simulate"
    return(SimulateTypeB.out)
  }
}


SimulateTypeC <- function(pars, seed = NULL, parents.distinct = FALSE,
                 plot = TRUE) {

# the parameter values
# (mu1, mu2, nu1, nu2, sigma1, sigma2) for the random variable Poisson
  pname <- c("mu1", "mu2", "nu1", "nu2", "sigma1", "sigma2")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu1 <- pa[1]
  mu2 <- pa[2]
  nu1 <- pa[3]
  nu2 <- pa[4]
  sigma1 <- pa[5]
  sigma2 <- pa[6]
# seeds : initial seeds for a sequence of uniform random numbers
  ix <- seed
  if (is.null(ix))
    ix <- -1
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# maximum number of parent points / maximum number of offspring points
  pmax <- 200
  omax <- 200

  z <- .Call("simC",
             as.integer(ix),
             as.double(ty),
             as.double(mu1),
             as.double(mu2),
             as.double(nu1),
             as.double(nu2),
             as.double(sigma1),
             as.double(sigma2),
             as.integer(pmax),
             as.integer(omax))

  if (z[[13L]] == -1 || z[[13L]] == -2) {
    stop(paste("too many parents (the default maximum number is ", pmax, ")"),
         call. = FALSE)
  } else if (z[[13L]] == -11 || z[[13L]] == -22) {
    stop(paste("too many offspring (the default maximum number is ", omax, ")"),
         call. = FALSE)
  } else {

    m1 <- z[[1L]]
    m2 <- z[[7L]]
    parents1.x <- z[[3L]][1:m1]
    parents1.y <- z[[4L]][1:m1]
    parents2.x <- z[[9L]][1:m2]
    parents2.y <- z[[10L]][1:m2]
    parents.xy <- array(0, dim = c( m1 + m2, 2 ))
    parents.xy[, 1] <- c(parents1.x, parents2.x)
    parents.xy[, 2] <- c(parents1.y, parents2.y)

    np1 <- sum(z[[2L]][1:m1])
    np2 <- sum(z[[8L]][1:m2])
    offspring1.x <- z[[5L]][1:np1]
    offspring1.y <- z[[6L]][1:np1]
    offspring2.x <- z[[11L]][1:np2]
    offspring2.y <- z[[12L]][1:np2]
    offspring.xy <- array(0, dim = c(np1 + np2, 2))
    offspring.xy[, 1] <- c(offspring1.x, offspring2.x)
    offspring.xy[, 2] <- c(offspring1.y, offspring2.y)

    if (plot == TRUE) {
      if (parents.distinct == FALSE) {
        plot.Simulate("C",  parents.xy, offspring.xy,
                      c(mu1, mu2, nu1, nu2, sigma1, sigma2))
      } else {
        parents1.xy <- array(c(parents1.x, parents1.y), dim = c(m1, 2))
        parents2.xy <- array(c(parents2.x, parents2.y), dim = c(m2, 2))
        offspring1.xy <- array(c(offspring1.x, offspring1.y), dim = c(np1, 2))
        offspring2.xy <- array(c(offspring2.x, offspring2.y), dim = c(np2, 2))
        plot.Simulate2("C", parents=list(xy1 = parents1.xy, xy2 = parents2.xy),
                       offspring=list(xy1 = offspring1.xy, xy2 = offspring2.xy),
                       c(mu1, mu2, nu1, nu2, sigma1, sigma2))
      }
    }

    SimulateTypeC.out <- list(parents = list(n = c(m1,m2), xy = parents.xy),
                              offspring = list(n = c(np1,np2), xy = offspring.xy))
    class(SimulateTypeC.out) <- "Simulate"
    return(SimulateTypeC.out)
  }
}


EstimateIP <- function(xy.points, pars, eps = 0.001, uplimit = 0.3, skip = 1,
              process.report = 0, plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)
# pa : the parameter values
# (mu, nu) for the random variable Poisson
# p : a decay order with respect to the distance
# c : scaling factor with respect to the distance
  pname <- c("mu", "nu", "p", "c")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  p <- pa[3]
  c <- pa[4]
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# parameter
  ipflg <- process.report
  if (process.report == 0 && plot == TRUE)
    ipflg <- 2
  if (process.report == 1 && plot == TRUE)
    ipflg <- 3

  n <- 4
  itmax <- 1000
  itmax1 <- 1
  if (ipflg > 1)
    itmax1 <- itmax + 1
  ipmax <- itmax * 2
  if (ipflg == 0 || ipflg == 2)
    ipmax <- 1

  z <- .Call("smplxIP",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.integer(skip),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(p),
             as.double(c),
             as.double(uplimit),
             as.double(eps),
             as.integer(itmax),
             as.integer(itmax1),
             as.integer(ipmax),
             as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array(z[[2L]], dim = c(ipmax, n))
  mple <- mples[nip, 1:n]
  names(mple) <- pname
  mples <- mples[1:nip, 1:n]
  if (nip != 1)
    mples <- data.frame(mples)
  names(mples) <- pname

  it1 <- z[[6L]]
  it2 <- 1
  if (ipflg == 2 || ipflg == 3)
    it2 <- it1
  xx <- array(z[[3L]], dim = c(n, itmax1))

  f <- z[[5L]][1:it2]
  para <- xx[1:n, 1:it2]
  std <- z[[4L]][1:it2]
  param <- data.frame(t(para))
  names(param) <- pname

  if (plot == TRUE)
    plot.Estimate("IP", param)

  if (process.report == 0) {
    EstimateIP.out <- list(pa.init = pa, mple = mple, process1 = NULL,
                           process2 = NULL)
  } else if (process.report == 1) {
    EstimateIP.out <- list(pa.init = pa, mple = mple,
                       process1 = list(cflg=ipri, logl = z[[1L]][1:nip],
                       mples = mples), process2 = NULL)
  } else if (process.report == 2) {
    EstimateIP.out <- list(pa.init = pa, mple = mple, process1 = NULL,
                      process2 = list(logl = f, stderr = std, pa.normal = param))
  } else if (process.report == 3) {
    EstimateIP.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = list(logl = f, stderr = std, pa.normal = param))
  }

  class(EstimateIP.out) <- "Estimate"
  return(EstimateIP.out)
}


EstimateThomas <- function(xy.points, pars, eps = 0.001, process.report = 0,
                  plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)
# pa : the parameter values
# (mu, nu, sigma) for the random variable Poisson
  pname <- c("mu", "nu", "sigma") 
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  sigma <- pa[3]
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# parameter
  ipflg <- process.report
  if (process.report == 0 && plot == TRUE)
    ipflg <- 2
  if (process.report == 1 && plot == TRUE)
    ipflg <- 3

  n <- 3
  itmax <- 1000
  itmax1 <- 1
  if (ipflg > 1)
    itmax1 <- itmax + 1
  ipmax <- itmax * 2
  if (ipflg == 0 || ipflg == 2)
    ipmax <- 1

  z <- .Call("smplxThom",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(sigma),
             as.double(eps),
             as.integer(itmax),
             as.integer(itmax1),
             as.integer(ipmax),
             as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array(z[[2L]], dim = c(ipmax, n))
  mple <- mples[nip, 1:n]
  names(mple) <- pname
  mples <- mples[1:nip, 1:n]
  if (nip != 1)
    mples <- data.frame(mples)
  names(mples) <- pname

  it1 <- z[[6L]]
  it2 <- 1
  if (ipflg == 2 || ipflg == 3)
    it2 <- it1
  xx <- array(z[[3L]], dim = c(n, itmax1))

  f <- z[[5L]][1:it2]
  para <- xx[1:n, 1:it2]
  std <- z[[4L]][1:it2]
  param <- data.frame(t(para))
  names(param) <- pname

  if (plot == TRUE)
    plot.Estimate("T", param)

  if (process.report == 0) {
    EstimateThomas.out <- list(pa.init = pa, mple = mple, process1 = NULL,
                               process2=NULL)
  } else if (process.report == 1) {
    EstimateThomas.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = NULL)
  } else if (process.report == 2) {
    EstimateThomas.out <-
      list(pa.init = pa, mple = mple, process1 = NULL,
           process2 = list(logl = f, stderr = std, pa.normal = param))
  } else if (process.report == 3) {
    EstimateThomas.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = list(logl = f, stderr = std, pa.normal = param))
  }

  class(EstimateThomas.out) <- "Estimate"
  return(EstimateThomas.out)
}


EstimateTypeA <- function(xy.points, pars, eps = 0.001, uplimit = 0.3, skip = 1,
                 process.report = 0, plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)
# pa : the parameter values
# (mu, nu, a, sigma1, sigma2) for the random variable Poisson
  pname <- c("mu", "nu", "a", "sigma1", "sigma2") 
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  a <- pa[3]
  sigma1 <- pa[4]
  sigma2 <- pa[5]
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# parameter
  ipflg <- process.report
  if (process.report == 0 && plot == TRUE)
    ipflg <- 2
  if (process.report == 1 && plot == TRUE)
    ipflg <- 3

  n <- 5
  itmax <- 1000
  itmax1 <- 1
  if (ipflg > 1)
    itmax1 <- itmax + 1
  ipmax <- itmax * 2
  if (ipflg == 0 || ipflg == 2)
    ipmax <- 1

  z <- .Call("smplxA", 
             as.double(x),
             as.double(y),
             as.integer(np),
             as.integer(skip),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(a),
             as.double(sigma1),
             as.double(sigma2),
             as.double(uplimit),
             as.double(eps),
             as.integer(itmax),
             as.integer(itmax1),
             as.integer(ipmax),
             as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array(z[[2L]], dim = c(ipmax, n))
  mple <- mples[nip, 1:n]
  names(mple) <- pname
  mples <- mples[1:nip, 1:n]
  if (nip != 1)
    mples <- data.frame(mples)
  names(mples) <- pname

  it1 <- z[[6L]]
  it2 <- 1
  if (ipflg == 2 || ipflg == 3)
    it2 <- it1
  xx <- array(z[[3L]], dim = c(n, itmax1))

  f <- z[[5L]][1:it2]
  para <- xx[1:n, 1:it2]
  std <- z[[4L]][1:it2]
  param <- data.frame(t(para))
  names(param) <- pname

  if (plot == TRUE)
    plot.Estimate("A", param)

  if (process.report == 0) {
    EstimateTypeA.out <- list(pa.init = pa, mple = mple, process1 = NULL,
                              process2 = NULL)
  } else if (process.report == 1) {
    EstimateTypeA.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = NULL)
  } else if (process.report == 2) {
    EstimateTypeA.out <-
      list(pa.init = pa, mple = mple, process1 = NULL,
           process2 = list(logl = f, stderr = std, pa.normal = param))
  } else if (process.report == 3) {
    EstimateTypeA.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = list(logl = f, stderr = std, pa.normal = param))
  }

  class(EstimateTypeA.out) <- "Estimate"
  return(EstimateTypeA.out)
}


EstimateTypeB <- function(xy.points, pars, eps = 0.001, process.report = 0,
                 plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

# pa : the parameter values
# (mu1, nu2, nu, sigma1, sigma2) for the random variable Poisson
  pname = c("mu1", "mu2", "nu", "sigma1", "sigma2")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu1 <- pa[1]
  mu2 <- pa[2]
  nu <- pa[3]
  sigma1 <- pa[4]
  sigma2 <- pa[5]
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# parameter
  ipflg <- process.report
  if (process.report == 0 && plot == TRUE)
    ipflg <- 2
  if (process.report == 1 && plot == TRUE)
    ipflg <- 3

  n <- 5
  itmax <- 1000
  itmax1 <- 1
  if (ipflg > 1)
    itmax1 <- itmax + 1
  ipmax <- itmax * 2
  if (ipflg == 0 || ipflg == 2)
    ipmax <- 1

  z <- .Call("smplxB",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.double(ty),
             as.double(mu1),
             as.double(mu2),
             as.double(nu),
             as.double(sigma1),
             as.double(sigma2),
             as.double(eps),
             as.integer(itmax),
             as.integer(itmax1),
             as.integer(ipmax),
             as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array(z[[2L]], dim = c(ipmax, n))
  pp <- mples[nip, 1:n]
  mple <- c(pp[1]*pp[3], pp[1]*(1-pp[3]), pp[2], pp[4], pp[5])
  names(mple) = pname
  mples <- mples[1:nip, 1:n]
  if (nip != 1)
    mples <- data.frame(mples)
  names(mples) <- c("mu", "nu", "a", "sigma1", "sigma2") 

  it1 <- z[[6L]]
  it2 <- 1
  if (ipflg == 2 || ipflg == 3)
    it2 <- it1

  xx <- array(z[[3L]], dim = c(n, itmax1))
  f <- z[[5L]][1:it2]
  para <- xx[1:n, 1:it2]
  std <- z[[4L]][1:it2]
  param <- data.frame(t(para))
  names(param) <- c("mu", "nu", "a", "sigma1", "sigma2") 

  if (plot == TRUE)
    plot.Estimate("B", param)

#  pini <- c(pa[1], pa[2], pa[3], pa[4], pa[5])

  if (process.report == 0) {
    EstimateTypeB.out <-
      list(pa.init = pa, mple = mple, process1 = NULL, process2 = NULL)
  } else if (process.report == 1) {
    EstimateTypeB.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = NULL)
  } else if (process.report == 2) {
    EstimateTypeB.out <-
      list(pa.init = pa, mple = mple, process1 = NULL,
           process2 = list(logl = f, stderr = std, pa.normal = param))
  } else if (process.report == 3) {
    EstimateTypeB.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = list(logl = f, stderr = std, pa.normal = param))
  }

  class(EstimateTypeB.out) <- "Estimate"
  return(EstimateTypeB.out)
}


EstimateTypeC <- function(xy.points, pars, eps = 0.001, process.report = 0,
                 plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

# (mu1, mu2, nu1, nu2, sigma1, sigma2) for the random variable Poisson
  pname <- c("mu1", "mu2", "nu1", "nu2", "sigma1", "sigma2")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu1 <- pa[1]
  mu2 <- pa[2]
  nu1 <- pa[3]
  nu2 <- pa[4]
  sigma1 <- pa[5]
  sigma2 <- pa[6]
# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
# parameter
  ipflg <- process.report
  if (process.report == 0 && plot == TRUE)
    ipflg <- 2
  if (process.report == 1 && plot == TRUE)
    ipflg <- 3

  n <- 5
  itmax <- 1000
  itmax1 <- 1
  if (ipflg > 1)
    itmax1 <- itmax + 1
  ipmax <- itmax * 2
  if (ipflg == 0 || ipflg == 2)
    ipmax <- 1

 z <- .Call("smplxC",
            as.double(x),
            as.double(y),
            as.integer(np),
            as.double(ty),
            as.double(mu1),
            as.double(mu2),
            as.double(nu1),
            as.double(nu2),
            as.double(sigma1),
            as.double(sigma2),
            as.double(eps),
            as.integer(itmax),
            as.integer(itmax1),
            as.integer(ipmax),
            as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array(z[[2L]], dim = c(ipmax, n))
  pp <- mples[nip, 1:n]
  mple <- c(pp[1] * pp[3] / pp[2], pp[2] * pp[4] / pp[5],  pp[2],
           (1 - pp[3]) * pp[1] * pp[5] / (pp[2] * pp[4]), pp[4], pp[5])
  names(mple) <- pname
  mples <- mples[1:nip, 1:n]
  if (nip != 1)
    mples <- data.frame(mples)
  names(mples) <- c("lambda", "nu1", "a", "sigma1", "sigma2")

  it1 <- z[[6L]]
  it2 <- 1
  if (ipflg == 2 || ipflg == 3)
    it2 <- it1
  xx <- array(z[[3L]], dim = c(n, itmax1))

  f <- z[[5L]][1:it2]
  para <- xx[1:n, 1:it2]
  std <- z[[4L]][1:it2]
  param <- data.frame(t(para))
  names(param) <- c("lambda", "nu1", "a", "sigma1", "sigma2") 

  if (plot == TRUE)
    plot.Estimate("C", param)

#  scllam <- mu1 * nu1 + mu2 * nu2
#  scla <- mu1 * nu1 / scllam
#  pini <- c(mu1, mu2, nu1, nu2, sigma1, sigma2)

  if (process.report == 0) {
    EstimateTypeC.out <-
      list(pa.init = pa, mple = mple, process1 = NULL, process2 = NULL)
  } else if (process.report == 1) {
    EstimateTypeC.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = NULL)
  } else if (process.report == 2) {
    EstimateTypeC.out <-
      list(pa.init = pa, mple = mple, process1 = NULL,
           process2 = list(logl = f, stderr = std, pa.normal = param))
  } else if (process.report == 3) {
    EstimateTypeC.out <-
      list(pa.init = pa, mple = mple,
           process1 = list(cflg = ipri, logl = z[[1L]][1:nip], mples = mples),
           process2 = list(logl = f, stderr = std, pa.normal = param))
  }

  class(EstimateTypeC.out) <- "Estimate"
  return(EstimateTypeC.out)
}


PalmIP <- function(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001,
          uplimit = 0.3, plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

  m <- 0
  mu <- rep(0, 2); nu <- rep(0, 2); p <- rep(0, 2); c <- rep(0, 2)
  pname <- c("mu", "nu", "p", "c")
  lname <- NULL

  if (is.null(pars1) == TRUE) {
      warning("Not enough input parameters. (pars1)", call. = FALSE)
  } else {
    pa1 <- ParsCheck(pars1, pname)
    if (is.null(pa1) == TRUE) {
      warning("Not enough input parameters. (pars1)", call. = FALSE)
    } else {
      m <- 1
      mu[1] <- pa1[1]
      nu[1] <- pa1[2]
      p[1] <- pa1[3]
      c[1] <- pa1[4]
      lname <- c(lname, "true")
    }
  }
  if (is.null(pars2) == TRUE) {
      warning("Not enough input parameters. (pars2)", call. = FALSE)
  } else {
    pa2 <- ParsCheck(pars2, pname)
    if (is.null(pa2)) {
      warning("Not enough input parameters. (pars2)", call. = FALSE)
    } else {
      m <- m + 1
      mu[m] <- pa2[1]
      nu[m] <- pa2[2]
      p[m] <- pa2[3]
      c[m] <- pa2[4]
      lname <- c(lname, "MPLE")
    }
  }

# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1

  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  z <- .Call("palmIP",
              as.double(x),
              as.double(y),
              as.integer(np),
              as.double(delta),
              as.double(ty),
              as.double(uplimit),
              as.double(mu),
              as.double(nu),
              as.double(p),
              as.double(c),
              as.integer(m),
              as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax) * delta
  if (m == 0) {
    palm1 <- NULL
  } else {
    palm1 <- array(z[[2L]], dim = c(jmax, m))
  }
  PalmIP.out <- list(r=r, np.palm = palm, palm.normal = palm1)

  if (plot == TRUE)  {
    plot.Palm("IP", r, palm, palm1, lname)
    return(invisible(PalmIP.out))
  } else {
    class(PalmIP.out) <- "Palm"
    return(PalmIP.out)
  }
}


PalmThomas <- function(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001,
              plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

  m <- 0
  mu <- rep(0, 2); nu <- rep(0, 2); sigma <- rep(0, 2)
  pname <- c("mu", "nu", "sigma")
  lname <- NULL

  if (is.null(pars1) == TRUE) {
    warning("Not enough input parameters. (pars1)", call. = FALSE)
  } else {
    pa1 <- ParsCheck(pars1, pname)
    if (is.null(pa1)) {
      warning("Not enough input parameters. (pars1)", call. = FALSE)
    } else {
      m <- 1
      mu[1] <- pa1[1]
      nu[1] <- pa1[2]
      sigma[1] <- pa1[3]
      lname <- c(lname, "true")
    }
  }
  if (is.null(pars2) == TRUE) {
    warning("Not enough input parameters. (pars2)", call. = FALSE)
  } else {
    pa2 <- ParsCheck(pars2, pname)
    if (is.null(pa2)) {
      warning("Not enough input parameters. (pars2)", call. = FALSE)
    } else {
      m <- m + 1
      mu[m] <- pa2[1]
      nu[m] <- pa2[2]
      sigma[m] <- pa2[3]
      lname <- c(lname, "MPLE")
    }
  }

# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1

  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  z <- .Call("palmT",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.double(delta),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(sigma),
             as.integer(m),
             as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax) * delta
  if (m == 0) {
    palm1 <- NULL
  } else {
    palm1 <- array(z[[2L]], dim = c(jmax, m))
  }

  PalmThomas.out <- list(r = r, np.palm = palm, palm.normal = palm1)

  if (plot == TRUE)  {
    plot.Palm("T", r, palm, palm1, lname)
    return(invisible(PalmThomas.out))
  } else {
    class(PalmThomas.out) <- "Palm"
    return(PalmThomas.out)
  }
}


PalmTypeA <- function(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001,
             uplimit = 0.3, plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

  m <- 0
  mu <- rep(0, 2); nu <- rep(0, 2); a <- rep(0,2)
  sigma1 <- rep(0, 2); sigma2 <- rep(0,2)
  pname <- c("mu", "nu", "a", "sigma1", "sigma2")
  lname <- NULL

  if (is.null(pars1) == TRUE) {
    warning("Not enough input parameters. (pars1)", call. = FALSE)
  } else {
    pa1 <- ParsCheck(pars1, pname)
    if (is.null(pa1)) {
      warning("Not enough input parameters. (pars1)", call. = FALSE)
    } else {
      m <- 1
      mu[1] <- pa1[1]
      nu[1] <- pa1[2]
      a[1] <- pa1[3]
      sigma1[1] <- pa1[4]
      sigma2[1] <- pa1[5]
      lname <- c(lname, "true")
    }
  }
  if (is.null(pars2) == TRUE) {
    warning("Not enough input parameters. (pars2)", call. = FALSE)
  } else {
    pa2 <- ParsCheck(pars2, pname)
    if (is.null(pa2)) {
      warning("Not enough input parameters. (pars2)", call. = FALSE)
    } else {
      m <- m + 1
      mu[m] <- pa2[1]
      nu[m] <- pa2[2]
      a[m] <- pa2[3]
      sigma1[m] <- pa2[4]
      sigma2[m] <- pa2[5]
      lname <- c(lname, "MPLE")
    }
  }

# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1

  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  z <- .Call("palmA",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.double(delta),
             as.double(ty),
             as.double(uplimit),
             as.double(mu),
             as.double(nu),
             as.double(a),
             as.double(sigma1),
             as.double(sigma2),
             as.integer(m),
             as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax) * delta
  if (m == 0) {
    palm1 <- NULL
  } else {
    palm1 <- array(z[[2L]], dim = c(jmax, m))
  }

  PalmTypeA.out <- list(r = r, np.palm = palm, palm.normal = palm1)

  if (plot == TRUE)  {
    plot.Palm("A", r, palm, palm1, lname)
    return(invisible(PalmTypeA.out))
  } else {
    class(PalmTypeA.out) <- "Palm"
    return(PalmTypeA.out)
  }
}


PalmTypeB <- function(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001,
             plot = TRUE){

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

  m <- 0
  mu <- rep(0, 2); nu <- rep(0, 2); a <- rep(0,2)
  sigma1 <- rep(0, 2); sigma2 <- rep(0,2)
  pname <- c("mu1", "mu2", "nu", "sigma1", "sigma2")
  lname <- NULL

  if (is.null(pars1) == TRUE) {
    warning("Not enough input parameters. (pars1)", call. = FALSE)
  } else {
    pa1 <- ParsCheck(pars1, pname)
    if (is.null(pa1)) {
      warning("Not enough input parameters. (pars1)", call. = FALSE)
    } else {
      m <- 1
      mu[1] <- pa1[1] + pa1[2]
      nu[1] <- pa1[3]
      a[1] <- pa1[1] / mu[1]
      sigma1[1] <- pa1[4]
      sigma2[1] <- pa1[5]
      lname <- c(lname, "true")
    }
  }
  if (is.null(pars2) == TRUE) {
    warning("Not enough input parameters. (pars2)", call. = FALSE)
  } else {
    pa2 <- ParsCheck(pars2, pname)
    if (is.null(pa2)) {
      warning("Not enough input parameters. (pars2)", call. = FALSE)
    } else {
      m <- m + 1
      mu[m] <- pa2[1] + pa2[2]
      nu[m] <- pa2[3]
      a[m] <- pa2[1] / mu[m]
      sigma1[m] <- pa2[4]
      sigma2[m] <- pa2[5]
      lname <- c(lname, "MPLE")
    }
  }

# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1

  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  z <- .Call("palmB",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.double(delta),
             as.double(ty),
             as.double(mu),
             as.double(nu),
             as.double(a),
             as.double(sigma1),
             as.double(sigma2),
             as.integer(m),
             as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax) * delta
  if (m == 0) {
    palm1 <- NULL
  } else {
    palm1 <- array(z[[2L]], dim = c(jmax, m))
  }

  PalmTypeB.out <- list(r = r, np.palm = palm, palm.normal = palm1)

  if (plot == TRUE)  {
    plot.Palm("B", r, palm, palm1, lname)
    return(invisible(PalmTypeB.out))
  } else {
    class(PalmTypeB.out) <- "Palm"
    return(PalmTypeB.out)
  }
}


PalmTypeC <- function(xy.points, pars1 = NULL, pars2 = NULL, delta = 0.001,
             plot = TRUE) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

  m <- 0
  lamb <- rep(0, 2); nu1 <- rep(0, 2); a <- rep(0,2)
  sigma1 <- rep(0, 2); sigma2 <- rep(0,2)
  pname <- c("mu1", "mu2", "nu1", "nu2", "sigma1", "sigma2")
  lname<- NULL

  if (is.null(pars1) == TRUE) {
    warning("Not enough input parameters. (pars1)", call. = FALSE)
  } else {
    pa1 <- ParsCheck(pars1, pname)
    if (is.null(pa1)) {
      warning("Not enough input parameters. (pars1)", call. = FALSE)
    } else {
      m <- 1
      lamb[1] <- pa1[1] * pa1[3] + pa1[2] * pa1[4]
      nu1[1] <- pa1[3]
      a[1] <- (pa1[1] * pa1[3])/ lamb[1]
      sigma1[1] <- pa1[5]
      sigma2[1] <- pa1[6]
      lname <- c(lname, "true")
    }
  }
  if (is.null(pars2) == TRUE) {
    warning("Not enough input parameters. (pars2)", call. = FALSE)
  } else {
    pa2 <- ParsCheck(pars2, pname)
    if (is.null(pa2)) {
      warning("Not enough input parameters. (pars2)", call. = FALSE)
    } else {
      m <- m + 1
      lamb[m] <- pa2[1] * pa2[3] + pa2[2] * pa2[4]
      nu1[m] <- pa2[3]
      a[m] <- (pa2[1] * pa2[3]) / lamb[m]
      sigma1[m] <- pa2[5]
      sigma2[m] <- pa2[6]
      lname <- c(lname, "MPLE")
    }
  }


# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1

  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  z <- .Call("palmC",
             as.double(x),
             as.double(y),
             as.integer(np),
             as.double(delta),
             as.double(ty),
             as.double(lamb),
             as.double(nu1),
             as.double(a),
             as.double(sigma1),
             as.double(sigma2),
             as.integer(m),
             as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax) * delta
  if (m == 0) {
    palm1 <- NULL
  } else {
    palm1 <- array(z[[2L]], dim = c(jmax, m))
  }

  PalmTypeC.out <- list(r = r, np.palm = palm, palm.normal = palm1)

  if (plot == TRUE)  {
    plot.Palm("C", r, palm, palm1, lname)
    return(invisible(PalmTypeC.out))
  } else {
    class(PalmTypeC.out) <- "Palm"
    return(PalmTypeC.out)
  }
}


##### Input Parameters Check for the Model #####

ParsCheck <- function(param, pname) {
  idsum <- c(1, 3, 6, 10, 15, 21)
  npa <- length(param)
  nn <-length(pname)
  if (npa != nn)
    return(NULL)
  pname.in <- names(param)
  if (is.null(pname.in) == TRUE) {
    return(param)
  } else {
    z <- pname.in[pname.in %in% pname]
    if (length(z) != nn) {
     return(NULL)
    } else { # length(z) == nn
      id <- rep(0,nn)
      for (i in 1:nn)
        id[i] <- which(pname.in == pname[i])
      if (sum(id) != idsum[nn])
        return(NULL)
      param1 <- c(param[id[1:nn]])
      return(param1)
    }
  }
}


##### Print and Plot  #####

print.Simulate <- function(x, ...) {
  np <- length(x$parents$n)
  if (np == 1) {
    np1 <- x$parents$n
    np2 <- x$offspring$n
    cat(sprintf("\nNumber of parents =  %8.1f\n", np1))
    cat(sprintf("Total number of offspring =%8.1f\n", np2))
  } else if (np == 2) {
    np1 <- x$parents$n[1]
    np2 <- x$parents$n[2]
    mp1 <- x$offspring$n[1]
    mp2 <- x$offspring$n[2]
    cat(sprintf("\nNumber of parents =  %8.1f\t%8.1f\tTotal%8.1f\n",
                np1, np2, np1 + np2))
    cat(sprintf("Number of offspring =%8.1f\t%8.1f\tTotal%8.1f\n",
                mp1, mp2, mp1 + mp2))
  }
}

plot.Simulate <- function(type, parents, offspring, pa){
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), pch = 20, cex = 0.75, xaxs = "i", yaxs = "i", pty = "s")

  if (type == "T") {
    mtitle <- "Simulation of Thomas model"
    stitle1 <- substitute(paste(mu == v1), list(v1 = pa[1]))
    stitle2 <- substitute(
      paste("(",mu,", ",nu,", ",sigma,")" == "(",v1,", ",v2,", ",v3,")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3]))
  } else if (type == "IP") {
    mtitle <- "Simulation of inverse-power type model"
    stitle1 <- substitute(paste(mu == v1), list(v1 = pa[1]))
    stitle2 <- substitute(
      paste("(",mu,", ",nu,", p, c)" == "(",v1,", ",v2,", ",v3,", ",v4,")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4]))
  } else if (type == "A") {
    mtitle <- "Simulation of Type A model"
    stitle1 <- substitute(paste(mu == v1), list(v1 = pa[1]))
    stitle2 <- substitute(
      paste("(",mu,", ",nu,", ",a,", ",sigma[1],", ",sigma[2],")"
         == "(",v1,", ",v2,", ",v3,", ",v4,", ",v5,")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4], v5 = pa[5]))
  } else if (type == "B") {
    mtitle <- "Simulation of Type B model"
    stitle1 <- substitute(
      paste("(", mu[1], ", ", mu[2], ")" == "(", v1, ", ", v2, ")"),
      list(v1=pa[1], v2=pa[2]))
    stitle2 <- substitute(
      paste("(", mu[1],", ",mu[2],", ",nu,", ",sigma[1],", ",sigma[2],")"
         == "(", v1, ", ", v2, ", ", v3, ", ", v4, ", ", v5, ")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4], v5 = pa[5]))
  } else if (type == "C") {
    mtitle <- "Simulation of Type C model"
    stitle1 <- substitute(
      paste("(", mu[1], ", ", mu[2],")" == "(", v1, ", ", v2, ")"),
      list(v1 = pa[1], v2 = pa[2]))
    stitle2 <- substitute(
      paste("(", mu[1], ", ", mu[2], ", ", nu[1], ", ", nu[2], ", ", sigma[1],
 ", ", sigma[2],")"
         == "(", v1, ", ", v2, ", ", v3, ", ", v4, ", ", v5, ", ", v6, ")"),
      list(v1=pa[1], v2=pa[2], v3=pa[3], v4=pa[4], v5=pa[5], v6=pa[6]))
  }

  plot(parents, xlim = c(0, 1), ylim = c(0, 1),
       main = "Simulation of parent points", sub = stitle1, xlab = "",
       ylab = "", cex = 0.8)
  plot(offspring, xlim = c(0,1), ylim = c(0,1), main = mtitle, sub = stitle2,
       xlab = "", ylab = "", cex = 0.7)

  par(mfrow = old.par$mfrow, pch = old.par$pch, cex = old.par$cex,
      xaxs = old.par$xaxs, yaxs = old.par$yaxs, pty = old.par$pty)
}

plot.Simulate2 <- function(type, parents, offspring, pa){
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(1,2), pch = 20, cex = 0.75, xaxs = "i", yaxs = "i", pty = "s")

  if (type == "B") {
    mtitle <- "Simulation of Type B model"
    stitle1 <- substitute(
      paste("(", mu[1], ", ", mu[2], ")" == "(", v1, ", ", v2, ")"),
      list(v1 = pa[1], v2 = pa[2]))
    stitle2 <- substitute(
      paste("(", mu[1],", ",mu[2],", ",nu,", ",sigma[1],", ",sigma[2],")"
         == "(", v1, ", ", v2, ", ", v3, ", ", v4, ", ", v5, ")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4], v5 = pa[5]))
  } else if (type == "C") {
    mtitle <- "Simulation of Type C model"
    stitle1 <- substitute(
      paste("(", mu[1], ", ", mu[2],")" == "(", v1, ", ", v2, ")"),
      list(v1 = pa[1], v2 = pa[2]))
    stitle2 <- substitute(
      paste("(", mu[1], ", ", mu[2], ", ", nu[1], ", ", nu[2], ", ", sigma[1], 
", ", sigma[2],")"
         == "(", v1, ", ", v2, ", ", v3, ", ", v4, ", ", v5, ", ", v6, ")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4],
           v5 = pa[5], v6 = pa[6]))
  }

  plot(parents$xy1, xlim = c(0, 1), ylim = c(0, 1),
       main = "Simulation of parent points", sub = stitle1,
       xlab = "", ylab = "", cex = 0.8, col = "brown2")
  par(new = TRUE)
  plot(parents$xy2, xlim = c(0, 1), ylim = c(0, 1), main = "", sub = "",
       xlab = "", ylab = "", cex = 0.8, col = "royalblue3")

  plot(offspring$xy1, xlim = c(0, 1), ylim = c(0, 1), main = mtitle,
       sub = stitle2, xlab = "", ylab = "", cex = 0.7, col = "brown2")
  par(new = TRUE)
  plot(offspring$xy2, xlim = c(0, 1), ylim = c(0, 1), main = "", sub = "",
       xlab = "", ylab = "", cex = 0.7, col = "royalblue3")

  par(mfrow = old.par$mfrow, pch = old.par$pch, cex = old.par$cex,
      xaxs = old.par$xaxs, yaxs = old.par$yaxs, pty = old.par$pty,
      new = old.par$new)
}


print.Palm <- function(x, ...) {
  cat("\n\t\t Palm intensity function\n")
  cat("\tr\t non-parametric\t     normalized Palm\n")
  n <- length(x$r)
  m <- dim(x$palm.normal)[2]
  for (i in 1:n) {
    cat(sprintf("\n %10.3f\t%12.5f\t", x$r[i], x$np.palm[i]))
    for (j in 1:m)
      cat(sprintf("\t%12.5f", x$palm.normal[i, j]))
  }
  cat("\n\n")
}

plot.Palm <- function(type, r, palm, palm1, lname) {
  if (is.null(palm1)) {
    ymax <- max(palm)*2
  } else {
    ymax <- max(palm, palm1) * 2
    m <-  dim(palm1)[2]
  }

 switch(type,
   "T" =  mtitle <- "Thomas model",
   "IP" = mtitle <- "Inverse-power type model",
   "A" = mtitle <- "Type A model",
   "B" = mtitle <- "Type B model",
   "C" = mtitle <- "Type C model",
 )

  plot(x = r, y = palm, pch = 20, ylim = c(0.5, ymax), log = "xy",
       main = mtitle, xlab = "r", ylab = "Palm intensity", cex.main = 2.0,
       cex.lab = 1.5)
  abline(h = 1)
  if (is.null(palm1)) {
    legend("topright", legend = c("non-parametric"), pch = 20, cex = 1.2)
  } else {
    lcol <- c(1, 3, 2, 4, 5, 6, 7, 8, 9)
    leg.txt <- c("non-parametric", lname)
    for (j in 1:m) {
      par(new = TRUE)
      plot(x = r, y = palm1[, j], type = "l", ylim = c(0.5, ymax),
           col = lcol[j + 1], log = "xy", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n")
    }
    legend("topright", legend = leg.txt[1:(m + 1)], col = lcol[1:(m + 1)],
           pch = c(20, rep(NA, m)), lty = c(NA, rep(1, m)), cex = 1.2)
  }
}


print.Estimate <- function(x, ...) {
  n <- length(x$pa.init)

  if (length(x$process1) != 0) {
    cflg <- x$process1$cflg
    logl <- x$process1$logl
    mples <- x$process1$mples
    pname1 <- names(mples)
    nn <- length(pname1)
    cat("\n\n Minimizing the negative Palm log-likelihood function\n")
    cat("\t\t-log L\t")
    for (i in 1:nn)
      cat(sprintf("\t%12s", pname1[i]))
    cat("\n") 
    ipri2 <- length(x$process1$logl)
    for( j in 1:ipri2 ) {
      if (cflg[j] == -1)
        cat(sprintf(" testfn = %18.10e", logl[j]))
      if (cflg[j] == 1)
        cat(sprintf(" update = %18.10e", logl[j]))
      for (i in 1:nn)
        cat(sprintf("\t%12.4e", mples[j, i]))
      cat("\n")
    }
  }

  if (length(x$process2) != 0) {
    logl <- x$process2$logl
    stderr <- x$process2$stderr
    pa <- x$process2$pa.normal
    pname2 <- names(pa)
    nn <- length(pname2)
    cat("\n\n
      Optimization procedure by the simplex with the normalized parameters\n")
    cat(" # \t\t -logL\t\t standard error")
    for (i in 1:nn)
      cat(sprintf("\t%12s", pname2[i]))
    cat("\n")
    it2 <- length(logl)
    for (j in 1:it2) {
      cat(sprintf(" %i\t %15.3f\t %12.4e", j - 1, logl[j], stderr[j]))
      for (i in 1:nn)
        cat(sprintf("\t %12.4e", pa[j, i]))
      cat("\n")
    }
  }

  pname <- names(x$mple)
  cat("\n\n Parameter")
  for (i in 1:n)
     cat(sprintf("\t%12s", pname[i]))
  cat("\n")
  cat(" Initial value")
  for ( i in 1:n)
    cat(sprintf("\t%12.4f", x$pa.init[i]))
  cat("\n")
  cat(" MPLE    ")
  for (i in 1:n)
    cat(sprintf("\t%12.4f", x$mple[i]))
  cat("\n\n")
}


plot.Estimate <- function(type, para) {

  pname <- names(para)
  n <- dim(para)[2]
  n1 <- n + 1

  old.par <- par(no.readonly = TRUE)
  new.mai <- old.par$mai
  new.mai[2] <- new.mai[1] * 1.2
  par(mai = new.mai, xaxs = "i")

 switch(type,
   "T" =  mtitle <- "Thomas model",
   "IP" = mtitle <- "Inverse-power type model",
   "A" = mtitle <- "Type A model",
   "B" = mtitle <- "Type B model",
   "C" = mtitle <- "Type C model",
 )

  plot(para[, 1], type = "l", ylim = c(0.25, 1.75), main = mtitle,
       xlab = "Number of iterations",
       ylab = "Values of the normalized parameters",
       col = 2, cex.main = 2.0, cex.lab = 1.5)

  for (i in 2:n) {
    par(new = TRUE)
    plot(para[, i], type = "l", ylim = c(0.25, 1.75), main = "",
         xlab = "", ylab = "", col = i + 1)
  }
  legend("topright", legend = pname, col = c(2:n1), lty = 1, cex = 1.2)

  par(mai = old.par$mai, xaxs = old.par$xaxs, new = old.par$new)
}

