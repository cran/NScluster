sim.cppm <- function(model = "Thomas", pars, seed = NULL) {

# seed : initial seeds for a sequence of uniform random numbers
  ix <- seed
  if (is.null(ix))
    ix <- -1

# ty : the variable for the standardized coordinates of points
#  in the rectangular region
  ty <- 1

# maximum number of parent points / maximum number of offspring points
##  pmax <- 100
##  omax <- 100
  pmax <- 500
  omax <- 500

  if (model == "IP") {

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
  pars1 <- c(mu, nu, p, c)

  z <- .Fortran(C_simip,
                as.integer(ix),
                as.double(ty),
                as.double(mu),
                as.double(nu),
                as.double(p),
                as.double(c),
                npts = integer(1),
                ncl = integer(pmax),
                x = double(pmax),
                y = double(pmax),
                xcl = double(pmax*omax),
                ycl = double(pmax*omax),
                as.integer(pmax),
                as.integer(omax),
                ier = integer(1))

  } else   if (model == "Thomas") {

# the parameter values
# (mu, nu, sigma) for the random variable Poisson
  pname <- c("mu", "nu", "sigma")
  pa <- ParsCheck(pars, pname)
  if (is.null(pa))
    stop("Not enough input parameters. (pars)", call. = FALSE)
  mu <- pa[1]
  nu <- pa[2]
  sigma <- pa[3]
  pars1 <- c(mu, nu, sigma)

  z <- .Fortran(C_simthom,
                as.integer(ix),
                as.double(ty),
                as.double(mu),
                as.double(nu),
                as.double(sigma),
                npts = integer(1),
                ncl = integer(pmax),
                x = double(pmax),
                y = double(pmax),
                xcl = double(pmax*omax),
                ycl = double(pmax*omax),
                as.integer(pmax),
                as.integer(omax),
                ier = integer(1))

  } else if (model == "TypeA") {

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
  pars1 <- c(mu, nu, a, sigma1, sigma2)

  z <- .Fortran(C_sima,
                as.integer(ix),
                as.double(ty),
                as.double(mu),
                as.double(nu),
                as.double(a),
                as.double(sigma1),
                as.double(sigma2),
                npts = integer(1),
                ncl = integer(pmax),
                x = double(pmax),
                y = double(pmax),
                xcl = double(pmax*omax),
                ycl = double(pmax*omax),
                as.integer(pmax),
                as.integer(omax),
                ier = integer(1))

  } else if (model == "TypeB") {

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
  pars1 <- c(mu1, mu2, nu, sigma1, sigma2)

  z <- .Fortran(C_simb,
                as.integer(ix),
                as.double(ty),
                as.double(mu1),
                as.double(mu2),
                as.double(nu),
                as.double(sigma1),
                as.double(sigma2),
                m1 = integer(1),
                ncl1 = integer(pmax),
                x1 = double(pmax),
                y1 = double(pmax),
                xx1 = double(pmax*omax),
                yy1 = double(pmax*omax),
                m2 = integer(1),
                ncl2 = integer(pmax),
                x2 = double(pmax),
                y2 = double(pmax),
                xx2 = double(pmax*omax),
                yy2 = double(pmax*omax),
                as.integer(pmax),
                as.integer(omax),
                ier = integer(1))

  } else if (model == "TypeC") {

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
  pars1 <- c(mu1, mu2, nu1, nu2, sigma1, sigma2)

# maximum number of parent points / maximum number of offspring points
  pmax <- pmax * 3
  omax <- omax * 3

  z <- .Fortran(C_simc,
                as.integer(ix),
                as.double(ty),
                as.double(mu1),
                as.double(mu2),
                as.double(nu1),
                as.double(nu2),
                as.double(sigma1),
                as.double(sigma2),
                m1 = integer(1),
                ncl1 = integer(pmax),
                x1 = double(pmax),
                y1 = double(pmax),
                xx1 = double(pmax*omax),
                yy1 = double(pmax*omax),
                m2 = integer(1),
                ncl2 = integer(pmax),
                x2 = double(pmax),
                y2 = double(pmax),
                xx2 = double(pmax*omax),
                yy2 = double(pmax*omax),
                as.integer(pmax),
                as.integer(omax),
                ier = integer(1))

  } else {
    stop("the model type is invalid.")
  }
  
  if (model != "TypeB"  && model != "TypeC") {  # IP, Thomas, TypeA

    ier <- z$ier
    if (ier == -1) {
      stop(paste("too many parents (the default maximum number is ", pmax, ")"),
           call. = FALSE)
    } else if (ier == -2) {
      stop(paste(
	       "too many offspring (the default maximum number is ", omax, ")"),
           call. = FALSE)
    } else {

      npts <- z$npts
      ncl <- z$ncl
      parents.x <- z$x[1:npts]
      parents.y <- z$y[1:npts]
      xcl <- array(z$xcl, dim = c(pmax, omax))
      ycl <- array(z$ycl, dim = c(pmax, omax))

      parents.xy <- array(c(parents.x, parents.y), dim = c(npts, 2))

      offspring.x <- NULL
      offspring.y <- NULL
      for (i in 1:npts) {
        offspring.x <- c(offspring.x, xcl[i, 1:ncl[i]])
        offspring.y <- c(offspring.y, ycl[i, 1:ncl[i]])
      }
      np <- sum(ncl[1:npts])
      offspring.xy <- array(c(offspring.x, offspring.y), dim = c(np, 2))

      out <- list(parents = list(n = npts, xy = parents.xy),
                    offspring = list(n = np, xy = offspring.xy))
    }

  } else {  # TypeB, TypeC 

    ier <- z$ier
    if (ier == -1 || ier == -2) {
      stop(paste("too many parents (the default maximum number is ", pmax, ")"),
           call. = FALSE)
    } else if (ier == -11 || ier == -22) {
      stop(paste(
	       "too many offspring (the default maximum number is ", omax, ")"),
           call. = FALSE)
    } else {

      m1 <- z$m1
      m2 <- z$m2
      parents1.x <- z$x1[1:m1]
      parents1.y <- z$y1[1:m1]
      parents2.x <- z$x2[1:m2]
      parents2.y <- z$y2[1:m2]
      parents.xy <- array(0, dim = c(m1 + m2, 2))
      parents.xy[, 1] <- c(parents1.x, parents2.x)
      parents.xy[, 2] <- c(parents1.y, parents2.y)
      np1 <- sum(z$ncl1[1:m1])
      offspring1.x <- z$xx1[1:np1]
      offspring1.y <- z$yy1[1:np1]
      np2 <- sum(z$ncl2[1:m2])
      offspring2.x <- z$xx2[1:np2]
      offspring2.y <- z$yy2[1:np2]
      offspring.xy <- array(0, dim = c(np1 + np2, 2))
      offspring.xy[, 1] <- c(offspring1.x, offspring2.x)
      offspring.xy[, 2] <- c(offspring1.y, offspring2.y)

      out <- list(parents = list(n = c(m1, m2), xy = parents.xy),
                  offspring = list(n = c(np1, np2), xy = offspring.xy))
    }
  }
  
  out <- c(list(model = model, pars = pars1), out)
  class(out) <- c("sim.cpp")
  out

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

##### S3method #####

print.sim.cpp <- function(x, ...) {
  np <- length(x$parents$n)
  if (np == 1) {
    np1 <- x$parents$n
    np2 <- x$offspring$n
    cat(sprintf("\nNumber of parent points = %8.1f\n", np1))
    cat(sprintf("Total number of offspring points = %8.1f\n", np2))
  } else if (np == 2) {
    np1 <- x$parents$n[1]
    np2 <- x$parents$n[2]
    mp1 <- x$offspring$n[1]
    mp2 <- x$offspring$n[2]
    cat(sprintf("\nNumber of parent points = %8.1f\t%8.1f\tTotal%8.1f\n",
                np1, np2, np1 + np2))
    cat(sprintf("Number of offspring points = %8.1f\t%8.1f\tTotal%8.1f\n",
                mp1, mp2, mp1 + mp2))
  }
}

plot.sim.cpp <- function(x, parents.distinct = FALSE, ...) {

  model <- x$model
  parents <- x$parents$xy
  offspring <- x$offspring$xy
  pa <- x$pars

  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), pch = 20, cex = 0.75, xaxs = "i", yaxs = "i", pty = "s")


  if (model == "Thomas") {
    mtitle <- "Simulation of Thomas model"
    stitle1 <- substitute(paste(mu == v1), list(v1 = pa[1]))
#    stitle1 <- substitute(paste(mu == v1),
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1)))
    stitle2 <- substitute(
      paste("(",mu,", ",nu,", ",sigma,")" == "(",v1,", ",v2,", ",v3,")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1),
#           v3 = format(round(pa[3], digits = 3), nsmall = 3)))

  } else if (model == "IP") {
    mtitle <- "Simulation of Inverse-power type model"
    stitle1 <- substitute(paste(mu == v1), list(v1 = pa[1]))
#    stitle1 <- substitute(paste(mu == v1),
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1)))
    stitle2 <- substitute(
      paste("(",mu,", ",nu,", p, c)" == "(",v1,", ",v2,", ",v3,", ",v4,")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1),
#           v3 = format(round(pa[3], digits = 1), nsmall = 1),
#           v4 = format(round(pa[4], digits = 3), nsmall = 3)))

  } else if (model == "TypeA") {
    mtitle <- "Simulation of Type A model"
    stitle1 <- substitute(paste(mu == v1), list(v1 = pa[1]))
#    stitle1 <- substitute(paste(mu == v1),
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1)))
    stitle2 <- substitute(
      paste("(",mu,", ",nu,", ",a,", ",sigma[1],", ",sigma[2],")"
         == "(",v1,", ",v2,", ",v3,", ",v4,", ",v5,")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4], v5 = pa[5]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1),
#           v3 = format(round(pa[3], digits = 1), nsmall = 1),
#           v4 = format(round(pa[4], digits = 3), nsmall = 3),
#           v5 = format(round(pa[5], digits = 3), nsmall = 3)))

  } else if (model == "TypeB") {
    mtitle <- "Simulation of Type B model"
    stitle1 <- substitute(
      paste("(", mu[1], ", ", mu[2], ")" == "(", v1, ", ", v2, ")"),
      list(v1=pa[1], v2=pa[2]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1)))
    stitle2 <- substitute(
      paste("(", mu[1],", ",mu[2],", ",nu,", ",sigma[1],", ",sigma[2],")"
         == "(", v1, ", ", v2, ", ", v3, ", ", v4, ", ", v5, ")"),
      list(v1 = pa[1], v2 = pa[2], v3 = pa[3], v4 = pa[4], v5 = pa[5]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1),
#           v3 = format(round(pa[3], digits = 1), nsmall = 1),
#           v4 = format(round(pa[4], digits = 3), nsmall = 3),
#           v5 = format(round(pa[5], digits = 3), nsmall = 3)))

  } else if (model == "TypeC") {
    mtitle <- "Simulation of Type C model"
    stitle1 <- substitute(
      paste("(", mu[1], ", ", mu[2],")" == "(", v1, ", ", v2, ")"),
      list(v1 = pa[1], v2 = pa[2]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1)))
    stitle2 <- substitute(
      paste("(", mu[1], ", ", mu[2], ", ", nu[1], ", ", nu[2], ", ", sigma[1],
 ", ", sigma[2],")"
         == "(", v1, ", ", v2, ", ", v3, ", ", v4, ", ", v5, ", ", v6, ")"),
      list(v1=pa[1], v2=pa[2], v3=pa[3], v4=pa[4], v5=pa[5], v6=pa[6]))
#      list(v1 = format(round(pa[1], digits = 1), nsmall = 1),
#           v2 = format(round(pa[2], digits = 1), nsmall = 1),
#           v3 = format(round(pa[3], digits = 1), nsmall = 1),
#           v4 = format(round(pa[4], digits = 1), nsmall = 1),
#           v5 = format(round(pa[5], digits = 3), nsmall = 3),
#           v6 = format(round(pa[6], digits = 3), nsmall = 3)))
  }

  if (model != "TypeB" && model !="TypeC")
    parents.distinct <- FALSE

  if (parents.distinct == FALSE) {

      plot(parents, xlim = c(0, 1), ylim = c(0, 1),
           main = "Simulation of parent points", sub = stitle1, xlab = "",
           ylab = "", cex = 0.8, ...)

      plot(offspring, xlim = c(0,1), ylim = c(0,1), main = mtitle,
           sub = stitle2, xlab = "", ylab = "", cex = 0.7, ...)

      par(mfrow = old.par$mfrow, pch = old.par$pch, cex = old.par$cex,
          xaxs = old.par$xaxs, yaxs = old.par$yaxs, pty = old.par$pty)

  } else {

      np1 <- x$parents$n[1]
      np <- np1 + x$parents$n[2]
      parents1 <- parents[1:np1, ]
      parents2 <- parents[(np1+1):np, ]

      plot(parents1, xlim = c(0, 1), ylim = c(0, 1),
           main = "Simulation of parent points", sub = stitle1,
       xlab = "", ylab = "", cex = 0.8, col = "brown2", ...)
       par(new = TRUE)
       plot(parents2, xlim = c(0, 1), ylim = c(0, 1), main = "", sub = "",
           xlab = "", ylab = "", cex = 0.8, col = "royalblue3", ...)

      mp1 <- x$offspring$n[1]
      mp <- mp1 + x$offspring$n[2]
      offspring1 <- offspring[1:mp1, ]
      offspring2 <- offspring[(mp1+1):mp, ]

      plot(offspring1, xlim = c(0, 1), ylim = c(0, 1), main = mtitle,
           sub = stitle2, xlab = "", ylab = "", cex = 0.7, col = "brown2", ...)
      par(new = TRUE)
      plot(offspring2, xlim = c(0, 1), ylim = c(0, 1), main = "", sub = "",
           xlab = "", ylab = "", cex = 0.7, col = "royalblue3", ...)

      par(mfrow = old.par$mfrow, pch = old.par$pch, cex = old.par$cex,
          xaxs = old.par$xaxs, yaxs = old.par$yaxs, pty = old.par$pty,
          new = old.par$new)

  }

}

