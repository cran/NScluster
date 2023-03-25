mple.cppm <- function(model = "Thomas", xy.points, pars = NULL, eps = 0.001,
                      uplimit = 0.3, skip = 1) {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)

# ty : the variable for the standardized coordinates of points
  ty <- 1
# parameter
  ipflg <- 3
  itmax <- 1000
  itmax1 <- itmax + 1
  ipmax <- itmax * 2

  switch(model,
    "IP" = pname <- c("mu", "nu", "p", "c"),
    "Thomas" =    pname <- c("mu", "nu", "sigma"),
    "TypeA" =   pname <- c("mu", "nu", "a", "sigma1", "sigma2"),
    "TypeB" =   pname <- c("mu1", "mu2", "nu", "sigma1", "sigma2"),
    "TypeC" =   pname <- c("mu1", "mu2", "nu1", "nu2", "sigma1", "sigma2"),
    stop("the model type is invalid.")
  )

  if (is.null(pars) == TRUE) {
    pa <- set.pars(model, xy.points)
  } else {
    pa <- ParsCheck(pars, pname)
    if (is.null(pa))
      stop("Not enough input parameters. (pars)", call. = FALSE)
  }

  if (model == "Thomas") {

# (mu, nu, sigma) for the random variable Poisson

    pname.org <- pname
    n <- 3
    mu <- pa[1]
    nu <- pa[2]
    sigma <- pa[3]

    z <- .Fortran(C_smplxthom,
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
                  fn = double(ipmax),
                  mple = double(n*ipmax),
                  xx = double(n*itmax1),
                  std = double(itmax1),
                  f = double(itmax1),
                  itr = integer(1),
                  nip = integer(1),
                  ipr = integer(ipmax),
                  as.integer(ipflg))

    nip <- z$nip
    mples <- array(z$mple, dim = c(ipmax, n))
    mple <- mples[nip, 1:n]

  } else if (model == "IP") {

# (mu, nu) for the random variable Poisson
# p : a decay order with respect to the distance
# c : scaling factor with respect to the distance

    pname.org <- pname
    n <- 4
    mu <- pa[1]
    nu <- pa[2]
    p <- pa[3]
    c <- pa[4]

    z <- .Fortran(C_smplxip,
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
                  fn = double(ipmax),
                  mple = double(n*ipmax),
                  xx = double(n*itmax1),
                  std = double(itmax1),
                  f = double(itmax1),
                  itr = integer(1),
                  nip = integer(1),
                  ipr = integer(ipmax),
                  as.integer(ipflg))

    nip <- z$nip
    mples <- array(z$mple, dim = c(ipmax, n))
    mple <- mples[nip, 1:n]

  } else if (model == "TypeA") {

# (mu, nu, a, sigma1, sigma2) for the random variable Poisson

    pname.org <- pname
    n <- 5
    mu <- pa[1]
    nu <- pa[2]
    a <- pa[3]
    sigma1 <- pa[4]
    sigma2 <- pa[5]

    z <- .Fortran(C_smplxa,
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
                  fn = double(ipmax),
                  mple = double(n*ipmax),
                  xx = double(n*itmax1),
                  std = double(itmax1),
                  f = double(itmax1),
                  itr = integer(1),
                  nip = integer(1),
                  ipr = integer(ipmax),
                  as.integer(ipflg))

    nip <- z$nip
    mples <- array(z$mple, dim = c(ipmax, n))
    mple <- mples[nip, 1:n]

  } else if (model == "TypeB") {

# (mu1, nu2, nu, sigma1, sigma2) for the random variable Poisson

    pname.org <- c("mu", "nu", "a", "sigma1", "sigma2") 
    n <- 5
    mu1 <- pa[1]
    mu2 <- pa[2]
    nu <- pa[3]
    sigma1 <- pa[4]
    sigma2 <- pa[5]

    z <- .Fortran(C_smplxb,
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
                  fn = double(ipmax),
                  mple = double(n*ipmax),
                  xx = double(n*itmax1),
                  std = double(itmax1),
                  f = double(itmax1),
                  itr = integer(1),
                  nip = integer(1),
                  ipr = integer(ipmax),
                  as.integer(ipflg))

    nip <- z$nip
    mples <- array(z$mple, dim = c(ipmax, n))
    pp <- mples[nip, 1:n]
    mple <- c(pp[1]*pp[3], pp[1]*(1-pp[3]), pp[2], pp[4], pp[5])

  } else if (model == "TypeC") {

# (mu1, mu2, nu1, nu2, sigma1, sigma2) for the random variable Poisson

    pname.org <- c("lambda", "nu1", "a", "sigma1", "sigma2")
    n <- 5
    mu1 <- pa[1]
    mu2 <- pa[2]
    nu1 <- pa[3]
    nu2 <- pa[4]
    sigma1 <- pa[5]
    sigma2 <- pa[6]

    z <- .Fortran(C_smplxc,
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
                  fn = double(ipmax),
                  mple = double(n*ipmax),
                  xx = double(n*itmax1),
                  std = double(itmax1),
                  f = double(itmax1),
                  itr = integer(1),
                  nip = integer(1),
                  ipr = integer(ipmax),
                  as.integer(ipflg))

    nip <- z$nip
    mples <- array(z$mple, dim = c(ipmax, n))
    pp <- mples[nip, 1:n]
    mple <- c(pp[1] * pp[3] / pp[2],
              pp[1] * (1 - pp[3]) * pp[4] / (pp[2] * pp[5]),
              pp[2], pp[2] * pp[5] / pp[4], pp[4], pp[5])

  }

  ipri <- z$ipr[1:nip]
  names(mple) <- pname
  mples <- mples[1:nip, 1:n]
  if (nip != 1)
    mples <- data.frame(mples)
  names(mples) <- pname.org

  it1 <- z$itr
  it2 <- 1
  if (ipflg == 2 || ipflg == 3)
    it2 <- it1
  xx <- array(z$xx, dim = c(n, itmax1))

  f <- z$f[1:it2]
  para <- xx[1:n, 1:it2]
  std <- z$std[1:it2]
  param <- data.frame(t(para))
  names(param) <- pname.org
  aic <- 2 * f[it2] + 2 * n

  out <- list(mple = mple, log.mpl = -f[it2], aic = aic, 
             process1 = list(cflg = ipri, logl = z$fn[1:nip], mples = mples),
             process2 = list(logl = f, stderr = std, pa.normal = param),
             input.val = list(model = model, points = xy.points, pa.init = pa,
                              eps = eps, uplimit = uplimit, skip = skip))

  class(out) <- "mple"
  out
}


set.pars <- function(model, xy.points) {
  switch(model,
   "IP" = pname <- c("mu", "nu", "p", "c"),
   "Thomas" =    pname <- c("mu", "nu", "sigma"),
   "TypeA" =   pname <- c("mu", "nu", "a", "sigma1", "sigma2"),
   "TypeB" =   pname <- c("mu1", "mu2", "nu", "sigma1", "sigma2"),
   "TypeC" =   pname <- c("mu1", "mu2", "nu1", "nu2", "sigma1", "sigma2")
  )

  n <- dim(xy.points)[1]
  np <- length(pname)
  pars <- rep(0, np)
  delta <- 0.001

  zs <- nonpara.palm(xy.points)
  p1 <- zs$peak1
  t <- zs$intsct
  p2 <- zs$peak2
  if (p1 == delta)
    p1 <- p1 + delta

  if (model == "Thomas") {
    pars[1] <- sqrt(n)
    pars[2] <- sqrt(n)
    pars[3] <- (p1 + t) / 3

  } else if (model == "IP") {
    pars[1] <- sqrt(n)
    pars[2] <- sqrt(n)
    pars[3] <- 1.1
    pars[4] <- p1

  } else if (model == "TypeA") {
    pars[1] <- sqrt(n)
    pars[2] <- sqrt(n)
    pars[3] <- 0.5
    pars[4] <- p1
    pars[5] <- p2

  } else if (model == "TypeB") {
    pars[1] <- sqrt(n) / 2
    pars[2] <- sqrt(n) / 2
    pars[3] <- sqrt(n)
    pars[4] <- p1
    pars[5] <- t

  } else if (model == "TypeC") {
      n2 <- sqrt(n / 2)
    pars[1] <- n2
    pars[2] <- n2
    pars[5] <- p1
    pars[6] <- (p1 + t) / 2
    pars[3] <- sqrt(2 * n) / (1 + pars[6] / pars[5])
    pars[4] <- sqrt(2 * n) - pars[3]
  }

  return(pars)
}


nonpara.palm <- function(xy.points, log = "xy") {

  x <- xy.points[, 1]
  y <- xy.points[, 2]
  np <- length(x)
  delta <- 0.001
  ty <- 1
  mu <- rep(0, 2); nu <- rep(0, 2); sigma <- rep(0, 2)
  m <- 0

  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  z <- .Fortran(C_palmt,
                as.double(x),
                as.double(y),
                as.integer(np),
                as.double(delta),
                as.double(ty),
                as.double(mu),
                as.double(nu),
                as.double(sigma),
                as.integer(m),
                as.integer(jmax),
                palm = double(jmax),
                palm1 = double(m*jmax))

  palm <- z$palm
  r <- rep(1:jmax) * delta

  n <- length(r)
  ip1 <- which.max(palm)

  i <- c(ip1+1:n)
  is1 <- min(which(palm[i] < 1))
  is1 <- ip1 + is1
  ip2 <- which.max(palm[(is1+1):n])
  ip2 <- is1 + ip2
  xx <- r[ip1:is1]
  yy <- palm[ip1:is1]
  nplm <- lm(yy ~xx)
  a <- coef(nplm)[2]
  b <- coef(nplm)[1]
  is1n <- as.integer(((1-b)/a)/delta)

  out <- list(peak1 = r[ip1], intsct = r[is1n], peak2 = r[ip2])
}


##### S3method #####

coef.mple <- function(object, ...) {
  object$mple
}

summary.mple <- function(object, ...) {
  model <-object$input.val$model
  switch(model,
    "Thomas" =  mtitle <- "Thomas model",
    "IP" = mtitle <- "Inverse-power type model",
    "TypeA" = mtitle <- "Type A model",
    "TypeB" = mtitle <- "Type B model",
    "TypeC" = mtitle <- "Type C model",
  )
  np <- dim(object$input.val$points)[1]
  cat(sprintf("%s\n", mtitle))
  cat(sprintf("The number of point pattern: %8i\n\n", np))

  pname <- names(object$mple)
  n <- length(pname)
  mples <- matrix(object$mple)
  mples <- round(mples, digits = 5)
  rownames(mples) <- names(object$mple)
  colnames(mples) <- "      MPLE"
  print(mples)

  log.mpl <- object$log.mpl
  cat(sprintf("\nLog(MPL): %15.3f\n", log.mpl))
  aic <- object$aic
  cat(sprintf("AIC:      %15.3f\n", aic))
}


print.mple <- function(x, print.level = 0, ...) {
  if (print.level < 0 || print.level > 3)
    stop("`print.level' is invalid.", call. = FALSE)

  if (print.level == 1 || print.level == 3) {
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
    for (j in 1:ipri2) {
      if (cflg[j] == -1)
        cat(sprintf(" testfn = %18.10e", logl[j]))
      if (cflg[j] == 1)
        cat(sprintf(" update = %18.10e", logl[j]))
      for (i in 1:nn)
        cat(sprintf("\t%12.4e", mples[j, i]))
      cat("\n")
    }
  }

  if (print.level > 1) {
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
  n <- length(pname)
  cat("\n\n Parameter")
  for (i in 1:n)
     cat(sprintf("\t%12s", pname[i]))
  cat("\n")
  cat(" Initial value")
  for (i in 1:n)
    cat(sprintf("\t%12.4f", x$input.val$pa.init[i]))
  cat("\n")
  cat(" MPLE    ")
  for (i in 1:n)
    cat(sprintf("\t%12.4f", x$mple[i]))
  cat("\n\n")

  invisible(x)
}

plot.mple <- function(x, ...) {

  model <- x$input.val$model
  para <- x$process2$pa.normal

  pname <- names(para)
  n <- dim(para)[2]
  n1 <- n + 1

  old.par <- par(no.readonly = TRUE)
  new.mai <- old.par$mai
  new.mai[2] <- new.mai[1] * 1.2
  par(mai = new.mai, xaxs = "i")

  switch(model,
    "Thomas" =  mtitle <- "Thomas model",
    "IP" = mtitle <- "Inverse-power type model",
    "TypeA" = mtitle <- "Type A model",
    "TypeB" = mtitle <- "Type B model",
    "TypeC" = mtitle <- "Type C model",
  )

  plot(para[, 1], type = "l", ylim = c(0.25, 1.75), main = mtitle,
       xlab = "Number of iterations",
       ylab = "Values of the normalized parameters",
       col = 2, cex.main = 2.0, cex.lab = 1.5, ...)

  for (i in 2:n) {
    par(new = TRUE)
    plot(para[, i], type = "l", ylim = c(0.25, 1.75), main = "",
         xlab = "", ylab = "", col = i + 1, ...)
  }
  legend("topright", legend = pname, col = c(2:n1), lty = 1, cex = 1.2)

  par(mai = old.par$mai, xaxs = old.par$xaxs, new = old.par$new)

  invisible(NULL)
}
