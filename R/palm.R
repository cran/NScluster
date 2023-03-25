palm.cppm <- function(mple, pars = NULL, delta = 0.001, uplimit = 0.3) {
  model <- mple$input.val$model
  xy <- mple$input.val$points
  mples <- mple$mple

  switch(model,
   "IP" = pname <- c("mu", "nu", "p", "c"),
   "Thomas" =    pname <- c("mu", "nu", "sigma"),
   "TypeA" =   pname <- c("mu", "nu", "a", "sigma1", "sigma2"),
   "TypeB" =   pname <- c("mu1", "mu2", "nu", "sigma1", "sigma2"),
   "TypeC" =   pname <- c("mu1", "mu2", "nu1", "nu2", "sigma1", "sigma2")
  )

  x <- xy[, 1]
  y <- xy[, 2]
  np <- length(x)

  m <- 0
  lname <- NULL

# ty : the variable for the standardized coordinates of points
# in the rectangular region
  ty <- 1
  r1 <- 0.5e0
  rmax <- r1 / delta
  jmax <- as.integer(rmax)

  pa <- ParsCheck(pars, pname)
  if (is.null(pars) == FALSE && is.null(pa) == TRUE) 
    warning("The true parameters (pars) is invalid.", call. = FALSE)

  if (model == "IP") {

    mu <- rep(0, 2); nu <- rep(0, 2); p <- rep(0, 2); c <- rep(0, 2)

    m <- 1
    mu[1] <- mples[1]
    nu[1] <- mples[2]
    p[1] <- mples[3]
    c[1] <- mples[4]
    lname <- c(lname, "MPLE")
    if (is.null(pa) == FALSE) {
      m <- 2
      mu[2] <- pa[1]
      nu[2] <- pa[2]
      p[2] <- pa[3]
      c[2] <- pa[4]
      lname <- c(lname, "true")
    }

    z <- .Fortran(C_xqgausip,
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
                  as.integer(jmax),
                  palm = double(jmax),
                  palm1 = double(m*jmax))

  } else if (model == "Thomas") {

    mu <- rep(0, 2); nu <- rep(0, 2); sigma <- rep(0, 2)

    m <- 1
    mu[1] <- mples[1]
    nu[1] <- mples[2]
    sigma[1] <- mples[3]
    lname <- c(lname, "MPLE")
    if (is.null(pa) == FALSE) {
      m <- 2
      mu[2] <- pa[1]
      nu[2] <- pa[2]
      sigma[2] <- pa[3]
      lname <- c(lname, "true")
    }

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

  } else if (model == "TypeA") {

    mu <- rep(0, 2); nu <- rep(0, 2); a <- rep(0,2)
    sigma1 <- rep(0, 2); sigma2 <- rep(0,2)

    m <- 1
    mu[1] <- mples[1]
    nu[1] <- mples[2]
    a[1] <- mples[3]
    sigma1[1] <- mples[4]
    sigma2[1] <- mples[5]
    lname <- c(lname, "MPLE")
    if (is.null(pa) == FALSE) {
      m <- 2
      mu[2] <- pa[1]
      nu[2] <- pa[2]
      a[2] <- pa[3]
      sigma1[2] <- pa[4]
      sigma2[2] <- pa[5]
      lname <- c(lname, "true")
    }

    z <- .Fortran(C_xqgausa,
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
                  as.integer(jmax),
                  palm = double(jmax),
                  palm1 = double(m*jmax))

  } else if (model == "TypeB") {

    mu <- rep(0, 2); nu <- rep(0, 2); a <- rep(0,2)
    sigma1 <- rep(0, 2); sigma2 <- rep(0,2)

    m <- 1
    mu[1] <- mples[1] + mples[2]
    nu[1] <- mples[3]
    a[1] <- mples[1] / mu[1]
    sigma1[1] <- mples[4]
    sigma2[1] <- mples[5]
    lname <- c(lname, "MPLE")
    if (is.null(pa) == FALSE) {
      m <- 2
      mu[2] <- pa[1] + pa[2]
      nu[2] <- pa[3]
      a[2] <- pa[1] / mu[2]
      sigma1[2] <- pa[4]
      sigma2[2] <- pa[5]
      lname <- c(lname, "true")
    }

    z <- .Fortran(C_palmb,
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
                  as.integer(jmax),
                  palm = double(jmax),
                  palm1 = double(m*jmax))

  } else if (model == "TypeC") {

    lamb <- rep(0, 2); nu1 <- rep(0, 2); a <- rep(0,2)
    sigma1 <- rep(0, 2); sigma2 <- rep(0,2)

    m <- 1
    lamb[1] <- mples[1] * mples[3] + mples[2] * mples[4]
    nu1[1] <- mples[3]
    a[1] <- (mples[1] * mples[3]) / lamb[1]
    sigma1[1] <- mples[5]
    sigma2[1] <- mples[6]
    lname <- c(lname, "MPLE")
    if (is.null(pa) == FALSE) {
      m <- 2
      lamb[2] <- pa[1] * pa[3] + pa[2] * pa[4]
      nu1[2] <- pa[3]
      a[2] <- (pa[1] * pa[3])/ lamb[2]
      sigma1[2] <- pa[5]
      sigma2[2] <- pa[6]
      lname <- c(lname, "true")
    }

    z <- .Fortran(C_palmc,
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
                  as.integer(jmax),
                  palm = double(jmax),
                  palm1 = double(m*jmax))

  }

  palm <- z$palm
  r <- rep(1:jmax) * delta
  palm1 <- array(z$palm1, dim = c(jmax, m))
  colnames(palm1) <- lname

  out <- list(r = r, np.palm = palm, norm.palm = palm1)
  out <- c(list(model = model), out)
  class(out) <- c("Palm")
  out
}




##### S3method #####

print.Palm <- function(x, ...) {
  lname <- colnames(x$norm.palm)
  cat("\n\t\t Palm intensity function\n")
  cat("\tr\t non-parametric\t     normalized Palm\n\n")
  n <- length(x$r)
  m <- dim(x$norm.palm)[2]

  if (m != 0) {
      cat("\t\t\t")
      for (j in 1:m)
      cat(sprintf("\t\t%s", lname[j]))
  }

  maxp <- which.max(x$np.palm)
  is <- 0
  for (i in 1:n) {
    cat(sprintf("\n %10.3f\t%12.5f", x$r[i], x$np.palm[i]))
    if (i > maxp)
      if (is == 0 && x$np.palm[i] < 1)
        is <- i
    for (j in 1:m)
      cat(sprintf("\t%12.5f", x$norm.palm[i, j]))
  }
  cat("\n\n")
  invisible(x)
}

plot.Palm <- function(x, ..., log = "xy") {

  model1 <- x$model
  r <- x$r
  palm <- x$np.palm
  palm1 <- x$norm.palm
  lname1 <- model1

  m1 <-  dim(palm1)[2]
  ymax <- max(palm, palm1) * 2

  argh <- list(...)
  n <- length(argh)
  nobj <- 0
  model2 <- NULL
  palm2 <- list()
  lname2 <- NULL


  if (n != 0)
    for (i in 1:n) {
      obj <- argh[[i]]
      if (is(obj) != "Palm") {
        warning("Additional object is invalid 'class' and was ignored.",
                call. = FALSE )
      } else if (all.equal(r, obj$r) != TRUE) {
        warning("Additional object is invalid 'r' and was ignored.",
                call. = FALSE )
      } else if (all.equal(palm, obj$np.palm) != TRUE) {
        warning("Additional object is invalid 'np.palm' and was ignored.",
                call. = FALSE )
      } else if (all.equal(palm1[, 1], obj$norm.palm[, 1]) == TRUE) {
        warning("Additional object is duplicated and was ignored",
                call. = FALSE )
      } else {
        if (nobj < 4) {
          model2 <- obj$model
          palm2 <- c(palm2, list(obj$norm.palm))
          lname2 <- c(lname2, model2)
          nobj <- nobj + 1
          ymax <- max(ymax, obj$norm.palm * 2)
        }
      }
   }

  plot(x = r, y = palm, pch = 20, ylim = c(0.5, ymax), log = log,
       main = "", xlab = "r", ylab = "Palm intensity", cex.lab = 1.5)
  abline(h = 1)

  lcol <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  leg.txt <- "non-parametric"
  for (j in 1:m1) {
    par(new = TRUE)
    plot(x = r, y = palm1[, j], type = "l", ylim = c(0.5, ymax),
         col = lcol[j + 1], log = log, xlab = "", ylab = "",
         xaxt = "n", yaxt = "n")
  }

  mm <- m1
  mm1 <- mm + 1
  if (nobj > 0) {
    for (i in 1:nobj) {
      m2 <-  dim(palm2[[i]])[2]
      par(new = TRUE)
      mm <- mm + 1
      mm1 <- mm1 + 1
      if (m2 == 1) {
        plot(x = r, y = palm2[[i]], type = "l", ylim = c(0.5, ymax),
             col = lcol[mm1], log = log, xlab = "", ylab = "",
             xaxt = "n", yaxt = "n")
      } else {
        plot(x = r, y = palm2[[i]][, 1], type = "l", ylim = c(0.5, ymax),
             col = lcol[mm1], log = log, xlab = "", ylab = "",
             xaxt = "n", yaxt = "n")
      }
    }
  }

  if (m1 == 2) {
    leg.txt <- c(leg.txt, "true")
    lcol <- c(1, 3, 2, 4, 5, 6, 7, 8, 9)
  }
  leg.txt <- c(leg.txt, lname1, lname2)

  legend("topright", legend = leg.txt[1:mm1], col = lcol[1:mm1],
         pch = c(20, rep(NA, mm)), lty = c(NA, rep(1, mm)), cex = 1.1)
}


