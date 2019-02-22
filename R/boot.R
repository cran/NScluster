boot.mple <- function(mple.out, n = 100, conf.level = 0.95, se = TRUE) {

  prog <- TRUE
  if (prog){
    pb <- txtProgressBar(min = 0, max = n, style = 3)
  }

  model <- mple.out$input.val$model
  pars <- mple.out$mple
  npar <- length(pars)
  boot.mples <- matrix(0, nrow = n, ncol = npar)
  colnames(boot.mples) <- names(pars)

  for (i in 1:n) {
    xy.points <- sim.cppm(model, pars, seed = NULL)$offspring$xy
    np <- dim(xy.points)[1]
    init.pars <- set.pars(model, xy.points)
    names(init.pars) <- names(pars)
    result <- mple.cppm(model, xy.points, init.pars, mple.out$input.val$eps,
                        mple.out$input.val$uplimit, mple.out$input.val$skip)
    boot.mples[i, ] <- result$mple
    if (prog){
      setTxtProgressBar(pb, i)
    }
  }

  prb1 <- (1 - conf.level) / 2
  prb2 <- (1 + conf.level) / 2
  interval <- t(apply(boot.mples, 2, quantile, probs = c(prb1, prb2),
                na.rm = TRUE))

  conf <- matrix(0, nrow = npar, ncol = 2)
  rownames(conf) <- names(pars)
  colnames(conf) <- c(paste(100 * prb1, "%"), paste(100 * prb2, "%"))
  conf <- interval

  if (se == TRUE) {
    n <- dim(boot.mples)[1]
    sd.mples <- apply(boot.mples, 2, sd)
    std.err <- sd.mples / sqrt(n)
    conf <- cbind(conf, std.err)
  }

  out <- list(mple = pars, boot.mples = boot.mples, confint = conf)
  class(out) <- "boot.mple"
  out
}


summary.boot.mple <- function(object, ...) {

  pars <- object$mple
  conf <- cbind(pars, object$confint)
  colnames(conf)[1] <- "MPLE"

  coef(list(coefficients = conf))

}
