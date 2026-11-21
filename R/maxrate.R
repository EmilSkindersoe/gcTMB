#' Maximum growth rate and its timing
#'
#' Estimates the time at which growth curves grow fastest, the inflection point, the maximum growth
#' rate there, and the fraction of the asymptote reached at that time, with standard errors by the
#' delta method. Confidence intervals by \code{confint}: Wald, profile likelihood or parametric
#' bootstrap intervals.
#'
#' @param x A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param curves Numeric matrix, list of matrices, or NULL. Linear functions to identify the parameters of the growth curves, with 5 rows (b1-b5) and one column for each model coefficient, see \code{\link{timeto}}; default is the identity matrix, for models without covariates in fform.
#' @param newdata Data frame or NULL. The values of the covariates in fform, one row for each growth curve, as an alternative to curves.
#'
#' @details The curves are the population curves, with all random effects at zero; the standard
#' errors include the uncertainty of the fixed effects. For each curve, \code{maxrate} gives
#' \itemize{
#'   \item \code{tmax}: the time of the maximum growth rate, the inflection point of the curve;
#'   \item \code{rate}: the maximum growth rate, the slope of the curve at \code{tmax}, in units of the
#'     response per unit of time;
#'   \item \code{fraction}: the fraction of the asymptote b1 reached at \code{tmax}, how much of the
#'     growth has occurred when the growth is fastest.
#' }
#' For the logistic (Richards), hill and symmetric weibull and gompertz curves these have closed
#' forms; for the gompertz and weibull curves with an asymmetry b5 different from 1 the inflection
#' point is the root of an equation, see the vignette. The hill and weibull curves have an inflection
#' point at a positive time only for b3 b5 > 1; otherwise they grow fastest at t = 0, with an
#' unbounded rate for the hill curve, and the results are NA. The linear model grows at the constant
#' rate b1, without an inflection point: its \code{tmax} and \code{fraction} are NA.
#'
#' A truncated curve (the default of \code{\link{gcTMB}}, at b4 = 0) is the curve without truncation,
#' f(t), shifted down to zero at its threshold b4 and rescaled, b1 (f(t) - f(b4))/(b1 - f(b4)) beyond b4. It has its maximum growth rate at the
#' inflection point of the curve without truncation if that is after the threshold, and otherwise
#' just after the threshold, where the slope jumps from zero: then \code{tmax} is b4 and
#' \code{fraction} zero.
#'
#' The standard errors are by the delta method, with numerical derivatives of the quantities with
#' respect to the parameters of the curves. Profile intervals are available for \code{tmax} (by b2,
#' unless the curve is at its threshold) and \code{rate} (by b1); for \code{fraction}, which depends
#' on the shape of the curve only, \code{confint} returns Wald intervals with a message for
#' \code{method="profile"}.
#'
#' @return A maxrateTMB object: a list with the estimates (\code{coef}), their covariance matrix
#'   (\code{vcov}), and what \code{confint} and \code{profile} need. The estimates are named
#'   \code{tmax}, \code{rate} and \code{fraction}, followed by the name of the curve for curves from
#'   \code{newdata} or named \code{curves}, e.g. \code{rate:trt=0.1}; several unnamed \code{curves}
#'   are numbered.
#' @seealso \code{\link{confint.maxrateTMB}} for confidence intervals, \code{\link{timeto}}
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25, model="logistic")
#' mr <- maxrate(fit)
#' mr
#' confint(mr, method="wald")
#'
#' # the Richards curve, with an estimated asymmetry
#' fitr <- gcTMB(height ~ DAS, data=FabaHeights25, fix=c(b5=FALSE), link=c(b5="log"))
#' maxrate(fitr)
#'
#' # a curve for each treatment level
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fit2 <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25)
#' maxrate(fit2, newdata=data.frame(trt=levels(FabaHeights25$trt)))
#' @export
maxrate <- function(x, curves=NULL, newdata=NULL){
  if (!inherits(x, "gcTMB")) stop("x must be a gcTMB object.")
  if (!is.null(newdata)){
    if (!is.null(curves)) stop("give either curves or newdata, not both.")
    # population curves: not named by the grouping variables of the random effects
    curves <- newdatacurves(x, newdata, groups=FALSE)
  }
  if (is.null(curves)) curves <- diag(5)
  if (!is.list(curves)) curves <- list(curves)
  # unnamed curves: no name for a single curve, numbered otherwise
  if (is.null(names(curves))) names(curves) <- if (length(curves) == 1) "" else seq_along(curves)
  curves <- checkcurves(x, curves, is.null(newdata))

  pv <- x$start
  nb <- sum(!x$fix)
  pv[!x$fix] <- x$sdrl$par.fixed[seq_len(nb)]
  qn <- c("tmax", "rate", "fraction")
  nq <- length(qn)
  nc <- length(curves)
  est <- numeric(nq*nc)
  grad <- matrix(0, nq*nc, nb)
  for (i in seq_len(nc)){
    L <- curves[[i]]
    lpv <- as.vector(L %*% pv)
    id <- (i - 1)*nq + seq_len(nq)
    est[id] <- maxratevalue(x, lpv)
    # derivatives w.r.t. the linear predictors of b1-b5, by central differences
    g <- sapply(1:5, function(j){
      if (all(L[j, !x$fix] == 0)) return(rep(0, nq))
      h <- 1e-5 * max(1, abs(lpv[j]))
      e <- replace(numeric(5), j, h)
      (maxratevalue(x, lpv + e) - maxratevalue(x, lpv - e)) / (2*h)
    })
    grad[id, ] <- (matrix(g, nrow=nq) %*% L)[, !x$fix, drop=FALSE]
  }
  vc <- fixedcov(x)[seq_len(nb), seq_len(nb), drop=FALSE]
  v <- grad %*% vc %*% t(grad)
  nm <- curvelabels(rep(qn, nc), rep(names(curves), each=nq))
  if (anyNA(est) && !islinear(modelindex(x$model))){
    nacurves <- unique(names(curves)[ceiling(which(is.na(est))/nq)])
    warning(if (identical(nacurves, "")) "the curve has" else paste("the curves", paste(nacurves, collapse=", "), "have"),
            " no inflection point at a positive time (b3 b5 <= 1); the results are NA.", call.=FALSE)
  }
  names(est) <- nm
  dimnames(v) <- list(nm, nm)
  res <- list(estimates=est, variance=v, fit=x, curves=curves,
              quantity=setNames(rep(qn, nc), nm), curve=setNames(rep(seq_len(nc), each=nq), nm))
  class(res) <- "maxrateTMB"
  res
}

# names of estimates: the label alone for an unnamed curve, otherwise label:curve
curvelabels <- function(label, curve){
  ifelse(curve == "", label, paste(label, curve, sep=":"))
}

# curves as a list of matrices with 5 rows and one column for each coefficient of the fit; a matrix
# without a row for the asymmetry b5 keeps b5 of the fit, when it is a single coefficient
checkcurves <- function(x, curves, hint=TRUE){
  ncf <- length(x$start)
  npi <- rep(1:5, vapply(x$Xs, ncol, integer(1)))
  lapply(curves, function(L){
    if (!is.matrix(L) || nrow(L) != 5 || ncol(L) != ncf){
      stop("each matrix in curves must have 5 rows (b1-b5) and one column for each of the ", ncf,
           " model coefficients: ", paste(x$coefnames, collapse=", "), ".",
           if (hint && ncf > 5) " Alternatively, give the covariate values of the curves in newdata.")
    }
    if (sum(npi == 5) == 1 && all(L[5, ] == 0)) L[5, npi == 5] <- 1
    L
  })
}

# tmax, rate and fraction of a curve of the fit with the parameters lpv (b1-b5 on the scale of the
# linear predictors)
maxratevalue <- function(fit, lpv){
  b <- curvepars(fit, lpv)
  maxrateR(modelindex(fit$model), b, fittrunc(fit))
}

# inflection point z of the curve without truncation for b1 = 1, as a function of b2, b3 and the
# asymmetry a = b5: the time t, the fraction G of the asymptote and the rate r there; NA without an
# inflection point at a positive time
maxrate0 <- function(mod, b2, b3, a){
  na <- c(t=NA, G=NA, r=NA)
  if (!all(is.finite(c(b2, b3, a))) || a <= 0) return(na)
  switch(mod,
    { # logistic, Richards: at z = 0
      c(t=b2, G=(1 + a)^(-1/a), r=(1 + a)^(-(1 + a)/a) / b3)
    },
    { # gompertz with a power a of the survival function: w = exp(-z), root of
      # (w - 1)(1 - exp(-w)) = (a - 1) w exp(-w)
      w <- if (a == 1) 1 else rootfind(function(w) (w - 1)*(-expm1(-w)) - (a - 1)*w*exp(-w), c(1e-8, 2), "upX")
      if (is.na(w)) return(na)
      S <- -expm1(-w)
      c(t=b2 - b3*log(w), G=-expm1(a*log(S)), r=a*S^(a - 1)*w*exp(-w) / b3)
    },
    { # hill, LL.5: v = (b2/t)^b3 = (b3 + 1)/(b3 a - 1)
      if (b3*a <= 1 || b2 <= 0) return(na)
      v <- (b3 + 1)/(b3*a - 1)
      tt <- b2 * v^(-1/b3)
      c(t=tt, G=(1 + v)^(-a), r=b3*a*v*(1 + v)^(-a - 1) / tt)
    },
    { # exponentiated weibull: s = (b3 - 1)/b3 (t/b2)^b3, root of (a - 1) s/(e^s - 1) - s + 1 - 1/b3 = 0
      cc <- (b3 - 1)/b3
      if (b3 <= 1 || b3*a <= 1 || b2 <= 0) return(na)
      s <- if (a == 1) cc else rootfind(function(s) (a - 1)*s/expm1(s) - s + 1 - 1/b3, c(1e-10, 2), "downX")
      if (is.na(s)) return(na)
      tt <- b2 * (s/cc)^(1/b3)
      E <- -expm1(-s)
      c(t=tt, G=E^a, r=b3*a*s*exp(-s)*E^(a - 1) / tt)
    },
    { # normal, a power a of F: z root of (a - 1) phi(z)/Phi(z) - z = 0, z = 0 for a = 1
      z <- if (a == 1) 0 else rootfind(function(z) (a - 1)*exp(dnorm(z, log=TRUE) - pnorm(z, log.p=TRUE)) - z, c(-1, 1), "downX")
      if (is.na(z)) return(na)
      c(t=b2 + b3*z, G=exp(a*pnorm(z, log.p=TRUE)), r=a*exp((a - 1)*pnorm(z, log.p=TRUE))*dnorm(z) / b3)
    },
    { # gamma with shape k = b3 and scale b2/(k - 1), a power a of F: x = t (k - 1)/b2, root of
      # (a - 1) p(x)/P(x) + (k - 1)/x - 1 = 0, x = k - 1 for a = 1
      k <- b3
      if (k <= 1 || k*a <= 1 || b2 <= 0) return(na)
      x <- if (a == 1) k - 1 else rootfind(function(x) (a - 1)*exp(dgamma(x, k, log=TRUE) - pgamma(x, k, log.p=TRUE)) + (k - 1)/x - 1,
                                           c(1e-8, k), "downX")
      if (is.na(x)) return(na)
      c(t=x*b2/(k - 1), G=exp(a*pgamma(x, k, log.p=TRUE)), r=a*exp((a - 1)*pgamma(x, k, log.p=TRUE))*dgamma(x, k)*(k - 1)/b2)
    })
}

# root of a monotone function in an interval, extended in the direction ext
rootfind <- function(f, int, ext){
  tryCatch(uniroot(f, int, extendInt=ext, tol=1e-15, maxiter=1000)$root, error=function(e) NA)
}

# slope dG/dt of the curve without truncation (b1 = 1) at time t
dGdt <- function(mod, t, b2, b3, a){
  m <- gcmodel(mod)
  z <- m$z(t, b2, b3)
  # dz/dt: the gamma time axis is linear, the other models for positive times are on the log time axis
  dz <- if (isTRUE(m$shapeF)) (b3 - 1)/b2 else if (m$logtime) b3/t else 1/b3
  if (m$powS){
    g <- a * exp((a - 1)*m$logS(z, b3)) * m$dens(z, b3)
  } else {
    zs <- z - m$shift(a)
    g <- exp(logG(m, z, a, b3)) * m$pow(a) * m$dens(zs, b3) / m$F(zs, b3)
  }
  g * dz
}

# tmax, rate and fraction of a curve with the parameters b (b1-b5 on their original scale); truncated
# at b4 (trunc=TRUE): the curve without truncation scaled by 1/S4 beyond b4, with S4 = 1 - G(z(b4))
maxrateR <- function(mod, b, trunc=FALSE){
  b <- as.numeric(b)
  # the linear model grows at the constant rate b1, without an inflection point
  if (islinear(mod)) return(c(tmax=NA, rate=b[1], fraction=NA))
  # chapmanrichards, beta and wykoff: their own closed forms, never truncated
  if (isstart(mod)){
    m0 <- gcmodel(mod)$mr(b[2], b[3], b[4], b[5], b1=b[1])
    return(c(tmax=m0[["t"]], rate=m0[["r"]], fraction=m0[["G"]]))
  }
  m0 <- maxrate0(mod, b[2], b[3], b[5])
  out <- c(tmax=m0[["t"]], rate=b[1]*m0[["r"]], fraction=m0[["G"]])
  if (!trunc) return(out)
  m <- gcmodel(mod)
  # a threshold up to zero does not truncate the curves on the log time axis
  if (m$logtime && b[4] <= 0) return(out)
  S4 <- exp(logS4(m, b[2], b[3], b[4], b[5]))
  if (!is.na(m0[["t"]]) && m0[["t"]] > b[4]){
    return(c(tmax=m0[["t"]], rate=b[1]*m0[["r"]]/S4, fraction=1 - (1 - m0[["G"]])/S4))
  }
  # the curve grows fastest just after its threshold
  c(tmax=b[4], rate=b[1]*dGdt(mod, b[4], b[2], b[3], b[5])/S4, fraction=0)
}


#' maxrateTMB Methods
#'
#' @param x,object A maxrateTMB object, the result of \code{\link{maxrate}}.
#' @param ... Not used.
#'
#' @return The estimates of the timing and value of the maximum growth rate, and their covariance
#'   matrix
#'
#' @name maxrateTMB-methods
#' @aliases print.maxrateTMB coef.maxrateTMB vcov.maxrateTMB
NULL

#' @rdname maxrateTMB-methods
#' @method print maxrateTMB
#' @export
print.maxrateTMB <- function(x, ...){
  cat("Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):\n\n")
  print(data.frame(Estimate=x$estimates, Std.Error=sqrt(diag(x$variance))), digits=3)
  invisible(x)
}

#' @rdname maxrateTMB-methods
#' @export
coef.maxrateTMB <- function(object, ...) object$estimates

#' @rdname maxrateTMB-methods
#' @export
vcov.maxrateTMB <- function(object, ...) object$variance


#' Confidence intervals for the maximum growth rate and its timing
#'
#' Wald, profile likelihood or parametric bootstrap intervals for the estimates of
#' \code{\link{maxrate}}.
#'
#' @param object A maxrateTMB object, the result of \code{\link{maxrate}}.
#' @param parm Character or integer vector, or NULL. Estimates, names or positions in \code{coef(object)}; default all
#' @param level Numeric between 0 and 1. The confidence level
#' @param method Character string. "profile" for profile likelihood intervals (default), "wald" for Wald intervals, or "boot" for parametric bootstrap percentile intervals
#' @param nsim Integer. Number of bootstrap simulations
#' @param seed Integer or NULL. An optional seed for the random number generator of the bootstrap
#' @param boot A bootstrap.gcTMB object or NULL. A bootstrap sample of the model from \code{\link{bootstrap}}, used instead of new simulations
#' @param ... Further arguments passed to \code{\link{profile.maxrateTMB}} or \code{\link{bootstrap}}.
#'
#' @details Profile intervals are available for \code{tmax} and \code{rate}, see \code{\link{maxrate}};
#' for \code{fraction}, and for \code{tmax} of a curve at its threshold, Wald intervals are returned
#' with a message. Wald intervals use t quantiles with small-sample inference
#' (\code{\link{gcTMBcontrol}}). Bootstrap fits whose curve has no inflection point are left out of
#' the bootstrap intervals.
#'
#' @return a matrix with the lower and upper confidence limits
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25)
#' mr <- maxrate(fit)
#' confint(mr, method="wald")
#' confint(mr, parm=c("tmax", "rate"))
#' @export
confint.maxrateTMB <- function(object, parm=NULL, level=0.95, method=c("profile", "wald", "boot"), nsim=500, seed=NULL, boot=NULL, ...){
  method <- if (!is.null(boot)) "boot" else match.arg(method)
  est <- coef(object)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown estimates in parm.")
  fit <- object$fit
  lab <- paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%")
  ci <- matrix(NA, length(parm), 2, dimnames=list(names(est)[parm], lab))
  se <- sqrt(diag(vcov(object)))[parm]
  z <- zquant(fit, level)
  wald <- cbind(est[parm] - z*se, est[parm] + z*se)
  if (method == "wald") { ci[] <- wald; return(ci) }
  if (method == "boot"){
    draws <- bootdraws(fit, getboot(fit, boot, nsim, seed, ...), FALSE)
    nb <- sum(!fit$fix)
    vals <- sapply(parm, function(i){
      L <- object$curves[[object$curve[[i]]]]
      q <- match(object$quantity[[i]], c("tmax", "rate", "fraction"))
      apply(draws, 2, function(par){
        pv <- fit$start
        pv[!fit$fix] <- par[seq_len(nb)]
        maxratevalue(fit, as.vector(L %*% pv))[q]
      })
    })
    vals <- matrix(vals, ncol=length(parm))
    ci[] <- t(apply(vals, 2, quantile, probs=c((1 - level)/2, (1 + level)/2), na.rm=TRUE, names=FALSE))
    return(ci)
  }
  # profile: tmax and rate, Wald intervals for the others
  prof <- profileable(object, parm)
  if (any(!prof)){
    message("Wald intervals for ", paste(names(est)[parm][!prof], collapse=", "),
            ": profile intervals are available for tmax and rate only, and not for tmax at the threshold of a truncated curve.")
    ci[!prof, ] <- wald[!prof, , drop=FALSE]
  }
  if (any(prof)){
    message("Computing profile confidence intervals ...")
    ci[prof, ] <- confint(profile(object, parm=parm[prof], level=level, ...), level=level)
  }
  ci
}

# estimates parm of a maxrateTMB object with a profile: tmax (not at the threshold of a truncated
# curve) and rate, with finite estimates
profileable <- function(object, parm){
  fit <- object$fit
  b4 <- vapply(parm, function(i){
    pv <- fit$start
    pv[!fit$fix] <- fit$sdrl$par.fixed[seq_len(sum(!fit$fix))]
    linkpars(linkindex(fit$link), as.vector(object$curves[[object$curve[[i]]]] %*% pv))[4]
  }, numeric(1))
  est <- coef(object)[parm]
  q <- object$quantity[parm]
  is.finite(est) & (q == "rate" | (q == "tmax" & !(fittrunc(fit) & abs(est - b4) < 1e-10)))
}

#' Profile likelihood of the maximum growth rate and its timing
#'
#' Profiles the negative log-likelihood over the time of the maximum growth rate (by b2) or the
#' maximum growth rate (by b1) of a curve, re-optimising the other coefficients.
#'
#' @param fitted A maxrateTMB object, the result of \code{\link{maxrate}}.
#' @param parm Character or integer vector, or NULL. Estimates to profile, names or positions in \code{coef(fitted)}; tmax and rate only. Default all of these.
#' @param level Numeric between 0 and 1. The highest confidence level the profiles are computed for.
#' @param maxit Integer. Maximum number of steps in each direction
#' @param progress Logical. Show a progress bar over the profiled estimates
#' @param ... Not used.
#'
#' @return a profile.maxrateTMB object, which has the confint and plot methods of
#'   \code{\link{profile.gcTMB}}
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25)
#' pr <- profile(maxrate(fit), parm="rate")
#' confint(pr)
#' @export
profile.maxrateTMB <- function(fitted, parm=NULL, level=0.99, maxit=100, progress=interactive(), ...){
  fit <- fitted$fit
  est <- coef(fitted)
  if (is.null(parm)) parm <- which(fitted$quantity != "fraction")
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown estimates in parm.")
  ok <- profileable(fitted, parm)
  if (any(!ok)) stop("profiles are only available for tmax and rate with finite estimates, and not for tmax at the threshold of a truncated curve: ",
                     paste(names(est)[parm][!ok], collapse=", "), ".")
  mod <- modelindex(fit$model)
  m <- gcmodel(mod)
  trunc <- fittrunc(fit)
  lnk <- linkindex(fit$link)
  pb <- progressbar(length(parm), progress)
  prof <- lapply(parm, function(i){
    q <- fitted$quantity[[i]]
    # the parameter of the curve computed from the profiled value: b2 for tmax, b1 for the rate
    ps <- if (q == "tmax") 2 else 1
    solve <- function(v, lpv){
      b <- linkpars(lnk, lpv)
      if (q == "tmax"){
        # tmax is b2 plus a shift for the curves on the time axis, and b2 times a factor on the log
        # time axis, given b3 and b5; beyond the threshold of a truncated curve
        if (trunc && v <= b[4]) return(NA)
        # chapmanrichards, beta and wykoff have their inflection point at b2, and the peak of
        # preecebaines is at b2 plus a shift that depends on b4 and b5 only
        if (isstart(mod)) return(suppressWarnings(linkfun(lnk[ps], v - (m$mr(b[2], b[3], b[4], b[5], b1=b[1])[["t"]] - b[2]))))
        u <- maxrate0(mod, if (m$logtime) 1 else 0, b[3], b[5])[["t"]]
        p <- if (m$logtime) v/u else v - u
      } else {
        # the rate is proportional to b1, or b1 itself for b1type="rate"; for preecebaines, to b1 - b3
        p <- if (fitrate(fit)) v else if (!is.null(m$mean)) b[3] + v / maxrateR(mod, replace(b, 1, b[3] + 1), trunc)[["rate"]] else
          v / maxrateR(mod, replace(b, 1, 1), trunc)[["rate"]]
      }
      suppressWarnings(linkfun(lnk[ps], p))
    }
    d <- profilesolved(fit, fitted$curves[[fitted$curve[[i]]]], ps, solve, est[[i]], sqrt(fitted$variance[i, i]),
                       names(est)[i], level, maxit)
    pb$tick()
    d
  })
  pb$close()
  names(prof) <- names(est)[parm]
  attr(prof, "estimates") <- est[parm]
  attr(prof, "cutoff") <- cutoffattr(fit, rep(TRUE, length(parm)))
  class(prof) <- c("profile.maxrateTMB", "profile.gcTMB")
  prof
}

# profile of a quantity of a curve (linear functions L of the coefficients) with the estimate t0 and
# standard error se: at each value, the coefficient of parameter ps with the largest weight in L is
# computed from the value and the other parameters by solvefun(value, lpv), on the scale of its link,
# and the other coefficients are re-optimised
profilesolved <- function(fit, L, ps, solvefun, t0, se, name, level, maxit){
  obj <- fit$obj
  # keep the state of the TMB object, as tmbprofile does
  keep <- c("last.par.best", "random.start", "value.best", "last.par", "inner.control", "tracemgc")
  old <- mget(keep, envir=obj$env, ifnotfound=list(NULL))
  on.exit(for (v in keep) assign(v, old[[v]], envir=obj$env))
  obj$env$inner.control$trace <- FALSE
  obj$env$tracemgc <- FALSE

  free <- which(!fit$fix)
  nf <- length(free)
  par0 <- fit$opt$par
  nllmin <- fit$opt$objective
  cutoff <- profcutoff(level, fit)
  Hfull <- tryCatch(solve(fit$sdrl$cov.fixed), error=function(e) NULL)
  lo <- if (is.null(fit$lower)) rep(-Inf, length(par0)) else fit$lower
  up <- if (is.null(fit$upper)) rep(Inf, length(par0)) else fit$upper
  oth <- setdiff(1:5, ps)

  w <- L[ps, free]
  if (all(w == 0)) stop(name, " does not depend on an estimated coefficient of b", ps, ".")
  kp <- which.max(abs(w))
  k <- free[kp]
  if (any(L[oth, k] != 0)) stop("the coefficient ", fit$coefnames[k], " of b", ps, " must not enter other parameters in curves.")

  pvfull <- function(par){
    pv <- fit$start
    pv[free] <- par[seq_len(nf)]
    pv
  }
  expand <- function(x, tt){
    par <- append(x, 0, after=kp - 1)
    pv <- pvfull(par)
    lpv <- as.vector(L %*% pv)
    p <- solvefun(tt, lpv)
    par[kp] <- (p - sum(L[ps, -k] * pv[-k])) / L[ps, k]
    list(par=par, lpv=lpv)
  }
  jacobian <- function(e, tt){
    h <- 1e-6
    d <- sapply(oth, function(j){
      hj <- replace(numeric(5), j, h)
      (solvefun(tt, e$lpv + hj) - solvefun(tt, e$lpv - hj)) / (2*h)
    })
    J <- diag(length(par0))[, -kp, drop=FALSE]
    J[kp, seq_len(nf - 1)] <- ((colSums(d * L[oth, free, drop=FALSE]) - L[ps, free]) / L[ps, k])[-kp]
    J
  }
  fn <- function(x, tt){
    e <- expand(x, tt)
    if (!isTRUE(e$par[kp] >= lo[kp] && e$par[kp] <= up[kp])) return(Inf)
    obj$fn(e$par)
  }
  gr <- function(x, tt){
    e <- expand(x, tt)
    as.vector(obj$gr(e$par) %*% jacobian(e, tt))
  }
  he <- if (is.null(Hfull)) NULL else function(x, tt){
    J <- jacobian(expand(x, tt), tt)
    t(J) %*% Hfull %*% J
  }
  step <- se / 5
  side <- function(dir){
    x <- par0[-kp]
    out <- NULL
    for (it in seq_len(maxit)){
      tt <- t0 + dir*it*step
      o <- tryCatch(suppressWarnings(nlminb(x, fn, gr, he, tt=tt, lower=lo[-kp], upper=up[-kp])), error=function(e) NULL)
      if (is.null(o) || !is.finite(o$objective)) break
      x <- o$par
      out <- rbind(out, data.frame(par=tt, deviance=2*(o$objective - nllmin)))
      if (out$deviance[nrow(out)] > cutoff) break
    }
    out
  }
  d <- rbind(side(-1), data.frame(par=t0, deviance=0), side(1))
  d[order(d$par), ]
}
