#' Likelihood profiles of a gcTMB model
#'
#' Profiles the (marginal) log-likelihood for single parameters by fixing the parameter at a sequence
#' of values and re-optimising all other parameters, using \code{\link[TMB]{tmbprofile}}. With random
#' effects, the profiled likelihood is the Laplace approximation of the marginal likelihood.
#'
#' The profiles are computed on the scale of \code{coef}: log standard deviations of the random effects
#' (log_sd), parameters of their correlation matrix (corr_par), and the log residual standard deviation
#' (log_sigma), log precision (log_phi), or log size (log_theta) or dispersion (log_phi) of the negative
#' binomial families. The confint method transforms the intervals of the log
#' standard deviations and log precision to the natural scale; correlation parameters are only
#' transformed for terms with two random effects, as their parameter determines the correlation
#' only then. If the profile of a log standard deviation does not reach the cutoff at low values,
#' the likelihood is flat towards zero and the lower limit of the standard deviation is 0.
#' Bounds given by \code{lower} and \code{upper} in \code{\link{gcTMB}} are used for the profiled
#' coefficient, but not when re-optimising the other coefficients.
#'
#' The cutoff of the deviance is the quantile of the chi-square distribution with one degree of
#' freedom. For gaussian and lognormal models without random effects and small-sample inference (see
#' \code{\link{gcTMB}}), the cutoff of the coefficients b1-b5 is F-based instead,
#' sum(w) log(1 + F/(n - p)) with the quantile F of the F distribution with 1 and n - p degrees of
#' freedom, which gives the intervals of a t-test of each coefficient, as profile intervals of nls.
#'
#' @param fitted A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param parm Character or integer vector, or NULL. Parameters to profile, names or positions in \code{coef(fitted)}. Default are the
#'   estimated coefficients b1-b5, without the variance parameters.
#' @param level Numeric between 0 and 1. The highest confidence level the profiles are computed for: each profile extends until
#'   the deviance from the maximum exceeds the cutoff (see Details). For confint and plot, the confidence level of the intervals and of the cutoff line.
#' @param trace Logical. Print the progress of \code{\link[TMB]{tmbprofile}}
#' @param progress Logical. Show a progress bar over the profiled parameters
#' @param object,x A profile.gcTMB object, the result of \code{profile}.
#' @param ... Further arguments passed to \code{\link[TMB]{tmbprofile}}, e.g. \code{ystep}; ignored by
#'   the confint and plot methods
#'
#' @return profile returns a profile.gcTMB object, a list with a data frame for each parameter with
#'   the parameter values (\code{par}) and the deviance from the maximum (\code{deviance}), i.e. twice
#'   the difference of the negative log-likelihoods. confint returns a matrix of confidence intervals.
#'   plot returns a ggplot object (package ggplot2) with a panel for each parameter, showing the
#'   deviance, the cutoff for \code{level} (dashed), and the
#'   estimate (dotted).
#' @export
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25)
#' pr <- profile(fit, progress=FALSE)
#' confint(pr)
#' plot(pr)
#' @name profile.gcTMB
#' @aliases profile.gcTMB confint.profile.gcTMB plot.profile.gcTMB
profile.gcTMB <- function(fitted, parm=NULL, level=0.99, trace=FALSE, progress=interactive(), ...){
  est <- coef(fitted)
  pn <- names(est)
  nb <- sum(!fitted$fix)
  if (is.null(parm)) parm <- seq_len(nb)
  if (is.character(parm)){
    id <- match(parm, pn)
    if (any(is.na(id))) stop("unknown parameters in parm: ", paste(parm[is.na(id)], collapse=", "))
    parm <- id
  }
  if (any(parm < 1 | parm > length(est))) stop("parm must select estimated parameters of coef(fitted).")

  # range of the profiled parameter: the bounds of gcTMB (infinite for the variance parameters)
  lower <- if (is.null(fitted$lower)) rep(-Inf, length(est)) else fitted$lower
  upper <- if (is.null(fitted$upper)) rep(Inf, length(est)) else fitted$upper

  nllmin <- fitted$opt$objective
  ytol <- profcutoff(level, fitted)/2 + 0.1
  pb <- progressbar(length(parm), progress)
  prof <- lapply(parm, function(i){
    p <- TMB::tmbprofile(fitted$obj, name=i, ytol=ytol, parm.range=c(lower[i], upper[i]), trace=trace, ...)
    d <- data.frame(par=p[[1]], deviance=2*(p$value - nllmin))
    pb$tick()
    d[order(d$par), ]
  })
  pb$close()
  names(prof) <- pn[parm]
  attr(prof, "estimates") <- est[parm]
  attr(prof, "corr2") <- corr2names(fitted)
  attr(prof, "cutoff") <- cutoffattr(fitted, parm <= nb)
  class(prof) <- "profile.gcTMB"
  prof
}

#' @rdname profile.gcTMB
#' @export
confint.profile.gcTMB <- function(object, parm=NULL, level=0.95, ...){
  if (is.null(parm)) parm <- names(object)
  if (is.numeric(parm)) parm <- names(object)[parm]
  est <- attr(object, "estimates")
  # limits: the parameter values at which the deviance reaches the cutoff, interpolated on each
  # side of the estimate; NA if the profile does not reach the cutoff
  logscale <- attr(object, "logscale")
  if (is.null(logscale)) logscale <- setNames(rep(FALSE, length(object)), names(object))
  limit <- function(d, e, side, nm){
    cutoff <- profcutoffs(object, level)[[nm]]
    ds <- d[is.finite(d$deviance) & (if (side == "lower") d$par <= e else d$par >= e), ]
    ds <- ds[order(abs(ds$par - e)), ]
    k <- which(ds$deviance >= cutoff)[1]
    if (is.na(k) || k == 1) return(NA)
    # the square root of the deviance is close to linear in the parameter, on the scale the
    # profile was computed on: spline through the points up to the cutoff (and one beyond), or
    # linear interpolation where the profile is not monotone
    tr <- if (logscale[[nm]]) log else identity
    itr <- if (logscale[[nm]]) exp else identity
    z <- sqrt(pmax(ds$deviance, 0))
    x <- tr(ds$par)
    use <- seq_len(min(k + 1, nrow(ds)))
    if (all(diff(z[use]) > 0)){
      itr(splinefun(z[use], x[use], method="monoH.FC")(sqrt(cutoff)))
    } else {
      itr(approx(z[c(k - 1, k)], x[c(k - 1, k)], xout=sqrt(cutoff))$y)
    }
  }
  ci <- t(sapply(parm, function(nm) c(limit(object[[nm]], est[[nm]], "lower", nm), limit(object[[nm]], est[[nm]], "upper", nm))))
  rn <- parm

  # a flat profile towards small log standard deviations: lower limit 0 for the standard deviation
  logsd <- grepl("^log_sd:", rn)
  ci[logsd & is.na(ci[, 1]), 1] <- -Inf
  if (any(is.na(ci))){
    warning("the profiles of ", paste(rn[apply(is.na(ci), 1, any)], collapse=", "),
            " do not reach the confidence level; compute them with a higher level in profile.", call.=FALSE)
  }

  natscale(ci, rn, attr(object, "corr2"), level)
}

# confidence intervals on the natural scale: standard deviations and precision from their logs,
# correlations of terms with two random effects from their parameter; rows named accordingly
natscale <- function(ci, rn, corr2, level){
  ci <- matrix(ci, ncol=2)
  tolog <- grepl("^log_sd:", rn) | rn %in% c("log_sigma", "log_phi", "log_theta")
  ci[tolog, ] <- exp(ci[tolog, ])
  rn[tolog] <- sub("^log_", "", rn[tolog])
  two <- rn %in% corr2
  ci[two, ] <- ci[two, ] / sqrt(1 + ci[two, ]^2)
  rn[two] <- sub("^corr_par:", "corr:", rn[two])
  dimnames(ci) <- list(rn, paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%"))
  ci
}

# correlation parameters of terms with two random effects, which determine the correlation
corr2names <- function(fit){
  unlist(lapply(fit$reterms, function(tm){
    if (length(tm$pars) == 2) paste("corr_par:", relabels(tm)[2], ",", relabels(tm)[1], "|", tm$group, sep="")
  }))
}

#' @rdname profile.gcTMB
#' @export
plot.profile.gcTMB <- function(x, level=0.95, ...){
  nm <- names(x)
  d <- do.call(rbind, lapply(nm, function(n) data.frame(parameter=n, par=x[[n]]$par, deviance=x[[n]]$deviance)))
  d$parameter <- factor(d$parameter, levels=nm)
  est <- data.frame(parameter=factor(nm, levels=nm), par=unname(attr(x, "estimates")))
  ggplot2::ggplot(d, ggplot2::aes(x=par, y=deviance)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(data=data.frame(parameter=factor(nm, levels=nm), cutoff=profcutoffs(x, level)[nm]),
                        ggplot2::aes(yintercept=cutoff), linetype="dashed") +
    ggplot2::geom_vline(data=est, ggplot2::aes(xintercept=par), linetype="dotted") +
    ggplot2::facet_wrap(~ parameter, scales="free_x") +
    ggplot2::labs(x=NULL, y="Deviance")
}

# columns of the profile data used in the aesthetics of plot.profile.gcTMB
utils::globalVariables(c("par", "deviance", "cutoff"))

# cutoff of the profile deviance for confidence level level: the chi-square quantile with one degree
# of freedom, or with small-sample inference (rdf of the fit) the F-based cutoff, equivalent to a t
# interval: the profile deviance of a gaussian model is sum(w) log(RSS/RSSmin), so
# (RSS - RSSmin)/(RSSmin/rdf) <= F(1, rdf) gives sum(w) log(1 + F(1, rdf)/rdf)
profcutoff <- function(level, fit=NULL){
  if (is.null(fit$rdf)) return(qchisq(level, 1))
  sum(fit$obj$env$data$w) * log(1 + qf(level, 1, fit$rdf)/fit$rdf)
}

# what a profile needs for its cutoffs: the fit's residual degrees of freedom and sum of weights,
# and which profiled parameters use them (coefficients and times, not variance parameters)
cutoffattr <- function(fit, small){
  if (is.null(fit$rdf)) return(NULL)
  list(rdf=fit$rdf, sumw=sum(fit$obj$env$data$w), small=small)
}

# cutoff for each profiled parameter
profcutoffs <- function(object, level){
  ca <- attr(object, "cutoff")
  cut <- setNames(rep(qchisq(level, 1), length(object)), names(object))
  if (!is.null(ca)) cut[ca$small] <- ca$sumw * log(1 + qf(level, 1, ca$rdf)/ca$rdf)
  cut
}

#' Likelihood profiles of times to reach a level of the asymptote
#'
#' Profiles the (marginal) log-likelihood for the times estimated by \code{\link{timeto}}. At each
#' value of a time, the model is reparameterised: the time b2 + b3 Finv(q) can be solved for b2, so one
#' coefficient of b2 is computed from the time and the other coefficients, which are re-optimised.
#'
#' The profiles are for the times of the population curves, without the random effects; the
#' variation between groups (\code{random=TRUE} in \code{\link{timeto}}) and the residual variation are
#' not included. Gaussian models without random effects use the F-based cutoff of
#' \code{\link{profile.gcTMB}} by default.
#'
#' @param fitted A timetoTMB object, the result of \code{\link{timeto}}.
#' @param parm Character or integer vector, or NULL. Times to profile, names or positions in \code{coef(fitted)}; default all
#' @param level Numeric between 0 and 1. The highest confidence level the profiles are computed for.
#' @param maxit Integer. Maximum number of steps in each direction
#' @param progress Logical. Show a progress bar over the profiled times
#' @param ... Not used.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25)
#' pr <- profile(timeto(fit, respLev=c(10, 90)), progress=FALSE)
#' confint(pr)
#' plot(pr)
#'
#' @return a profile.timetoTMB object, which has the confint and plot methods of
#'   \code{\link{profile.gcTMB}}
#' @export
profile.timetoTMB <- function(fitted, parm=NULL, level=0.99, maxit=100, progress=interactive(), ...){
  fit <- fitted$fit
  if (is.null(fit)) stop("compute the times again with timeto() to profile them.")
  est <- coef(fitted)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown times in parm.")
  if (any(timetoflag(fitted, "conditional")[parm])) stop("profiles are not available for the times of fitted groups.")
  if (anyNA(est[parm])) stop("the times ", paste(names(est)[parm][is.na(est[parm])], collapse=", "), " are NA.")
  absolute <- identical(fitted$type, "absolute")

  obj <- fit$obj
  # keep the state of the TMB object, as tmbprofile does
  keep <- c("last.par.best", "random.start", "value.best", "last.par", "inner.control", "tracemgc")
  old <- mget(keep, envir=obj$env, ifnotfound=list(NULL))
  on.exit(for (v in keep) assign(v, old[[v]], envir=obj$env))
  obj$env$inner.control$trace <- FALSE
  obj$env$tracemgc <- FALSE

  mod <- modelindex(fit$model)
  lnk <- linkindex(fit$link)
  free <- which(!fit$fix)
  nf <- length(free)
  par0 <- fit$opt$par
  nllmin <- fit$opt$objective
  cutoff <- profcutoff(level, fit)
  Hfull <- tryCatch(solve(fit$sdrl$cov.fixed), error=function(e) NULL)
  # bounds of the fit
  lo <- if (is.null(fit$lower)) rep(-Inf, length(par0)) else fit$lower
  up <- if (is.null(fit$upper)) rep(Inf, length(par0)) else fit$upper
  # the parameter computed from the time: b2
  ps <- 2
  oth <- setdiff(1:5, ps)

  pb <- progressbar(length(parm), progress)
  prof <- lapply(parm, function(i){
    L <- fitted$linfct[[i]]
    q <- fitted$q[[i]]
    # the coefficient of b2 computed from the time: free, with the largest weight
    w2 <- L[ps, free]
    if (all(w2 == 0)) stop("the time ", names(est)[i], " does not depend on an estimated coefficient of b2.")
    kp <- which.max(abs(w2))
    k <- free[kp]
    if (any(L[3:5, k] != 0)) stop("the coefficient ", fit$coefnames[k], " of b2 must not enter b3 in linfct.")
    if (absolute && L[1, k] != 0) stop("the coefficient ", fit$coefnames[k], " of b2 must not enter b1 in linfct.")

    # full parameter vector from the reduced one (without coefficient k)
    pvfull <- function(par){
      pv <- fit$start
      pv[free] <- par[seq_len(nf)]
      pv
    }
    # b2 on the scale of its coefficients, from b3 (and b1 for absolute times) in lpv, at time tt
    b2fun <- function(tt, lpv){
      b <- linkpars(lnk, lpv)
      qf <- timeq(fit, q, lpv, absolute)
      linkfun(lnk[ps], tt - b[3] * gcmodel(mod)$Finv(qf))
    }
    expand <- function(x, tt){
      par <- append(x, 0, after=kp - 1)
      pv <- pvfull(par)
      lpv <- as.vector(L %*% pv)
      b2 <- b2fun(tt, lpv)
      par[kp] <- (b2 - sum(L[ps, -k] * pv[-k])) / L[ps, k]
      list(par=par, lpv=lpv)
    }
    # derivatives of the full parameter vector w.r.t. the reduced one: coefficient k depends on the
    # other coefficients of b2 and, through the time, on b3, and for absolute times on b1
    jacobian <- function(e, tt){
      h <- 1e-6
      d <- sapply(oth, function(j){
        hj <- replace(numeric(5), j, h)
        (b2fun(tt, e$lpv + hj) - b2fun(tt, e$lpv - hj)) / (2*h)
      })
      J <- diag(length(par0))[, -kp, drop=FALSE]
      J[kp, seq_len(nf - 1)] <- ((colSums(d * L[oth, free, drop=FALSE]) - L[ps, free]) / L[ps, k])[-kp]
      J
    }
    fn <- function(x, tt){
      e <- expand(x, tt)
      # the computed coefficient within its bounds, and finite: absolute times need a response
      # between zero and the asymptote
      if (!isTRUE(e$par[kp] >= lo[kp] && e$par[kp] <= up[kp])) return(Inf)
      obj$fn(e$par)
    }
    gr <- function(x, tt){
      e <- expand(x, tt)
      as.vector(obj$gr(e$par) %*% jacobian(e, tt))
    }
    # the Hessian at the estimates, transformed to the reduced parameters, as a fixed approximation:
    # far fewer iterations than quasi-Newton steps from scratch
    he <- if (is.null(Hfull)) NULL else function(x, tt){
      J <- jacobian(expand(x, tt), tt)
      t(J) %*% Hfull %*% J
    }

    # steps of a fifth of the standard error
    t0 <- est[[i]]
    step <- sqrt(fitted$variance[i, i]) / 5
    side <- function(dir){
      x <- par0[-kp]
      out <- NULL
      for (it in seq_len(maxit)){
        tt <- t0 + dir*it*step
        # the profile ends where the parameterisation fails, e.g. an absolute time whose response
        # approaches the asymptote
        o <- tryCatch(suppressWarnings(nlminb(x, fn, gr, he, tt=tt, lower=lo[-kp], upper=up[-kp])), error=function(e) NULL)
        if (is.null(o) || !is.finite(o$objective)) break
        x <- o$par
        out <- rbind(out, data.frame(par=tt, deviance=2*(o$objective - nllmin)))
        if (out$deviance[nrow(out)] > cutoff) break
      }
      out
    }
    d <- rbind(side(-1), data.frame(par=t0, deviance=0), side(1))
    pb$tick()
    d[order(d$par), ]
  })
  pb$close()
  names(prof) <- names(est)[parm]
  attr(prof, "estimates") <- est[parm]
  attr(prof, "cutoff") <- cutoffattr(fit, rep(TRUE, length(parm)))
  class(prof) <- c("profile.timetoTMB", "profile.gcTMB")
  prof
}

# time of entry i of the timetoTMB object for the curve parameters lpv (b1-b5 on the scale of the linear
# predictors); NA if the curve does not reach the level
timevalue <- function(object, i, lpv){
  fit <- object$fit
  b <- linkpars(linkindex(fit$link), lpv)
  q <- timeq(fit, object$q[[i]], lpv, identical(object$type, "absolute"))
  if (!is.finite(q) || q <= 0 || q >= 1) return(NA)
  gctimeR(modelindex(fit$model), q, b[2], b[3])
}

#' Confidence intervals for gcTMB models and times to reach a level of the asymptote
#'
#' Profile likelihood intervals (default), Wald intervals, or parametric bootstrap intervals for the
#' coefficients and variance parameters of a gcTMB model, or for the times estimated by
#' \code{\link{timeto}}.
#'
#' Standard deviations of the random effects, the residual standard deviation and the precision are
#' reported on their natural scale; Wald intervals for them are computed on the log scale. The
#' correlation of a term with two random effects is reported as a correlation; the parameters of
#' larger correlation matrices on the unconstrained scale of \code{coef}.
#'
#' Profile and bootstrap intervals for the times are for the population curves, and bootstrap
#' intervals also for the curves of fitted groups (\code{newdata} in \code{\link{timeto}}), with the
#' random effects of the groups drawn from their conditional distribution at each bootstrap fit; for
#' these, profile intervals are not available and Wald intervals are returned with a message. Times
#' computed with \code{random=TRUE} for random effect terms without a fitted group, or with
#' \code{residual=TRUE}, include further variation in their standard errors; for these, Wald intervals
#' are returned with a message.
#' Bootstrap intervals need many fits of the model; with random effects this can take several minutes.
#' Computing profiles can take a while for models with random effects; see \code{\link{profile.gcTMB}}
#' and \code{\link{profile.timetoTMB}}, whose results can be reused with \code{confint}.
#'
#' @param object A gcTMB object (\code{\link{gcTMB}}) or a timetoTMB object (\code{\link{timeto}}).
#' @param parm Character or integer vector, or NULL. Parameters or times, names or positions in \code{coef(object)}; default all
#' @param level Numeric between 0 and 1. The confidence level
#' @param method Character string. "profile" for profile likelihood intervals, "wald" for Wald intervals, or "boot" for
#'   percentile intervals of a parametric bootstrap: responses are simulated from the fitted model
#'   (with new random effects, see \code{\link{simulate.gcTMB}}) and the model is fitted again to each,
#'   starting at the estimates. Fits that do not converge are left out, with a message.
#' @param nsim Integer. Number of bootstrap simulations
#' @param seed Integer or NULL. An optional seed for the random number generator of the bootstrap
#' @param boot A bootstrap.gcTMB object or NULL. A bootstrap sample of the model from \code{\link{bootstrap}}, used instead of new
#'   simulations; implies \code{method="boot"}
#' @param ... Further arguments passed to \code{\link{profile.gcTMB}}, \code{\link{profile.timetoTMB}} or
#'   \code{\link{bootstrap}}, e.g. \code{progress=FALSE}
#'
#' @return a matrix with the lower and upper confidence limits
#' @export
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25)
#' confint(fit)
#' confint(fit, method="wald")
#' confint(timeto(fit, respLev=c(10, 50, 90)))
#' \donttest{
#' # parametric bootstrap, reused for the times
#' bt <- bootstrap(fit, nsim=200, seed=1, progress=FALSE)
#' confint(fit, boot=bt)
#' confint(timeto(fit, respLev=c(10, 50, 90)), boot=bt)
#' }
#' @name confint.gcTMB
#' @aliases confint.gcTMB confint.timetoTMB
confint.gcTMB <- function(object, parm=NULL, level=0.95, method=c("profile", "wald", "boot"), nsim=500, seed=NULL, boot=NULL, ...){
  method <- if (!is.null(boot)) "boot" else match.arg(method)
  est <- coef(object)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown parameters in parm.")
  if (method == "profile"){
    message("Computing profile confidence intervals ...")
    return(confint(profile(object, parm=parm, level=level, ...), level=level))
  }
  if (method == "boot") return(confint(getboot(object, boot, nsim, seed, ...), parm=parm, level=level))
  se <- sqrt(diag(vcov(object)))[parm]
  # t quantiles for the coefficients b1-b5 with small-sample inference
  z <- ifelse(parm <= sum(!object$fix), zquant(object, level), qnorm((1 + level)/2))
  natscale(cbind(est[parm] - z*se, est[parm] + z*se), names(est)[parm], corr2names(object), level)
}

#' @rdname confint.gcTMB
#' @export
confint.timetoTMB <- function(object, parm=NULL, level=0.95, method=c("profile", "wald", "boot"), nsim=500, seed=NULL, boot=NULL, ...){
  method <- if (!is.null(boot)) "boot" else match.arg(method)
  est <- coef(object)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown times in parm.")
  if (method %in% c("profile", "boot")){
    fit <- object$fit
    if (is.null(fit)){
      message("Wald intervals: compute the times again with timeto() for ", method, " intervals.")
      method <- "wald"
    } else if (any(timetoflag(object, "newvar")[parm]) || isTRUE(object$residual)){
      message("Wald intervals: ", method, " intervals are only available for population times, ",
              "timeto(..., random=FALSE, residual=FALSE), or times of fitted groups.")
      method <- "wald"
    } else if (method == "profile" && any(timetoflag(object, "conditional")[parm])){
      message("Wald intervals: profile intervals are not available for the times of fitted groups; ",
              "use method=\"boot\" for bootstrap intervals.")
      method <- "wald"
    }
  }
  if (method == "profile"){
    message("Computing profile confidence intervals ...")
    return(confint(profile(object, parm=parm, level=level, ...), level=level))
  }
  if (method == "boot"){
    fit <- object$fit
    # parameters of the bootstrap fits, with random effects of the fitted groups drawn from their
    # conditional distribution at each fit for the times of fitted groups
    draws <- bootdraws(fit, getboot(fit, boot, nsim, seed, ...), any(timetoflag(object, "conditional")[parm]))
    r <- fit$obj$env$random
    nb <- sum(!fit$fix)
    # times of each bootstrap fit
    tms <- sapply(parm, function(i){
      apply(draws, 2, function(par){
        pv <- fit$start
        pv[!fit$fix] <- par[seq_len(nb)]
        lpv <- as.vector(object$linfct[[i]] %*% pv)
        if (isTRUE(object$conditional[[i]])) lpv <- lpv + as.vector(object$uweights[[i]] %*% par[r])
        timevalue(object, i, lpv)
      })
    })
    ci <- t(apply(matrix(tms, ncol=length(parm)), 2, quantile, probs=c((1 - level)/2, (1 + level)/2), na.rm=TRUE, names=FALSE))
    dimnames(ci) <- list(names(est)[parm], paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%"))
    return(ci)
  }
  se <- sqrt(diag(vcov(object)))[parm]
  z <- if (is.null(object$fit)) qnorm((1 + level)/2) else zquant(object$fit, level)
  ci <- cbind(est[parm] - z*se, est[parm] + z*se)
  dimnames(ci) <- list(names(est)[parm], paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%"))
  ci
}

# logical flags of the times of a timetoTMB object: conditional (on the random effects of fitted
# groups) or newvar (with the variance of the random effects of a new group); for objects without
# them, newvar for all times with random=TRUE and random effects
timetoflag <- function(object, flag){
  f <- object[[flag]]
  if (!is.null(f)) return(f)
  n <- length(object$estimates)
  if (flag == "newvar") rep(isTRUE(object$random) && length(object$fit$reterms) > 0, n) else rep(FALSE, n)
}
