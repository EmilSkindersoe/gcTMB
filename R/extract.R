#' Extract the components of a gcTMB fit
#'
#' gcTMB has its own generics fixef, ranef and VarCorr and does not need nlme or lme4. They have the
#' same names as the generics of nlme (re-exported by lme4), and work together with them: whichever
#' package is attached last, fixef, ranef and VarCorr work both for gcTMB fits and for the fits of
#' nlme and lme4 (the gcTMB generics pass other objects on to nlme, if it is installed). The results
#' for gcTMB fits have the format of lme4: \code{ranef} a list of data frames, one for each term, also
#' for a single term (nlme returns a data frame for a single level of grouping), and \code{VarCorr} a
#' list of covariance matrices (nlme returns a character table).
#'
#' Extractors for the estimated parts of a model: the fixed effects, the predicted random effects,
#' the covariance matrices of the random effects, the residual standard deviation, and the parameters
#' of the residual autocorrelation and of the variance function.
#'
#' @param object,x A gcTMB object, the result of \code{\link{gcTMB}}; for the generics
#'   \code{fixef}, \code{ranef} and \code{VarCorr}, a fitted model of any class.
#' @param sigma Not used; for compatibility with the generic of nlme.
#' @param ... Not used.
#'
#' @return
#' \code{fixef} returns the coefficients of the parameters b1-b5, named as in \code{coef}: the
#' estimated coefficients and the coefficients fixed by the user (e.g. an asymmetry b5 fixed at a
#' value other than 1), without the variance parameters and without the defaults b4 = 0 and b5 = 1.
#'
#' \code{ranef} returns a list with a data frame for each random effect term, named by its grouping
#' factor: a row for each level and a column for each random coefficient (e.g. b1, or b1:x for a
#' random slope of b1 in x). The random effects are the conditional modes, on the scale of the
#' linear predictors of the parameters (the scale of the link).
#'
#' \code{VarCorr} returns a list with the covariance matrix of the random effects of each term, named
#' as in \code{ranef}, with the standard deviations and the correlation matrix as attributes
#' \code{stddev} and \code{correlation}; the residual standard deviation is the attribute \code{sc}.
#'
#' \code{sigma} returns the residual standard deviation of a gaussian model (of log(y) for lognormal; the
#' maximum likelihood estimate, while \code{nlme::gnls} reports it with n - p degrees of freedom),
#' or the dispersion parameter of the other families as printed by \code{print.gcTMB} (precision,
#' size or dispersion), and NA for poisson and binomial models. The residual variance is
#' \code{sigma(object)^2}. For gaussian models without random effects and small-sample inference (see
#' \code{\link{gcTMB}}), it is the unbiased estimate.
#'
#' \code{corpars} returns the parameter of the autocorrelation of the residuals within groups (see
#' \code{\link{corStruct}}), rho, the correlation of residuals one unit of time apart, and an empty vector
#' for a model with independent residuals; it is \code{coef(fit$modelStruct$corStruct,
#' unconstrained=FALSE)} of nlme. \code{varpars} returns the parameters of the variance
#' function of the residuals, as a named vector; the current models have a constant residual variance,
#' so it returns an empty vector, as the extractor for a variance function once it is added.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' fixef(fit)
#' head(ranef(fit)$plant_ID)
#' VarCorr(fit)
#' sigma(fit)^2
#' corpars(fit)
#' varpars(fit)
#'
#' @name gcTMB-extract
#' @aliases fixef ranef VarCorr fixef.default ranef.default VarCorr.default fixef.gcTMB ranef.gcTMB VarCorr.gcTMB sigma.gcTMB
#' @importFrom stats ave var aggregate
#' @importFrom utils head tail
NULL

#' @rdname gcTMB-extract
#' @export
fixef <- function(object, ...) UseMethod("fixef")

#' @rdname gcTMB-extract
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' @rdname gcTMB-extract
#' @export
VarCorr <- function(x, sigma=1, ...) UseMethod("VarCorr")

#' @rdname gcTMB-extract
#' @export
fixef.default <- function(object, ...) othergeneric("fixef", object, ...)

#' @rdname gcTMB-extract
#' @export
ranef.default <- function(object, ...) othergeneric("ranef", object, ...)

#' @rdname gcTMB-extract
#' @export
VarCorr.default <- function(x, sigma=1, ...) othergeneric("VarCorr", x, sigma=sigma, ...)

# objects of other classes, e.g. fits of nlme or lme4 when gcTMB masks their generics: the generic
# of nlme, where nlme and lme4 register their methods
othergeneric <- function(name, object, ...){
  if (!requireNamespace("nlme", quietly=TRUE)) stop("no ", name, " method for an object of class ", class(object)[1], ".", call.=FALSE)
  getExportedValue("nlme", name)(object, ...)
}

#' @rdname gcTMB-extract
#' @export
#' @rawNamespace if (getRversion() >= "3.6.0") S3method(nlme::fixef, gcTMB)
fixef.gcTMB <- function(object, ...){
  b <- object$start
  b[!object$fix] <- object$estimates[seq_len(sum(!object$fix))]
  names(b) <- object$coefnames
  keep <- !object$fix | object$coefnames %in% names(userfixed(object))
  b[keep]
}

#' @rdname gcTMB-extract
#' @export
#' @rawNamespace if (getRversion() >= "3.6.0") S3method(nlme::ranef, gcTMB)
ranef.gcTMB <- function(object, ...){
  u <- object$sdrl$par.random
  uoff <- 0
  out <- lapply(object$reterms, function(tm){
    nl <- length(tm$levels)
    np <- length(tm$pars)
    U <- matrix(u[uoff + seq_len(nl*np)], nrow=nl, dimnames=list(tm$levels, relabels(tm)))
    uoff <<- uoff + nl*np
    as.data.frame(U, check.names=FALSE)
  })
  setNames(out, retermnames(object))
}

#' @rdname gcTMB-extract
#' @export
#' @rawNamespace if (getRversion() >= "3.6.0") S3method(nlme::VarCorr, gcTMB)
VarCorr.gcTMB <- function(x, sigma=1, ...){
  out <- lapply(recov(x), function(vc){
    nm <- vc$labels
    cv <- vc$cov
    cr <- vc$corr
    dimnames(cv) <- dimnames(cr) <- list(nm, nm)
    structure(cv, stddev=setNames(vc$sd, nm), correlation=cr)
  })
  structure(setNames(out, retermnames(x)), sc=sigma.gcTMB(x))
}

#' @rdname gcTMB-extract
#' @importFrom stats sigma
#' @export
sigma.gcTMB <- function(object, ...){
  s <- residpar(object)
  if (is.null(s)) NA_real_ else s
}

#' @rdname gcTMB-extract
#' @export
corpars <- function(object){
  ct <- cortypeindex(object$correlation)
  if (ct == 0) return(setNames(numeric(0), character(0)))
  v <- object$obj$env$parList(object$opt$par)$cor_par
  c(rho=corrho(ct, v))
}

#' @rdname gcTMB-extract
#' @export
varpars <- function(object){
  setNames(numeric(0), character(0))
}

# names of the random effect terms: the grouping factors, made unique for several terms of a factor
retermnames <- function(object){
  make.unique(vapply(object$reterms, function(tm) tm$group, character(1)))
}

# level of the grouping factor of random effect term t for each observation of the fit
obsgroups <- function(object, t, data=object$data){
  termgroups(object$reterms[[t]], data)
}

# level of the grouping factor of the random effect term tm for each row of data
termgroups <- function(tm, data){
  bar <- str2lang(tm$name)
  bar[[2]] <- 1
  rt <- mkReTrms(list(bar), model.frame(subbars(as.formula(call("~", bar))), data=data, na.action=na.pass))
  as.character(rt$flist[[1]])
}

#' Model selection criteria
#'
#' Information criteria and goodness-of-fit criteria of one or more gcTMB fits, for the choice
#' between growth models fitted to the same data.
#'
#' @param ... gcTMB objects, or a single list of gcTMB objects. The rows of the result are named by
#'   the names of the arguments or of the list, or numbered.
#' @param which Character vector or NULL. The criteria to return, in this order, e.g.
#'   \code{c("AIC", "R2", "BIC")}; NULL (default) returns all columns. A goodness-of-fit criterion
#'   without suffix, e.g. "R2", gives both its conditional and marginal version (R2.c and R2.m);
#'   "R2.m" gives only the marginal one.
#'
#' @details
#' With L the maximised log-likelihood, k the number of estimated parameters (free coefficients and
#' variance parameters, the df of \code{logLik}) and n the number of observations:
#' AIC = -2L + 2k; AICc = AIC + 2k(k + 1)/(n - k - 1); BIC = -2L + log(n)k;
#' BICc = BIC + k(k + 2)log(n)/(n - k - 2); evBIC = BIC + 2L/n;
#' HQIC = -2L + 2k log(log(n)); HQICc = -2L + 2k log(log(n)) n/(n - k - 2); and the bridge
#' criterion BC = -2L + n^(2/3) (1 + 1/2 + ... + 1/k). All are smaller for better models.
#'
#' The goodness-of-fit criteria come in a conditional version (suffix .c), with the curves of the
#' groups including their random effects, and a marginal version (suffix .m), with the population
#' curve of the fixed effects. The conditional versions use the residuals r = y - f(t) of the
#' observations; the marginal versions use the mean observation of each time and population curve
#' (the curves of the different values of the covariates in fform), r = mean(y) - f(t), for each
#' observation. RMSE is the root mean squared residual, ME the mean residual (the bias), and MAPE
#' 100 times the mean absolute residual relative to the (mean) observation. The root trapezoidal
#' squared error RTSE approximates the root mean integrated squared error between the observed and
#' fitted curves by the trapezoidal rule: the integral of r^2 over the times of each curve divided by
#' its time span, averaged over the curves, and the square root. The curves are the levels of the
#' random effect term with the most levels; the marginal RTSE uses the population curves, and the
#' conditional RTSE is NA without random effects.
#'
#' R2 is the coefficient of determination of Nakagawa and Schielzeth (2013): with the variance of the
#' marginal fitted values vf, the variance of the conditional minus the marginal fitted values vr
#' (zero without random effects) and the variance of the conditional residuals ve,
#' R2.m = vf/(vf + vr + ve) and R2.c = (vf + vr)/(vf + vr + ve). CCC is the concordance correlation
#' coefficient of Vonesh et al. (1996) between the observations and the conditional (CCC.c) or
#' marginal (CCC.m) fitted values. The adjusted versions R2adj and CCCadj are 1 - n/(n - p)(1 - R2)
#' with p the number of free coefficients. R2 and CCC are larger for better models.
#'
#' For binomial and beta-binomial responses the residuals are on the scale of proportions.
#'
#' @return A data frame with a row for each fit and the columns logLik, df, nobs, the information
#' criteria AIC, AICc, BIC, BICc, evBIC, HQIC, HQICc and BC, and the goodness-of-fit criteria
#' RMSE, ME, MAPE, RTSE, R2, R2adj, CCC and CCCadj, each with suffix .c (conditional) and .m
#' (marginal).
#'
#' @references
#' Nakagawa, S. and Schielzeth, H. (2013). A general and simple method for obtaining R2 from
#' generalized linear mixed-effects models. Methods in Ecology and Evolution 4, 133-142.
#'
#' Vonesh, E. F., Chinchilli, V. M. and Pu, K. (1996). Goodness-of-fit in generalized nonlinear
#' mixed-effects models. Biometrics 52, 572-587.
#'
#' @examples
#' fl <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' fg <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25, model="gompertz")
#' criteria(logistic=fl, gompertz=fg)
#' criteria(logistic=fl, gompertz=fg, which=c("AIC", "R2", "BIC"))
#'
#' @export
criteria <- function(..., which=NULL){
  fits <- list(...)
  if (length(fits) == 1 && is.list(fits[[1]]) && !inherits(fits[[1]], "gcTMB")) fits <- fits[[1]]
  if (!all(vapply(fits, inherits, logical(1), what="gcTMB"))) stop("the arguments must be gcTMB objects.")
  nm <- names(fits)
  if (is.null(nm)) nm <- rep("", length(fits))
  nm[nm == ""] <- as.character(seq_along(fits))[nm == ""]
  out <- do.call(rbind, lapply(fits, function(f) as.data.frame(as.list(criteria1(f)))))
  rownames(out) <- nm
  if (!is.null(which)){
    if (!is.character(which)) stop("which must be a character vector of criteria, e.g. c(\"AIC\", \"R2\").")
    # a name without suffix selects both the conditional and the marginal version
    cols <- lapply(which, function(w) if (w %in% names(out)) w else intersect(paste(w, c("c", "m"), sep="."), names(out)))
    bad <- which[lengths(cols) == 0]
    if (length(bad) > 0) stop("unknown criteria: ", paste(bad, collapse=", "), ". Available: ",
                              paste(unique(sub("\\.[cm]$", "", names(out))), collapse=", "), ".")
    out <- out[, unique(unlist(cols)), drop=FALSE]
  }
  out
}

# criteria of one fit, as a named vector
criteria1 <- function(object){
  ll <- logLik(object)
  L <- as.numeric(ll)
  k <- attr(ll, "df")
  n <- nobs(object)
  p <- sum(!object$fix)
  ic <- c(logLik=L, df=k, nobs=n,
          AIC=-2*L + 2*k,
          AICc=-2*L + 2*k + 2*k*(k + 1)/(n - k - 1),
          BIC=-2*L + log(n)*k,
          BICc=-2*L + log(n)*k + k*(k + 2)*log(n)/(n - k - 2),
          evBIC=-2*L + log(n)*k + 2*L/n,
          HQIC=-2*L + 2*k*log(log(n)),
          HQICc=-2*L + 2*k*log(log(n))*n/(n - k - 2),
          BC=-2*L + n^(2/3)*sum(1/seq_len(k)))

  t <- as.vector(object$obj$env$data$x)
  rnd <- length(object$reterms) > 0
  fc <- fitted(object, random=TRUE)
  fm <- fitted(object, random=FALSE)
  y <- fc + residuals(object, type="response", random=TRUE)

  # population curves: the different parameter values of the fixed effects
  pars <- predict(object, newdata=object$data, random=FALSE, type="parameters")
  pop <- apply(signif(pars, 12), 1, paste, collapse=" ")
  cell <- paste(pop, t)
  ybar <- ave(y, cell)

  rc <- y - fc
  rm <- ybar - fm
  # integral of r^2 over the times of each curve by the trapezoidal rule, over its time span
  rtse <- function(r, curve){
    ise <- vapply(split(seq_along(r), curve), function(i){
      d <- unique(data.frame(t=t[i], r2=r[i]^2))
      d <- aggregate(r2 ~ t, data=d, FUN=mean)
      if (nrow(d) < 2) return(NA_real_)
      sum(diff(d$t) * (head(d$r2, -1) + tail(d$r2, -1)) / 2) / diff(range(d$t))
    }, numeric(1))
    sqrt(mean(ise, na.rm=TRUE))
  }
  curve <- if (rnd){
    t0 <- which.max(vapply(object$reterms, function(tm) length(tm$levels), integer(1)))
    obsgroups(object, t0)
  }

  vf <- var(fm)
  vr <- if (rnd) var(fc - fm) else 0
  ve <- var(rc)
  r2m <- vf/(vf + vr + ve)
  r2c <- (vf + vr)/(vf + vr + ve)
  ccc <- function(y, f) 1 - sum((y - f)^2)/(sum((y - mean(y))^2) + sum((f - mean(f))^2) + length(y)*(mean(y) - mean(f))^2)
  adj <- function(r) 1 - n/(n - p)*(1 - r)
  cc <- ccc(y, fc)
  cm <- ccc(y, fm)

  c(ic,
    RMSE.c=sqrt(mean(rc^2)), RMSE.m=sqrt(mean(rm^2)),
    ME.c=mean(rc), ME.m=mean(rm),
    MAPE.c=100*mean(abs(rc/y)), MAPE.m=100*mean(abs(rm/ybar)),
    RTSE.c=if (rnd) rtse(rc, curve) else NA_real_, RTSE.m=rtse(rm, pop),
    R2.c=r2c, R2.m=r2m, R2adj.c=adj(r2c), R2adj.m=adj(r2m),
    CCC.c=cc, CCC.m=cm, CCCadj.c=adj(cc), CCCadj.m=adj(cm))
}
