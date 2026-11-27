#' gcTMB Methods
#'
#' @param x,object A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param newdata Data frame or NULL. The data for the predictions, with the time and the covariates of the model. If NULL, the original data is used.
#'   For simulate, the times, covariates and grouping variables of the observations to simulate, e.g. new times or new individuals; each level of a grouping variable is one curve, with random effects shared by its rows. With random=TRUE all levels get new random effects; with random=FALSE the levels in the fitted data keep their fitted random effects and the other levels get new ones. For binomial responses the number of trials is taken from the response in newdata (one for a binary response).
#' @param random Logical. If TRUE, the random effects are included in the predictions, fitted values and residuals (conditional residuals). TRUE corresponds to the innermost level of \code{predict} and \code{fitted} of nlme, and FALSE to \code{level=0} (\code{re.form=NA} of lme4). For simulate, if TRUE, new random effects are drawn for each simulation, otherwise the fitted random effects are kept.
#'   For prediction intervals, random=TRUE gives the interval of a new observation in the group of each row of newdata, and random=FALSE the interval of a new observation in a new group, including the variance of the random effects. Rows of newdata with a level that is not in the fitted data are treated as new groups for that random effect term.
#' @param se.fit Logical. If TRUE, predict returns the standard errors of the predicted means, by the delta method; with random=TRUE they include the uncertainty of the predicted random effects.
#' @param interval Character string. The type of interval of predict: "none", "confidence" for confidence intervals of the mean, computed on the log scale for poisson and lognormal and on the logit scale for binomial and beta responses, or "prediction" for prediction intervals of a new observation with weight 1. Gaussian and lognormal models without random effects use t quantiles with n - p degrees of freedom and the unbiased residual variance, see \code{\link{gcTMB}}. For the gaussian and lognormal (on the log scale) families without boot, the prediction intervals use the normal (or t) distribution, and the variance of the random effects of new groups is added by the delta method. Otherwise they are quantiles of new observations simulated at nsim draws from the normal approximation of the fixed and random effects, with the variance parameters kept at their estimates, or at the bootstrap fits (see boot); with new random effects for new groups; for a binomial response as proportions, with the number of trials from the response in newdata (one for a binary response).
#' @param level Numeric between 0 and 1. Confidence level of the intervals of predict
#' @param digits Integer. Number of significant digits in the printed output
#' @param n Integer or NULL. The maximum number of coefficients printed; NULL (default) for the option \code{gcTMB.print.max}, 100 unless set, e.g. \code{options(gcTMB.print.max = 200)}; \code{Inf} for all. When coefficients are left out, a note at the end says how many, and how to show or extract all of them
#' @param type Character string. For residuals, the type of residuals: "pearson" (default), the response residuals divided by the standard deviation of the fitted distribution and multiplied by the square root of the weights; "response", the observed minus the fitted values, as proportions for a binomial response; or "quantile", randomised quantile residuals (Dunn and Smyth, 1996), which are standard normal under the model. The types "response", "pearson" and "normalized" are those of \code{residuals} of nlme. "normalized" (gaussian and lognormal) are the standardised residuals decorrelated by the autocorrelation of the fit (see \code{\link{corStruct}}): each standardised residual minus rho^dt times the previous one of its group, divided by sqrt(1 - rho^(2 dt)), independent standard normal under the model, as in nlme; without autocorrelation they are the pearson residuals; for the other families they are the pearson residuals, with a warning. For the gaussian family, the quantile residuals are multiplied by the square root of the weights, as precision weights with standard deviation sigma/sqrt(w); for the other families they do not use the weights, as a weighted log-likelihood is not the density of a distribution of the observation. The quantile residuals are random for the count families; use set.seed for reproducible values. For the lognormal family, the pearson and quantile residuals are the standardised residuals of log(y), (log(y) - log(f))/sigma, multiplied by the square root of the weights.
#'   For predict, the type of prediction: "response" (default), the mean response at the times of newdata, or "parameters", the parameters b1-b3 of the curve of each row of newdata (with its covariates in fform and, with random=TRUE, the random effects of its groups); newdata then needs no time. The standard errors of the parameters are those of the estimates, by the delta method, and their confidence intervals are computed on the scale of the link and back-transformed; with random=TRUE, they include the uncertainty of the predicted random effects of the fitted groups. Prediction intervals of the parameters add the variance of the random effects of the terms without a fitted level, the parameters of a new group, and are not available with boot.
#' @param nsim Integer. Number of response vectors to simulate; for predict, the number of parameter draws for simulated prediction intervals
#' @param seed Integer or NULL. Seed for the random number generator; NULL (default) sets no seed: the draws continue from the current state of the random number generator
#' @param boot A bootstrap.gcTMB object or NULL. A bootstrap sample of the model from \code{\link{bootstrap}}. With NULL (default), the standard errors and confidence intervals use the delta method. With a bootstrap sample, predict gives the standard deviation of the bootstrap predictions as se.fit, percentile confidence intervals, and prediction intervals simulated at the bootstrap fits, including their variance parameters. The bootstrap responses have new random effects, so for random=TRUE the random effects of the fitted groups are drawn, for each bootstrap fit, from the Laplace approximation of their conditional distribution given the original data at the parameters of that fit.
#' @param scale Character string. For type="parameters", the scale of the parameters: "original" (default) or "link", the scale of the coefficients, see the argument link of \code{\link{gcTMB}}.
#' @param ... Not used.
#'
#' @return Estimated model coefficients, variance-covariance matrix, or predictions.
#' predict returns a vector of predicted means, or, with se.fit=TRUE or an interval, a data frame with
#' the predicted means (fit), their standard errors (se.fit), and the limits of the intervals (lwr, upr).
#' With type="parameters", it returns a matrix of the parameters used by the model, a row for each row of
#' newdata, or, with se.fit=TRUE or an interval, a data frame with the row of newdata (row), the
#' parameter, and fit, se.fit, lwr and upr.
#' coef and vcov include the variance parameters (unlike \code{coef} of nlme and lme4, coef does not return the coefficients of each group; see \code{predict(type="parameters")}): for each random effect term the log standard deviations
#' (log_sd:b1|group, or log_sd:b1:x|group for a random slope of b1 in x) and, for correlated random effects, the parameters of the correlation matrix
#' (corr_par:b2,b1|group, unconstrained elements of the lower triangle of its Cholesky factor, by row);
#' and the log residual standard deviation (log_sigma, gaussian; of log(y) for lognormal), log precision (log_phi, beta), log size
#' (log_theta, nbinom2), log dispersion (log_phi, nbinom1) or log precision (log_phi, betabinomial); and
#' for an estimated autocorrelation the parameter cor_par:AR1 or cor_par:CAR1, or cor_par:range and cor_par:nugget, on the unconstrained scale of
#' nlme (\code{coef} of an nlme corStruct), see \code{\link{corStruct}}.
#' summary returns a summary.gcTMB object with the coefficient table (Estimate, Std. Error, and the
#' Wald test of a coefficient of zero: the z value and its p-value, or the t value with small-sample
#' inference, see \code{\link{gcTMB}}), the variance components, log-likelihood, AIC and BIC, sample
#' sizes and convergence information. Whether the family suits the data, e.g. overdispersion of
#' counts, is assessed with the quantile residuals (\code{residuals(fit, type="quantile")} and
#' \code{\link{plot.gcTMB}}), which are standard normal under the model for every family. The p-values
#' of the coefficient table are not adjusted for the number of tests. With treatment contrasts in
#' fform, the coefficients of a parameter are the differences of the levels to the reference level;
#' the tests of the differences between curves adjusted for their number are done by
#' \code{\link{partest}} (by default each level against a reference, by Dunnett's method), and the
#' test that a parameter is the same for all curves by \code{\link{anova.gcTMB}}. logLik returns the log-likelihood at the
#' estimates, with the number of free coefficients and variance parameters as degrees of freedom; it is the
#' maximum likelihood, comparable with nlme fits with \code{method="ML"}, not with the REML fits of \code{nlme::lme} and \code{nlme::gls}.
#' df.residual returns the residual degrees of freedom n - p of gaussian and lognormal models without random
#' effects (small-sample inference, see \code{\link{gcTMB}}), and NULL otherwise; as the \code{df}
#' argument of \code{multcomp::glht} for the times of \code{\link{timeto}}, it gives multivariate t instead of normal
#' quantiles.
#' simulate returns a data frame with nsim columns of simulated responses at the estimated
#' parameters (or those set by \code{\link{update.gcTMB}}), a row for each observation of the fit or row of newdata; for a binomial response given as cbind(successes, failures), the number of successes.
#' Its attribute \code{ranef} has the random effects of each simulation, a list named like the
#' columns with an element in the format of \code{\link{ranef}} (for random=FALSE, the fitted random
#' effects of the fitted levels), and its attribute \code{mean} the true means of the simulated
#' responses, a data frame like the simulations: the curve f(t) of each observation with its random
#' effects (the median for lognormal; for binomial and beta-binomial responses the proportion).
#'
#' @examples
#' fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
#' summary(fit)
#' coef(fit)
#' logLik(fit)
#' # the population curve with confidence and prediction intervals
#' nd <- data.frame(DAS=c(20, 40, 60))
#' predict(fit, newdata=nd, random=FALSE, interval="confidence")
#' predict(fit, newdata=nd, random=FALSE, interval="prediction")
#' # the curve of plant 6, and its parameters
#' predict(fit, newdata=data.frame(DAS=c(20, 40, 60), Plant_ID=6))
#' predict(fit, newdata=data.frame(Plant_ID=6), type="parameters")
#' head(residuals(fit, type="quantile"))
#' dim(simulate(fit, nsim=5, seed=1))
#' # three new plants at new times
#' nd <- expand.grid(DAS=seq(10, 80, 10), Plant_ID=c("new1", "new2", "new3"))
#' sims <- simulate(fit, nsim=2, seed=1, newdata=nd)
#' head(sims)
#' # the simulated random effects of the first simulation, and the true curves
#' attr(sims, "ranef")$sim_1
#' head(attr(sims, "mean"))
#'
#' @name gcTMB-methods
#' @aliases print.gcTMB summary.gcTMB print.summary.gcTMB logLik.gcTMB nobs.gcTMB df.residual.gcTMB fitted.gcTMB residuals.gcTMB predict.gcTMB coef.gcTMB vcov.gcTMB simulate.gcTMB
#' 
#' @rdname gcTMB-methods
#' @method print gcTMB
#' @export
print.gcTMB <- function(x, digits=max(3L, getOption("digits") - 3L), n=NULL, ...){
  if (!is.null(x$call)){
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("Model:  ", x$model, trunctext(x), "\n", sep="")
  cat("Family: ", x$family, ", ", linktext(x$link, x$model, fittrunc(x)), "\n", sep="")

  cat("\nCoefficients:\n")
  cf <- coef(x)[seq_len(sum(!x$fix))]
  k <- printrows(n, length(cf))
  print.default(format(cf[seq_len(k)], digits=digits), print.gap=2L, quote=FALSE)
  printomitted(k, length(cf), "coefficients", get="coef(x) returns all of them")
  printfixed(userfixed(x), digits)
  printre(recov(x), digits)
  printresid(x$family, residpar(x), digits)
  printvar(x, digits)
  printcor(x, digits)
  if (isTRUE(x$modified)) cat("\nParameters set by update(): not estimates, without standard errors.\n")
  if (!is.null(x$fiterror)) cat("\nThe fit failed: ", paste(x$fiterror, collapse="; "), "\n", sep="")
  invisible(x)
}

#' @rdname gcTMB-methods
#' @method summary gcTMB
#' @export
summary.gcTMB <- function(object, ...){
  dat <- object$obj$env$data
  m <- gcmodel(object$model)
  vcs <- recov(object)
  ll <- logLik(object)

  # coefficient table of the free coefficients
  ests <- object$ssdrl
  coefs <- ests[rownames(ests) %in% paste("b", 1:5, sep=""), , drop=FALSE]
  rownames(coefs) <- object$coefnames[!object$fix]
  # NaN for a Hessian that is not positive definite, with a warning in the printed summary
  coefs[, 2] <- suppressWarnings(sqrt(diag(fixedcov(object))))[seq_len(nrow(coefs))]
  # unadjusted Wald tests of a coefficient of zero: t tests with small-sample inference, z tests
  # otherwise; adjusted tests of the differences between curves are in partest
  stat <- coefs[, 1]/coefs[, 2]
  stat[!is.finite(coefs[, 2]) | coefs[, 2] <= 0] <- NA
  pval <- if (is.null(object$rdf)) 2*pnorm(-abs(stat)) else 2*pt(-abs(stat), object$rdf)
  coefs <- cbind(coefs[, 1:2, drop=FALSE], stat, pval)
  colnames(coefs) <- c("Estimate", "Std. Error", if (is.null(object$rdf)) c("z value", "Pr(>|z|)") else c("t value", "Pr(>|t|)"))

  # standard deviations at the boundary (singular fit): the likelihood is flat in the log standard
  # deviation, which gives a very large or no standard error, independent of the scale of the parameter
  ses <- ests[rownames(ests) == "theta", 2]
  toff <- 0
  singular <- unlist(lapply(vcs, function(vc){
    np <- length(vc$pars)
    se <- ses[toff + seq_len(np)]
    toff <<- toff + np + np*(np - 1)/2
    nm <- paste(vc$group, vc$labels)
    # parameters set by update(): no standard errors
    if (isTRUE(object$modified) || isFALSE(object$control$sdreport)) se <- rep(0, np)
    nm[!is.finite(se) | se > 100 | vc$sd < 1e-8]
  }))

  # number of levels of each grouping factor, once per factor
  ngroups <- sapply(object$reterms, function(tm) length(tm$levels))
  names(ngroups) <- sapply(object$reterms, function(tm) tm$group)
  ngroups <- ngroups[!duplicated(names(ngroups))]

  # gradient projected on the bounds: zero for coefficients at a bound that the likelihood pushes
  # against
  gr <- tryCatch(as.vector(object$obj$gr(object$opt$par)), error=function(e) rep(NaN, length(object$opt$par)))
  nb <- sum(!object$fix)
  if (!is.null(object$lower)){
    lo <- object$lower
    up <- object$upper
    gr[object$opt$par <= lo & gr > 0] <- 0
    gr[object$opt$par >= up & gr < 0] <- 0
  }
  res <- list(call=object$call, model=object$model, truncated=fittrunc(object), trunctext=trunctext(object), family=object$family, link=object$link, b3inverse=object$link[["b3"]] == "inverse",
              nobs=nobs(object), ntrials=if (object$family %in% c("binomial", "betabinomial")) sum(dat$bn), ngroups=ngroups,
              logLik=ll, AIC=AIC(object), BIC=BIC(object),
              coefficients=coefs, fixed=userfixed(object), varcor=vcs, residpar=residpar(object), rdf=object$rdf,
              singular=singular, atbound=atbound(object),
              convergence=object$opt$convergence, message=object$opt$message,
              maxgrad=max(abs(gr)), fiterror=object$fiterror, pdHess=object$sdrl$pdHess, modified=isTRUE(object$modified),
              nosdreport=isFALSE(object$control$sdreport), correlation=object$correlation, corgroup=object$corgroup,
              reterms=object$reterms, corrho=corpars(object), variance=object$variance, varpars=varpars(object))
  class(res) <- "summary.gcTMB"
  res
}

#' @rdname gcTMB-methods
#' @method print summary.gcTMB
#' @export
print.summary.gcTMB <- function(x, digits=max(3L, getOption("digits") - 3L), n=NULL, ...){
  if (!is.null(x$call)){
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("Model:  ", x$model, x$trunctext, "\n", sep="")
  cat("Family: ", x$family, ", ", linktext(x$link, x$model, fittrunc(x)), "\n", sep="")

  cat("\nNumber of observations: ", x$nobs, sep="")
  if (!is.null(x$ntrials)) cat(", number of trials: ", x$ntrials, sep="")
  if (length(x$ngroups) > 0) cat(", groups: ", paste(names(x$ngroups), x$ngroups, sep=" ", collapse="; "), sep="")
  cat("\n")
  cat("Log-likelihood: ", format(as.numeric(x$logLik), digits=digits), " (df = ", attr(x$logLik, "df"), ")",
      ", AIC: ", format(x$AIC, digits=digits), ", BIC: ", format(x$BIC, digits=digits), "\n", sep="")

  cat("\nCoefficients:\n")
  k <- printrows(n, nrow(x$coefficients))
  printCoefmat(x$coefficients[seq_len(k), , drop=FALSE], digits=digits)
  printomitted(k, nrow(x$coefficients), "coefficients", get="summary(fit)$coefficients returns the whole table")
  # covariates in fform: the adjusted tests of the differences between the curves are in partest
  if (!all(grepl(":\\(Intercept\\)$", rownames(x$coefficients))))
    cat("Unadjusted Wald tests; tests of the differences between curves adjusted for their number: partest().\n")
  if (length(x$atbound) > 0) cat("At a bound: ", paste(names(x$atbound), "=", format(x$atbound, digits=digits), collapse=", "),
                                  "; the standard errors and Wald intervals do not account for the bound.\n", sep="")
  if (isTRUE(x$b3inverse)) cat("The coefficients of b3 are on the inverse scale: the rate 1/b3.\n")
  for (lk in c("log", "logit")){
    pl <- names(x$link)[x$link == lk]
    if (length(pl) > 0) cat("The coefficients of ", paste(pl, collapse=", "), " are on the ", lk, " scale.\n", sep="")
  }
  for (p in names(x$link)[x$link == "logm1"]) cat("The coefficients of ", p, " are on the scale log(", p, " - 1).\n", sep="")
  printfixed(x$fixed, digits)
  printre(x$varcor, digits)
  if (length(x$singular) > 0){
    cat("Standard deviations close to zero (singular fit): ", paste(x$singular, collapse=", "), "\n", sep="")
  }
  printresid(x$family, x$residpar, digits, x$rdf)
  printvar(x, digits)
  printcor(x, digits)

  if (isTRUE(x$modified)){
    cat("\nParameters set by update(): not estimates, without standard errors.\n")
    return(invisible(x))
  }
  cat("\nOptimiser: ", x$message, ", maximum absolute gradient ", format(x$maxgrad, digits=2), "\n", sep="")
  # the warnings of the fit, repeated
  if (!is.null(x$fiterror)){
    cat("The fit failed: ", paste(x$fiterror, collapse="; "), "\n", sep="")
    return(invisible(x))
  }
  if (!is.finite(as.numeric(x$logLik))) cat("Warning: the negative log-likelihood is not finite at the estimates.\n")
  else if (!identical(as.integer(x$convergence), 0L)) cat("Warning: the optimiser did not converge: ", x$message, "\n", sep="")
  if (isTRUE(x$nosdreport)) cat("Standard errors not computed (control sdreport=FALSE).\n")
  else if (!isTRUE(x$pdHess)) cat("Warning: the Hessian is not positive definite; the standard errors are not reliable.\n")
  invisible(x)
}

#' @rdname gcTMB-methods
#' @export
logLik.gcTMB <- function(object, ...){
  # degrees of freedom: free coefficients and variance parameters
  structure(-object$opt$objective, df=length(object$opt$par), nobs=nobs(object), class="logLik")
}

#' @rdname gcTMB-methods
#' @export
nobs.gcTMB <- function(object, ...){
  length(object$obj$env$data$y)
}

#' @rdname gcTMB-methods
#' @export
df.residual.gcTMB <- function(object, ...){
  object$rdf
}

#' @rdname gcTMB-methods
#' @export
fitted.gcTMB <- function(object, random=TRUE, ...){
  predict(object, random=random)
}

#' @rdname gcTMB-methods
#' @export
residuals.gcTMB <- function(object, type=c("pearson", "response", "quantile", "normalized"), random=TRUE, ...){
  type <- match.arg(type)
  dat <- object$obj$env$data
  mu <- predict(object, random=random)
  cnt <- as.vector(dat$y)
  n <- as.vector(dat$bn)
  phi <- residpar(object)
  # binomial and beta-binomial: on the scale of proportions
  y <- if (object$family %in% c("binomial", "betabinomial")) cnt/n else cnt

  if (type == "response") return(y - mu)
  # square roots of the precision weights, divided by the standard deviation multiplier of the
  # variance function (gaussian and lognormal)
  sw <- sqrt(as.vector(dat$w))
  if (!is.null(object$variance)) sw <- sw/fitvarmult(object, random=TRUE)

  # normalized residuals: the standardised residuals of the gaussian and lognormal families, decorrelated
  # by the autocorrelation of the fit (Pinheiro and Bates, 2000); the pearson residuals otherwise
  if (type == "normalized"){
    if (!(object$family %in% c("gaussian", "lognormal"))){
      warning("normalized residuals are only available for the gaussian and lognormal families; ",
              "the pearson residuals of the ", object$family, " family are returned instead.", call.=FALSE)
      type <- "pearson"
    } else {
      r <- if (object$family == "lognormal") (log(y) - log(mu))/phi else (y - mu)/phi
      return(normalizeresid(object, r * sw))
    }
  }

  # lognormal: standardised residuals of log(y), for both pearson and quantile residuals
  if (object$family == "lognormal") return((log(y) - log(mu))/phi * sw)

  if (type == "pearson"){
    v <- switch(object$family,
                gaussian=phi^2,
                poisson=mu,
                binomial=mu*(1 - mu)/n,
                betabinomial=mu*(1 - mu)/n*(1 + (n - 1)/(phi + 1)),
                beta=mu*(1 - mu)/(1 + phi),
                nbinom2=mu*(1 + mu/phi),
                nbinom1=mu*(1 + phi))
    r <- (y - mu)/sqrt(v) * sw
    # a zero count at a mean of zero (a truncated curve before its threshold) has no variation
    r[mu == 0 & y == 0 & object$family != "gaussian"] <- 0
    return(r)
  }

  # randomised quantile residuals: normal quantiles of the fitted distribution function at the
  # observations, drawn uniformly between its left and right limits for the discrete families
  # gaussian: weights as precision weights, with standard deviation phi/sqrt(w), as for pearson
  if (object$family == "gaussian") return((y - mu)/phi * sw)
  if (object$family == "beta") return(qnorm(pbeta(y, mu*phi, (1 - mu)*phi)))
  if (object$family == "poisson"){
    lo <- ppois(cnt - 1, mu)
    hi <- ppois(cnt, mu)
  } else if (object$family %in% c("nbinom2", "nbinom1")){
    size <- nbsize(object$family, mu, phi)
    lo <- pnbinom(cnt - 1, size=size, mu=mu)
    hi <- pnbinom(cnt, size=size, mu=mu)
  } else if (object$family == "betabinomial"){
    lo <- pbetabinom(cnt - 1, n, mu, phi)
    hi <- pbetabinom(cnt, n, mu, phi)
  } else {
    lo <- pbinom(cnt - 1, n, mu)
    hi <- pbinom(cnt, n, mu)
  }
  qnorm(runif(length(y), lo, hi))
}

#' @rdname gcTMB-methods
#' @export
predict.gcTMB <- function(object, newdata=NULL, random=TRUE, se.fit=FALSE, interval=c("none", "confidence", "prediction"), level=0.95, boot=NULL, nsim=1000, seed=NULL,
                           type=c("response", "parameters"), scale=c("original", "link"), ...){
  interval <- match.arg(interval)
  type <- match.arg(type)
  scale <- match.arg(scale)
  if (is.null(newdata)) newdata <- object$data
  if (type == "parameters") return(predictpars(object, newdata, random, se.fit, interval, level, boot, nsim, seed, scale))
  pd <- preddata(object, newdata, random, covs=interval == "prediction")
  fit <- predmean(pd, pd$eta)
  if (!se.fit && interval == "none") return(fit)
  if (!is.null(seed)) set.seed(seed)

  # standard errors and confidence intervals: delta method, or the curves at the bootstrap draws
  probs <- c((1 - level)/2, (1 + level)/2)
  if (is.null(boot)){
    se <- predictse(object, pd)
  } else {
    draws <- bootdraws(object, getboot(object, boot, nsim, seed), random)
    bc <- drawcurves(object, pd, draws)
    se <- apply(bc, 1, sd, na.rm=TRUE)
  }
  res <- data.frame(fit=fit, se.fit=se)
  z <- zquant(object, level)
  if (interval == "confidence"){
    res[c("lwr", "upr")] <- if (is.null(boot)) predci(fit, se, object$family, z)
                            else t(apply(bc, 1, quantile, probs=probs, na.rm=TRUE, names=FALSE))
  }

  if (interval == "prediction"){
    if (object$family == "gaussian" && is.null(boot)){
      # new observation: uncertainty of the mean, residual variance, and the variance of the random
      # effects of the terms without a fitted level (new groups, or all terms for random=FALSE)
      sdp <- sqrt(se^2 + (residpar(object)*fitvarmult(object, newdata, random, pd=pd, mu=fit))^2 + newgroupvar(object, pd))
      res$lwr <- fit - z*sdp
      res$upr <- fit + z*sdp
    } else if (object$family == "lognormal" && is.null(boot)){
      # the same on the log scale, with the standard errors of the log median by the delta method
      sdp <- sqrt((se^2 + newgroupvar(object, pd))/fit^2 + (residpar(object)*fitvarmult(object, newdata, random, pd=pd, mu=fit))^2)
      res$lwr <- fit*exp(-z*sdp)
      res$upr <- fit*exp(z*sdp)
    } else {
      # simulated new observations at parameter draws: from the bootstrap fits, or from the normal
      # approximation of the fixed and random effects
      trials <- predtrials(object, newdata)
      if (is.null(boot)) draws <- jointdraws(object, nsim)
      res[c("lwr", "upr")] <- predsim(object, pd, draws, trials, probs, vardraws=!is.null(boot), newdata=newdata)
    }
  }
  res
}

# design of the predictions for newdata: times, design matrices, linear predictors of b1-b5 (with the
# random effects of the fitted levels if random=TRUE), the random effects design Zp for the TMB
# template, and for each term the rows without a fitted level (newgroup)
preddata <- function(object, newdata, random, covs=FALSE){
  form <- object$formula
  form[2] <- NULL
  x <- model.frame(form, data=newdata, na.action=na.pass)[, 1]
  n <- length(x)
  Xs <- lapply(object$fterms, function(ft) designmatrix(ft, newdata))
  npi <- rep(1:5, sapply(Xs, ncol))
  ests <- object$ssdrl
  rn <- rownames(ests)
  cest <- object$start
  cest[!object$fix] <- ests[rn %in% paste("b", 1:5, sep=""), 1]
  eta <- sapply(1:5, function(i) as.vector(Xs[[i]] %*% cest[npi == i]))
  eta <- matrix(eta, nrow=n)

  # random effects of the fitted levels; levels not in the fitted data get zero
  rl <- relevels(object, newdata, use=random, covs=covs)
  Zl <- rl$Zl
  newgroup <- rl$newgroup
  re <- ests[rn %in% "u", 1]
  uoff <- 0
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    np <- length(tm$pars)
    U <- matrix(re[uoff + seq_len(nl*np)], nrow=nl)
    uoff <- uoff + nl*np
    eta <- addre(eta, tm, Zl[[t]], U)
    # random slopes with missing covariates
    if (random) eta[rowSums(is.na(rl$X[[t]])) > 0, unique(tm$pars)] <- NA
  }
  Zp <- if (length(Zl) > 0) do.call(cbind, unlist(Zl, recursive=FALSE)) else Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(n, 0))
  list(x=x, Xs=Xs, eta=eta, Zp=Zp, Zl=Zl, X=rl$X, newgroup=newgroup,
       mod=modelindex(object$model), trunc=fittrunc(object), b1rate=fitrate(object), lnk=linkindex(object$link))
}

# design of the random effects of the fitted levels in newdata. For each term: Zl, a design block for
# warning for the rows of newdata that give the grouping variables of a random effect term, but with a
# missing level or a level not in the fitted data (rl from relevels with skip=TRUE): their curves have no
# random effects of the term; what follows for the result is in msg
newgroupwarning <- function(object, newdata, rl, msg){
  for (t in seq_along(object$reterms)){
    gv <- all.vars(object$reterms[[t]]$formula)
    unknown <- all(gv %in% names(newdata)) & rl$newgroup[, t]
    if (any(unknown)){
      warning("levels of ", paste(gv, collapse=":"), " missing or not in the fitted data in rows ", paste(which(unknown), collapse=", "),
              " of newdata: ", msg, call.=FALSE)
    }
  }
  invisible()
}

# each column of the term (intercept, slopes), with a row for each row of newdata and a column for
# each fitted level; X, the columns of the term in each row (1 for the intercept, the covariates of
# the slopes), for the random effects of a new group; and newgroup, the rows without a fitted level.
# Rows with a level not in the fit or a missing level, all rows for use=FALSE, and with skip=TRUE all
# rows of terms whose grouping variables are not in newdata have no fitted level. Rows with missing
# covariates of the slopes have missing X; without the covariates in newdata, an error for use=TRUE
# or covs=TRUE, and missing X otherwise.
relevels <- function(object, newdata, use=TRUE, skip=FALSE, covs=TRUE){
  n <- nrow(newdata)
  Zl <- list()
  Xl <- list()
  newgroup <- matrix(TRUE, n, length(object$reterms))
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    nc <- length(tm$cnms)
    Zl[[t]] <- rep(list(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(n, nl))), nc)
    Xl[[t]] <- matrix(NA_real_, n, nc)
    cmiss <- setdiff(all.vars(tm$xterms$terms), names(newdata))
    if (length(cmiss) > 0){
      if (use || covs) stop("newdata needs the covariates of the random slopes: ", paste(cmiss, collapse=", "), ".")
      next
    }
    # columns of the term with the factor levels and contrasts of the fit
    X <- designmatrix(tm$xterms, newdata)
    # level of each row from the intercept of the term; rows without a fitted level get a level
    # that is not in the fit
    bar <- str2lang(tm$name)
    bar[[2]] <- 1
    nd <- newdata
    for (v in all.vars(bar[[3]])){
      if (use && !skip && !(v %in% names(nd))) stop("newdata needs the grouping variable ", v, " of the random effects, or random=FALSE.")
      g <- if (use && v %in% names(nd)) as.character(nd[[v]]) else rep(NA_character_, n)
      g[is.na(g)] <- " new"
      nd[[v]] <- g
    }
    rt <- mkReTrms(list(bar), model.frame(subbars(as.formula(call("~", bar))), data=nd))
    lid <- match(as.character(rt$flist[[1]]), tm$levels)
    ok <- which(!is.na(lid))
    J <- Matrix::sparseMatrix(i=ok, j=lid[ok], x=1, dims=c(n, nl))
    # rows with missing covariates: missing X, and no random effects in the blocks
    Zl[[t]] <- lapply(seq_len(nc), function(c) Matrix::Diagonal(x=ifelse(is.na(X[, c]), 0, X[, c])) %*% J)
    Xl[[t]] <- X
    newgroup[, t] <- is.na(lid)
  }
  list(Zl=Zl, X=Xl, newgroup=newgroup)
}

# add the random effects U of term tm (a column for each random coefficient, a row for each level)
# to the linear predictors eta of b1-b5, with the design blocks Zb of the term
addre <- function(eta, tm, Zb, U){
  for (k in seq_along(tm$pars)) eta[, tm$pars[k]] <- eta[, tm$pars[k]] + as.vector(Zb[[tm$zcol[k]]] %*% U[, k])
  eta
}

# weights of the random coefficients of term tm on the linear predictors of their parameters for a
# new group, from the columns X of the term: a column for each coefficient
reloads <- function(tm, X){
  X[, tm$zcol, drop=FALSE]
}

# mean function at the linear predictors eta of b1-b5
predmean <- function(pd, eta){
  b <- linkpars(pd$lnk, eta)
  as.vector(gcmeanR(pd$mod, pd$x, b[, 1], b[, 2], b[, 3], b[, 4], b[, 5], trunc=isTRUE(pd$trunc), b1rate=isTRUE(pd$b1rate)))
}

# standard errors of the predicted means by the delta method of TMB: the mean at the new times is
# reported by the template, at the estimates; with random effects, sdreport includes the uncertainty
# of the predicted random effects
predictse <- function(object, pd){
  env <- object$obj$env
  dat <- env$data
  dat$xp <- pd$x
  for (i in 1:5) dat[[paste("Xp", i, sep="")]] <- pd$Xs[[i]]
  # TMB reads sparse matrices in triplet form
  dat$Zp <- as(as(pd$Zp, "TsparseMatrix"), "generalMatrix")
  rnd <- !is.null(env$random)
  obj <- MakeADFun(dat, env$parList(object$opt$par), map=env$map, random=if (rnd) "u",
                   DLL="gcTMB_TMBExports", silent=TRUE)
  # the likelihood is unchanged: the Hessian of the fit is reused
  H <- if (isTRUE(object$sdrl$pdHess)) solve(fixedcov(object))
  sdr <- sdreport(obj, par.fixed=object$opt$par, hessian.fixed=H, getReportCovariance=FALSE)
  s <- summary(sdr, "report")
  unname(s[rownames(s) == "mu", 2])
}

# parameters b1-b5 of the curves in the rows of newdata (predict with type="parameters"), on their
# original scale or on the scale of their links: a matrix with a column for each parameter used by the
# model, or with se.fit or an interval, a data frame with a row for each row of newdata and parameter.
# Standard errors on the link scale by the delta method, as the linear predictors are linear in the
# coefficients and random effects; intervals on the link scale, back-transformed
predictpars <- function(object, newdata, random, se.fit, interval, level, boot, nsim, seed, scale){
  # the parameters do not depend on the time
  for (v in setdiff(all.vars(object$formula[[3]]), names(newdata))) newdata[[v]] <- NA
  pd <- preddata(object, newdata, random, covs=interval == "prediction")
  n <- length(pd$x)
  used <- setdiff(1:5, unusedpars(pd$mod, isTRUE(pd$trunc)))
  # b4 and b5 only if they are not at their defaults
  for (k in 4:5) if (atdefault(object, k)) used <- setdiff(used, k)
  pnames <- paste("b", used, sep="")
  out <- function(eta) if (scale == "link") eta else linkpars(pd$lnk, eta)
  fit <- out(pd$eta)[, used, drop=FALSE]
  colnames(fit) <- pnames
  if (!se.fit && interval == "none"){
    # b1 as the maximum growth rate: the asymptote of the curves, on the original scale
    if (fitrate(object) && !islinear(pd$mod) && scale == "original"){
      b <- linkpars(pd$lnk, pd$eta)
      fit <- cbind(fit, asymptote=b[, 1] / rate1(pd$mod, b[, 2], b[, 3], b[, 4], b[, 5], isTRUE(pd$trunc)))
    }
    return(fit)
  }
  if (!is.null(seed)) set.seed(seed)

  z <- zquant(object, level)
  if (is.null(boot)){
    se <- predparse(object, pd)
    lo <- pd$eta - z*se
    hi <- pd$eta + z*se
    if (interval == "prediction"){
      # new groups: the variance of the random effects of the terms without a fitted level
      sdp <- sapply(1:5, function(k) sqrt(se[, k]^2 + newgroupvar(object, pd, outer(rep(1, n), 1:5 == k))))
      sdp <- matrix(sdp, nrow=n)
      lo <- pd$eta - z*sdp
      hi <- pd$eta + z*sdp
    }
    # standard errors on the original scale: times the derivative of the inverse link
    if (scale == "original") se <- se * sapply(1:5, function(k) dlinkinv(pd$lnk[k], pd$eta[, k]))
  } else {
    if (interval == "prediction") stop("prediction intervals of the parameters are not available with boot.")
    draws <- bootdraws(object, getboot(object, boot, nsim, seed), random)
    # linear predictors at the bootstrap fits, an array of rows, parameters and fits
    bl <- lapply(seq_len(ncol(draws)), function(j) draweta(object, pd, draws[, j]))
    be <- array(unlist(bl), c(n, 5, length(bl)))
    bo <- array(unlist(lapply(bl, out)), dim(be))
    se <- apply(bo, 1:2, sd, na.rm=TRUE)
    # percentile intervals: the quantiles of the linear predictors, back-transformed below
    probs <- c((1 - level)/2, (1 + level)/2)
    lo <- apply(be, 1:2, quantile, probs=probs[1], na.rm=TRUE, names=FALSE)
    hi <- apply(be, 1:2, quantile, probs=probs[2], na.rm=TRUE, names=FALSE)
  }
  se <- matrix(se, nrow=n)
  lo <- matrix(lo, nrow=n)
  hi <- matrix(hi, nrow=n)
  se[is.na(pd$eta)] <- NA
  res <- data.frame(row=rep(seq_len(n), length(used)), parameter=rep(pnames, each=n), fit=as.vector(fit),
                    se.fit=as.vector(se[, used]))
  if (interval != "none"){
    res$lwr <- as.vector(out(lo)[, used])
    res$upr <- as.vector(out(hi)[, used])
  }
  res
}

# standard errors of the linear predictors of b1-b5 in pd, a column for each parameter, by the delta
# method: they are linear in the free coefficients and, with random=TRUE, in the random effects of the
# fitted levels, with the joint covariance of both
predparse <- function(object, pd){
  env <- object$obj$env
  np <- length(env$par)
  nb <- sum(!object$fix)
  r <- env$random
  usere <- !is.null(r) && ncol(pd$Zp) > 0 && any(pd$Zp != 0)
  id <- c(seq_len(nb), if (usere) r)
  C <- if (usere) jointcov(object) else fixedcov(object)[seq_len(nb), seq_len(nb), drop=FALSE]
  # derivatives of the linear predictors w.r.t. each parameter in id; draweta includes the fixed coefficients
  e0 <- draweta(object, pd, numeric(np))
  J <- lapply(id, function(j) draweta(object, pd, replace(numeric(np), j, 1)) - e0)
  se <- sapply(1:5, function(k){
    Jk <- matrix(sapply(J, function(Jj) Jj[, k]), nrow=length(pd$x))
    sqrt(rowSums((Jk %*% C) * Jk))
  })
  matrix(se, nrow=length(pd$x))
}

# linear predictors of b1-b5 at a parameter vector par (as obj$env$par), with the random effects of
# the fitted levels in pd (none for random=FALSE)
draweta <- function(object, pd, par){
  pl <- object$obj$env$parList(par=par)
  eta <- matrix(sapply(1:5, function(i) as.vector(pd$Xs[[i]] %*% pl[[paste("b", i, sep="")]])), ncol=5)
  uoff <- 0
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    np <- length(tm$pars)
    U <- matrix(pl$u[uoff + seq_len(nl*np)], nrow=nl)
    uoff <- uoff + nl*np
    eta <- addre(eta, tm, pd$Zl[[t]], U)
  }
  eta
}

# predicted means at the parameter draws, one column per draw
drawcurves <- function(object, pd, draws){
  bc <- apply(draws, 2, function(par) predmean(pd, draweta(object, pd, par)))
  # apply drops to a vector for a single prediction
  matrix(bc, nrow=length(pd$x))
}

# parameter vectors (as obj$env$par, one per column) of the converged bootstrap fits. With random
# effects and random=TRUE, the random effects of the fitted levels are drawn from their conditional
# distribution given the original data at the parameters of each bootstrap fit: the bootstrap
# responses have new random effects, so the random effects of the bootstrap fits are not those of
# the fitted groups. Otherwise the random effects are not used.
bootdraws <- function(object, boot, random){
  est <- bootconverged(boot)
  env <- object$obj$env
  r <- env$random
  if (is.null(r)) return(matrix(t(est), ncol=nrow(est)))
  draws <- matrix(0, length(env$par), nrow(est))
  draws[-r, ] <- t(est)
  if (random){
    obj <- MakeADFun(env$data, env$parList(object$opt$par), map=env$map, random="u", DLL=env$DLL, silent=TRUE)
    for (k in seq_len(nrow(est))) draws[r, k] <- conddraw(obj, est[k, ])
    ok <- !is.na(draws[r[1], ])
    if (sum(ok) < nrow(est)) message(nrow(est) - sum(ok), " of ", nrow(est), " conditional distributions of the random effects at the bootstrap fits failed; these fits are left out.")
    if (sum(ok) < 2) stop("fewer than 2 bootstrap fits are left for bootstrap intervals.")
    draws <- draws[, ok, drop=FALSE]
  }
  draws
}

# one draw of the random effects of obj from the Laplace approximation of their conditional
# distribution given the data, at the fixed parameters par: normal with the conditional mode as the
# mean and the Hessian of the negative log-likelihood in the random effects as the precision; NA if
# the inner optimisation fails
conddraw <- function(obj, par){
  r <- obj$env$random
  u <- tryCatch({
    if (!is.finite(suppressWarnings(obj$fn(par)))) stop("not finite")
    lp <- obj$env$last.par
    R <- chol(as.matrix(obj$env$spHess(lp, random=TRUE)))
    lp[r] + backsolve(R, rnorm(length(r)))
  }, error=function(e) NULL)
  if (is.null(u) || any(!is.finite(u))) rep(NA, length(r)) else u
}

# draws from the normal approximation of all parameters: the covariance of the fixed parameters,
# or with random effects the joint precision of the fixed and random effects from sdreport
jointdraws <- function(object, nsim){
  env <- object$obj$env
  r <- env$random
  if (is.null(r)){
    m <- object$opt$par
    R <- tryCatch(chol(fixedcov(object)), error=function(e) NULL)
    if (is.null(R)) stop("the covariance matrix of the estimates is not positive definite.")
    return(m + t(R) %*% matrix(rnorm(length(m)*nsim), nrow=length(m)))
  }
  sdr <- sdreport(object$obj, par.fixed=object$opt$par, getJointPrecision=TRUE)
  m <- numeric(length(env$par))
  m[-r] <- sdr$par.fixed
  m[r] <- sdr$par.random
  R <- tryCatch(chol(as.matrix(sdr$jointPrecision)), error=function(e) NULL)
  if (is.null(R)) stop("the joint precision matrix of the parameters is not positive definite.")
  m + backsolve(R, matrix(rnorm(length(m)*nsim), nrow=length(m)))
}

# correlation matrix of a random effect term from its unconstrained parameters, as UNSTRUCTURED_CORR
# in TMB: lower triangle of L filled by row, correlation of L L'
corrmatrix <- function(cpar, np){
  L <- diag(np)
  L[upper.tri(L)] <- cpar
  L <- t(L)
  S <- L %*% t(L)
  S / sqrt(diag(S) %o% diag(S))
}

# number of trials of the new observations of a binomial or beta-binomial model: from the response in
# newdata, or one for a binary response
predtrials <- function(object, newdata){
  if (!(object$family %in% c("binomial", "betabinomial"))) return(NULL)
  y <- tryCatch(model.response(model.frame(object$formula, data=newdata, na.action=na.pass)), error=function(e) NULL)
  if (is.matrix(y) && ncol(y) == 2) return(rowSums(y))
  if (all(object$obj$env$data$bn == 1)) return(rep(1, nrow(newdata)))
  stop("prediction intervals of a binomial model need the number of trials: give the response, e.g. cbind(successes, failures), in newdata.")
}

# quantiles of new observations simulated at each parameter draw (columns of draws): conditional on
# the drawn random effects of the fitted levels for random=TRUE, with new random effects for the
# rows without a fitted level; binomial observations as proportions. The variance parameters are taken
# from the draws for vardraws=TRUE (bootstrap fits), and otherwise kept at their estimates: their
# normal approximation is poor, e.g. for standard deviations close to zero
predsim <- function(object, pd, draws, trials, probs, vardraws=FALSE, newdata=NULL){
  env <- object$obj$env
  n <- length(pd$x)
  pest <- env$parList(object$opt$par)
  sims <- apply(draws, 2, function(par){
    eta <- draweta(object, pd, par)
    if (vardraws) pest <- env$parList(par=par)
    toff <- 0
    for (t in seq_along(object$reterms)){
      tm <- object$reterms[[t]]
      np <- length(tm$pars)
      sdt <- exp(pest$theta[toff + seq_len(np)])
      Ct <- corrmatrix(pest$theta[toff + np + seq_len(np*(np - 1)/2)], np)
      toff <- toff + np + np*(np - 1)/2
      new <- which(pd$newgroup[, t])
      if (length(new) > 0){
        E <- matrix(rnorm(length(new)*np), ncol=np) %*% chol(Ct) %*% diag(sdt, nrow=np)
        E <- E * reloads(tm, pd$X[[t]][new, , drop=FALSE])
        for (k in seq_len(np)) eta[new, tm$pars[k]] <- eta[new, tm$pars[k]] + E[, k]
      }
    }
    mu <- predmean(pd, eta)
    phi <- exp(pest$log_sigma)
    # the standard deviation multiplier of the variance function (gaussian, lognormal)
    if (!is.null(object$variance) && object$family %in% c("gaussian", "lognormal")){
      phi <- phi*varmult(object$variance, pest$var_par, object$varlevels, newdata, mu, function(){
        h <- ratestep(pd$x)
        p1 <- p2 <- pd
        p1$x <- pd$x + h
        p2$x <- pd$x - h
        (predmean(p1, eta) - predmean(p2, eta))/(2*h)
      })
    }
    # draws outside the valid range of the parameters give NA
    suppressWarnings(switch(object$family,
                            gaussian=rnorm(n, mu, phi),
                            lognormal=rlnorm(n, log(mu), phi),
                            poisson=rpois(n, mu),
                            binomial=rbinom(n, trials, mu)/trials,
                            betabinomial=rbinom(n, trials, rbeta(n, mu*phi, (1 - mu)*phi))/trials,
                            beta=rbeta(n, mu*phi, (1 - mu)*phi),
                            rnbinom(n, size=nbsize(object$family, mu, phi), mu=mu)))
  })
  sims <- matrix(sims, nrow=n)
  bad <- mean(!is.finite(sims))
  if (bad > 0) warning(format(100*bad, digits=2), "% of the simulated observations are not finite and are left out; ",
                       "parameter draws may fall outside the valid parameter range.", call.=FALSE)
  # discrete responses: quantiles of the simulated values
  type <- if (object$family %in% c("poisson", "binomial", "betabinomial", "nbinom2", "nbinom1")) 1 else 7
  t(apply(sims, 1, quantile, probs=probs, na.rm=TRUE, names=FALSE, type=type))
}

# confidence intervals of the means, computed on the log scale for the count families (and the medians of
# lognormal responses) and the logit
# scale for binomial and beta, so that they stay inside the range of the mean
predci <- function(fit, se, family, z){
  tr <- switch(family,
               gaussian=list(f=identity, d=function(m) 1, i=identity),
               poisson=, nbinom2=, nbinom1=, lognormal=list(f=log, d=function(m) 1/m, i=exp),
               list(f=qlogis, d=function(m) 1/(m*(1 - m)), i=plogis))
  eta <- tr$f(fit)
  se.eta <- se*tr$d(fit)
  ci <- cbind(lwr=tr$i(eta - z*se.eta), upr=tr$i(eta + z*se.eta))
  # means at the boundary, e.g. an asymptote fixed at 0, without uncertainty
  pt <- se == 0 & !is.na(se)
  ci[pt, ] <- fit[pt]
  ci
}

# variance of the mean from the random effects of the terms without a fitted level in each row, by
# the delta method with numerical derivatives of the mean w.r.t. the linear predictors of b1-b5, or of a
# function with the derivatives g (a column for each of b1-b5)
newgroupvar <- function(object, pd, g=NULL){
  n <- length(pd$x)
  v <- numeric(n)
  if (length(object$reterms) == 0 || !any(pd$newgroup)) return(v)
  if (is.null(g)) g <- sapply(1:5, function(k){
    h <- 1e-5*pmax(1, abs(pd$eta[, k]))
    e1 <- e2 <- pd$eta
    e1[, k] <- e1[, k] + h
    e2[, k] <- e2[, k] - h
    (predmean(pd, e1) - predmean(pd, e2)) / (2*h)
  })
  g <- matrix(g, nrow=n)
  vcs <- recov(object)
  for (t in seq_along(vcs)){
    gt <- g[, vcs[[t]]$pars, drop=FALSE] * reloads(object$reterms[[t]], pd$X[[t]])
    vt <- rowSums((gt %*% vcs[[t]]$cov) * gt)
    v <- v + ifelse(pd$newgroup[, t], vt, 0)
  }
  v
}

#' @rdname gcTMB-methods
#' @export
coef.gcTMB <- function(object, ...){
    est <- object$estimates
    names(est) <- parnames(object)
    return(est)
}

#' @rdname gcTMB-methods
#' @export
vcov.gcTMB <- function(object, ...){
    vc <- fixedcov(object)
    rownames(vc) <- colnames(vc) <- names(coef(object))
    return(vc)
}

#' @rdname gcTMB-methods
#' @export
simulate.gcTMB <- function(object, nsim=1, seed=NULL, random=TRUE, newdata=NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(newdata)){
    return(simnew(object, newdata, nsim, random, seed))
  }
  obj <- object$obj
  par <- fitpar(object)
  if (!is.null(object$rform) && !random){
    # keep the fitted random effects: rebuild the objective with simre = 0
    dat <- obj$env$data
    dat$simre <- 0
    obj <- MakeADFun(dat, obj$env$parList(par=par), map=obj$env$map, random="u", DLL="gcTMB_TMBExports", silent=TRUE)
  }
  out <- lapply(seq_len(nsim), function(i) obj$simulate(par))
  sims <- lapply(out, function(o) o$y)
  nbad <- sum(sapply(sims, function(y) any(!is.finite(y))))
  if (nbad > 0) warning(nbad, " of ", nsim, " simulations contain non-finite values; simulated random effects may fall outside the valid parameter range.")
  simout(sims, lapply(out, function(o) o$fsim),
         lapply(out, function(o) ulist(object, o$usim, lapply(object$reterms, function(tm) tm$levels))), seed)
}

# the full parameter vector of the fit (as obj$env$par): the estimates, and the random effects at their
# conditional modes. Not obj$env$last.par.best, which TMB updates whenever the objective is evaluated at a
# better point (e.g. by a profile or another optimisation after a fit that stopped short), so that it
# need not be the fitted model.
fitpar <- function(object){
  env <- object$obj$env
  par <- env$par
  r <- env$random
  if (is.null(r)) return(object$opt$par)
  par[-r] <- object$opt$par
  par[r] <- object$sdrl$par.random
  par
}

# simulations as a data frame with nsim columns, with the attributes mean (the true means, a data
# frame like the simulations), ranef (for each simulation the simulated random effects, as ranef) and
# seed
simout <- function(sims, means, res, seed){
  nm <- paste("sim", seq_along(sims), sep="_")
  sims <- as.data.frame(sims, col.names=nm)
  attr(sims, "mean") <- as.data.frame(means, col.names=nm)
  attr(sims, "ranef") <- setNames(res, nm)
  attr(sims, "seed") <- seed
  sims
}

# random effects u (as in the template: by term, then by coefficient, then by level) as a list of data
# frames like ranef, with the levels lev of each term
ulist <- function(object, u, lev){
  uoff <- 0
  out <- lapply(seq_along(object$reterms), function(t){
    tm <- object$reterms[[t]]
    nl <- length(lev[[t]])
    np <- length(tm$pars)
    U <- matrix(u[uoff + seq_len(nl*np)], nrow=nl, dimnames=list(lev[[t]], relabels(tm)))
    uoff <<- uoff + nl*np
    as.data.frame(U, check.names=FALSE)
  })
  setNames(out, retermnames(object))
}


# responses simulated at the rows of newdata, nsim columns: the curves of the fixed effects plus random
# effects for each level of the grouping factors in newdata, one draw per level shared by its rows.
# random=TRUE draws new random effects for all levels; random=FALSE keeps the fitted random effects of
# the levels in the fit and draws new ones for the other levels. Observations with weight 1;
# binomial responses as numbers of successes, with the trials from the response in newdata.
simnew <- function(object, newdata, nsim, random, seed=NULL){
  newdata <- as.data.frame(newdata)
  pd <- preddata(object, newdata, random=FALSE)
  n <- length(pd$x)
  if (anyNA(pd$x)) stop("newdata needs the time of every row.")
  vcs <- recov(object)
  re <- if (!random) ranef.gcTMB(object)
  terms <- lapply(seq_along(object$reterms), function(t){
    tm <- object$reterms[[t]]
    g <- obsgroups(object, t, newdata)
    if (anyNA(g)) stop("newdata needs the grouping variables of the random effects, without missing values.")
    lv <- unique(g)
    load <- reloads(tm, designmatrix(tm$xterms, newdata))
    if (anyNA(load)) stop("newdata needs the covariates of the random slopes, without missing values.")
    # square root of the covariance matrix, also for zero variances
    e <- eigen(vcs[[t]]$cov, symmetric=TRUE)
    R <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow=length(e$values))
    fitted <- if (!random) match(lv, rownames(re[[t]]))
    list(tm=tm, id=match(g, lv), lv=lv, nl=length(lv), load=load, R=R, fitted=fitted, U=if (!random) as.matrix(re[[t]]))
  })
  trials <- predtrials(object, newdata)
  phi <- residpar(object)
  if (!is.null(object$rdf)) phi <- exp(object$sdrl$par.fixed[["log_sigma"]])
  out <- lapply(seq_len(nsim), function(i){
    eta <- pd$eta
    us <- numeric(0)
    for (tt in terms){
      np <- length(tt$tm$pars)
      U <- matrix(rnorm(tt$nl*np), ncol=np) %*% t(tt$R)
      if (!random){
        ok <- which(!is.na(tt$fitted))
        U[ok, ] <- tt$U[tt$fitted[ok], , drop=FALSE]
      }
      for (k in seq_len(np)) eta[, tt$tm$pars[k]] <- eta[, tt$tm$pars[k]] + tt$load[, k]*U[tt$id, k]
      us <- c(us, as.vector(U))
    }
    mu <- predmean(pd, eta)
    # residuals with the autocorrelation of the fit within groups (independent without one)
    if (object$family %in% c("gaussian", "lognormal")) z <- corsimresid(object, newdata) * fitvarmult(object, newdata, pd=pd, eta=eta, mu=mu)
    y <- suppressWarnings(switch(object$family,
                            gaussian=mu + phi*z,
                            lognormal=mu*exp(phi*z),
                            poisson=rpois(n, mu),
                            binomial=rbinom(n, trials, mu),
                            betabinomial=rbinom(n, trials, rbeta(n, mu*phi, (1 - mu)*phi)),
                            beta=rbeta(n, mu*phi, (1 - mu)*phi),
                            rnbinom(n, size=nbsize(object$family, mu, phi), mu=mu)))
    list(y=y, mu=mu, u=us)
  })
  sims <- lapply(out, function(o) o$y)
  nbad <- sum(sapply(sims, function(y) any(!is.finite(y))))
  if (nbad > 0) warning(nbad, " of ", nsim, " simulations contain non-finite values; simulated random effects may fall outside the valid parameter range.")
  lev <- lapply(terms, function(tt) tt$lv)
  simout(sims, lapply(out, function(o) o$mu), lapply(out, function(o) ulist(object, o$u, lev)), seed)
}

# standard deviations, correlation and covariance matrix of each random effect term
recov <- function(object){
  ests <- object$ssdrl
  rn <- rownames(ests)
  thetas <- ests[rn %in% "theta", 1]
  corrs <- ests[rn %in% "corr", 1]
  toff <- 0
  coff <- 0
  lapply(object$reterms, function(tm){
    np <- length(tm$pars)
    sdm <- exp(thetas[toff + seq_len(np)])
    cr <- matrix(corrs[coff + seq_len(np*np)], nrow=np)
    toff <<- toff + np + np*(np-1)/2
    coff <<- coff + np*np
    list(pars=tm$pars, labels=relabels(tm), group=tm$group, sd=sdm, corr=cr, cov=diag(sdm, nrow=np) %*% cr %*% diag(sdm, nrow=np))
  })
}


# free coefficients b1-b5 at their lower or upper bound, named
atbound <- function(object){
  nb <- sum(!object$fix)
  est <- object$opt$par[seq_len(nb)]
  lo <- if (is.null(object$lower)) rep(-Inf, nb) else object$lower[seq_len(nb)]
  up <- if (is.null(object$upper)) rep(Inf, nb) else object$upper[seq_len(nb)]
  # nlminb returns active bounds exactly
  at <- est <= lo | est >= up
  setNames(est[at], object$coefnames[!object$fix][at])
}

# links of a fit, as text: the identity link for all parameters, or the parameters with other links
# and the identity link for the others
linktext <- function(link, model, trunc=FALSE){
  used <- setdiff(1:5, unusedpars(modelindex(model), trunc))
  link <- link[used]
  other <- link[link != "identity"]
  if (length(other) == 0) return("identity link for all parameters")
  paste(paste(other, " link for ", names(other), sep="", collapse=", "), ", identity link for the other parameters", sep="")
}

# TRUE for a fit with the asymmetry b5 fixed at the symmetric curve, b5 = 1 for all curves (the default)
symmetric <- function(object) atdefault(object, 5)

# TRUE for a fit with all coefficients of parameter k (4 or 5) fixed at its default, b4 = 0 (the
# truncation at t = 0) or b5 = 1 (the symmetric curve), for all curves
atdefault <- function(object, k){
  npi <- rep(1:5, vapply(object$Xs, ncol, integer(1)))
  ik <- npi == k
  if (!all(object$fix[ik])) return(FALSE)
  b <- linkfun(linkindex(object$link[[k]]), as.vector(object$Xs[[k]] %*% object$start[ik]), inverse=TRUE)
  all(abs(b - c(0, 1)[k - 3]) < 1e-8)
}

# description of the truncation of a fit for print and summary, and of b1 as the maximum growth rate
trunctext <- function(object){
  tt <- if (!fittrunc(object) || isstart(modelindex(object$model))) NULL else if (atdefault(object, 4)) ", truncated at b4 = 0" else ", truncated at the threshold b4"
  if (fitrate(object) && !islinear(modelindex(object$model))) tt <- paste(tt, "; b1 is the maximum growth rate", sep="")
  tt
}

# coefficients fixed by the user, named; the parameters the model does not use are always fixed
userfixed <- function(object){
  npi <- as.integer(sub("^b([1-5]):.*", "\\1", object$coefnames))
  # b4 = 0 and b5 = 1, fixed by default, are not listed
  sym <- (npi == 4 & atdefault(object, 4)) | (npi == 5 & atdefault(object, 5))
  ufix <- object$fix & !(npi %in% unusedpars(modelindex(object$model), fittrunc(object))) & !sym
  setNames(object$start[ufix], object$coefnames[ufix])
}

printfixed <- function(fixed, digits){
  if (length(fixed) > 0){
    cat("Fixed: ", paste(names(fixed), "=", format(fixed, digits=digits), collapse=", "), "\n", sep="")
  }
}

# random effects: one row per parameter of each term, with standard deviation and correlations
# variance components (standard deviations of the random effects, parameters of the variance function and
# the correlation structure) for printing: two decimals each, and two significant digits
# for those that would round to 0.00 without being zero (e.g. on a log link)
dec2text <- function(s){
  vapply(s, function(v){
    if (!is.finite(v)) return(format(v))
    if (v != 0 && abs(v) < 0.005) formatC(v, digits=2, format="g") else formatC(v, digits=2, format="f")
  }, character(1))
}

printre <- function(vcs, digits){
  if (length(vcs) == 0) return(invisible())
  rows <- lapply(vcs, function(vc){
    np <- length(vc$pars)
    cr <- matrix("", np, max(np - 1, 0))
    for (i in seq_len(np)[-1]) cr[i, seq_len(i - 1)] <- formatC(vc$corr[i, seq_len(i - 1)], digits=2, format="f")
    cbind(Groups=c(vc$group, rep("", np - 1)), Name=vc$labels,
          Std.Dev.=dec2text(unname(vc$sd)), cr)
  })
  ncr <- max(sapply(rows, ncol)) - 3
  rows <- lapply(rows, function(r) cbind(r, matrix("", nrow(r), ncr - (ncol(r) - 3))))
  tab <- do.call(rbind, rows)
  colnames(tab) <- c("Groups", "Name", "Std.Dev.", if (ncr > 0) c("Corr", rep("", ncr - 1)))
  rownames(tab) <- rep("", nrow(tab))
  cat("\nRandom effects:\n")
  print(tab, quote=FALSE, right=FALSE)
}

# residual standard deviation (gaussian, of log(y) for lognormal), precision (beta, betabinomial), size theta (nbinom2) or dispersion phi
# (nbinom1); NULL for poisson and binomial; the unbiased estimate for small-sample inference
residpar <- function(object){
  if (object$family %in% c("poisson", "binomial")) return(NULL)
  unname(exp(object$sdrl$par.fixed["log_sigma"])) * sqrt(varscale(object))
}

# unbiased over maximum likelihood residual variance of gaussian (and lognormal) models without random effects: the
# ML estimate divides the weighted residual sum of squares by the sum of the weights, the unbiased
# one by the residual degrees of freedom rdf; 1 without small-sample inference
varscale <- function(object){
  if (is.null(object$rdf)) return(1)
  sum(object$obj$env$data$w) / object$rdf
}

# covariance of the estimated parameters: the inverse Hessian, with the block of the coefficients
# b1-b5 scaled by varscale for small-sample inference
fixedcov <- function(object){
  vc <- object$sdrl$cov.fixed
  s <- varscale(object)
  if (s != 1){
    d <- rep(1, nrow(vc))
    d[seq_len(sum(!object$fix))] <- sqrt(s)
    vc <- vc * outer(d, d)
  }
  vc
}

# quantile of two-sided intervals: t with the residual degrees of freedom for small-sample inference,
# otherwise normal
zquant <- function(object, level){
  if (is.null(object$rdf)) qnorm((1 + level)/2) else qt((1 + level)/2, object$rdf)
}

printresid <- function(family, residpar, digits, rdf=NULL){
  if (family == "gaussian") cat("\nResidual standard deviation: ", format(residpar, digits=digits),
                                if (!is.null(rdf)) paste(" on", rdf, "degrees of freedom"), "\n", sep="")
  if (family == "lognormal") cat("\nResidual standard deviation of log(y): ", format(residpar, digits=digits),
                                 if (!is.null(rdf)) paste(" on", rdf, "degrees of freedom"), "\n", sep="")
  if (family == "beta") cat("\nPrecision: ", format(residpar, digits=digits), "\n", sep="")
  if (family == "nbinom2") cat("\nNegative binomial size (theta): ", format(residpar, digits=digits), "\n", sep="")
  if (family == "betabinomial") cat("\nBeta-binomial precision (phi): ", format(residpar, digits=digits), "\n", sep="")
  if (family == "nbinom1") cat("\nNegative binomial dispersion (phi): ", format(residpar, digits=digits), "\n", sep="")
}

# names of the estimated parameters: the free coefficients; for each random effect term the log
# standard deviations (log_sd:b1|grp) and the parameters of the correlation matrix
# (corr_par:b2,b1|grp, for the pairs of the lower triangle by row, on the unconstrained scale of TMB);
# the log residual standard deviation (gaussian) or log precision (beta)
parnames <- function(object){
  renames <- unlist(lapply(object$reterms, function(tm){
    pn <- relabels(tm)
    np <- length(pn)
    pairs <- unlist(lapply(seq_len(np)[-1], function(i) paste(pn[i], pn[seq_len(i - 1)], sep=",")))
    c(paste("log_sd:", pn, "|", tm$group, sep=""),
      if (np > 1) paste("corr_par:", pairs, "|", tm$group, sep=""))
  }))
  # the autocorrelation, unless fixed: cor_par:AR1 or cor_par:CAR1, or cor_par:range (and cor_par:nugget)
  # of a spatial structure, on the unconstrained scale of nlme
  ct <- cortypeindex(object$correlation)
  cn <- if (ct > 0 && !isTRUE(attr(object$correlation, "fixed"))){
    if (ct >= 3) c("cor_par:range", if (corhasnugget(object$correlation)) "cor_par:nugget") else paste("cor_par:", c("AR1", "CAR1")[ct], sep="")
  }
  # the free parameters of the variance function, between log_sigma and cor_par as in the template
  vn <- if (!is.null(object$variance)){
    vsu <- varsetup(object$variance, object$data)
    vsu$names[!is.na(vsu$map)]
  }
  c(object$coefnames[!object$fix], renames,
    switch(object$family, gaussian=, lognormal="log_sigma", beta="log_phi", nbinom2="log_theta", nbinom1="log_phi", betabinomial="log_phi"),
    vn, cn)
}

# distribution function of the beta-binomial distribution with n trials, mean mu and precision phi
pbetabinom <- function(q, n, mu, phi){
  mapply(function(q, n, a, b){
    if (q < 0) return(0)
    k <- 0:min(q, n)
    sum(exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)))
  }, q, n, mu*phi, (1 - mu)*phi)
}

# size parameter of R's negative binomial functions for mean mu: theta for nbinom2, mu/phi for nbinom1
nbsize <- function(family, mu, phi){
  if (family == "nbinom2") rep(phi, length(mu)) else mu/phi
}
