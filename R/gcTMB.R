#' Fit a growth curve
#'
#' Fits a nonlinear growth curve over time by maximum likelihood with the Template Model Builder
#' (TMB), optionally with covariates on the parameters and random effects, which are integrated out by
#' the Laplace approximation.
#'
#' @param form Formula, with the names of the response and the time variable, e.g. \code{weight ~ day}
#' @param fform Formula, list of formulas, or NULL. Fixed effects for the model parameters: a formula of the form b1 + b2 ~ Treatment, or a list of such formulas, e.g. list(b1 ~ Treatment, b2 ~ Treatment + Site). Parameters not listed have a single intercept (~ 1). A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~ 1), is also accepted.
#' @param rform Formula, list of formulas, or NULL. Random effects: a formula of the form b1 + b2 ~ 1 | individual for random intercepts, or b1 ~ x | individual for random intercepts and slopes in a covariate x, or a list of such formulas, e.g. list(b1 + b2 ~ 1 | plant, b1 ~ 1 | plot). The random effects within a formula are correlated across the listed parameters and their intercepts and slopes; different formulas, the levels of nested factors (1 | plot/plant), and intercepts and slopes with || (b1 ~ x || individual) are independent terms. The random slopes of a parameter usually go with a fixed effect of the covariate in fform.
#' @param data Data frame. Observations with missing values in the variables of the model or in the weights are left out.
#' @param start Numeric vector or NULL. Starting values for each model coefficient. NA entries are replaced by the starting values of \code{\link{gcTMBstart}}. Alternatively, a vector named by parameter, e.g. c(b1=10), gives the same starting value to all coefficients of each named parameter.
#' @param fix Logical vector or NULL. One entry for each model coefficient. If TRUE, a model coefficient is fixed to its value in start. Alternatively, a vector named by parameter, e.g. c(b1=TRUE), fixes all coefficients of the named parameters; their values must then be given in start, e.g. start=c(b1=10). The threshold b4 and the asymmetry b5 are fixed by default, at b4 = 0 (the curve starts at zero at t = 0) and b5 = 1 (the symmetric curve) unless given in start; fix=c(b4=FALSE) or fix=c(b5=FALSE), or an unnamed fix vector, estimates them. Without truncation, b4 is not used and always fixed, at 0 unless given in start.
#' @param weights Numeric vector or NULL. Weights, one for each row of data.
#' @param lower Numeric, a single value or a vector. Lower bounds for the model coefficients b1-b5, a single value or a vector with the same length as start. Bounds of fixed coefficients are ignored; variance parameters are unbounded.
#' @param upper Numeric, a single value or a vector. Upper bounds for the model coefficients b1-b5, a single value or a vector with the same length as start.
#' @param family Character string. The distribution of the response, one of "gaussian", "lognormal", "binomial", "betabinomial", "beta", "poisson", "nbinom2", and "nbinom1". For the lognormal family, log(y) is normal with mean log f(t) and standard deviation sigma: the curve f is the median of a positive response, with a constant coefficient of variation; the mean is f(t) exp(sigma^2/2). The negative binomial families are for overdispersed counts, with variance mu (1 + mu/theta) for nbinom2 and mu (1 + phi) for nbinom1, and the beta-binomial family for overdispersed binomial counts, with variance n mu (1 - mu) (1 + (n - 1)/(phi + 1)) for n trials, as in glmmTMB.
#' @param model Character string. The growth model: "logistic" (default), "gompertz", "hill", "weibull", "normal" or "gamma", each with an asymmetry parameter b5, or "linear" for a straight line, see Details.
#' @param truncated Logical. If TRUE (default), the curve is truncated at a threshold time b4, by default fixed at 0: the response is zero up to b4 and then grows, see Details. FALSE for the curves without truncation.
#' @param link Named character vector or NULL. Links of the parameters b1-b5: NULL for the defaults (the identity link for all parameters), or links named by parameter, e.g. \code{link=c(b1="log", b3="inverse")}. The asymptote b1 has the link "identity" (default), "log" or "logit"; the times b2 and b4 (the threshold of a truncated curve) "identity" (default) or "log"; the asymmetry b5 "identity" (default) or "log", to keep it positive; and the scale b3 "identity" (default), "inverse", with which its coefficients are the growth rate 1/b3, or "log". The coefficients, and start, fix, lower, upper and curves of \code{\link{timeto}}, are on the scale of the links, see Details.
#' @param gsl Numeric vector or NULL. Lower bounds for a grid search, same length as start; NA for the default range. If gsl, gsu and gsn are all NULL, no grid search is performed. See \code{\link{gcTMBstart}}.
#' @param gsu Numeric vector or NULL. Upper bounds for a grid search, same length as start; NA for the default range.
#' @param gsn Integer, a single number or a vector, or NULL. Number of grid points for each parameter, a single number or a vector with the same length as start. \code{gsn=10} alone runs a grid search over the default ranges.
#' @param control List. Control parameters, the arguments of \code{\link{gcTMBcontrol}}, e.g. \code{list(iter.max=2000)}. See \code{\link{gcTMBcontrol}} for details.
#'
#' @details The growth models are of the form f(t) = b1 F(z(t)), with a distribution function F
#' increasing from 0 to 1 and a transformed time z(t):
#' \itemize{
#'   \item logistic: f(t) = b1 / (1 + exp(-(t - b2)/b3)), with the logistic F(z) = 1/(1 + exp(-z)) and
#'     z = (t - b2)/b3; symmetric around b2, where the curve reaches half of its asymptote;
#'   \item gompertz: f(t) = b1 exp(-exp(-(t - b2)/b3)), F(z) = exp(-exp(-z)) and z = (t - b2)/b3;
#'     it reaches exp(-1), about 37\%, of its asymptote at b2 and approaches the asymptote more slowly
#'     than it leaves zero;
#'   \item hill: f(t) = b1 / (1 + (b2/t)^b3), the logistic F of z = b3 (log t - log b2), a logistic
#'     curve on the log time axis (log-logistic); b2 is the time at which the curve reaches half of its
#'     asymptote, and b3 the slope on the log time axis;
#'   \item weibull: f(t) = b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3)), F(z) = 1 - exp(-exp(z)) and
#'     z = b3 (log t - log b2) + log((b3 - 1)/b3); b2 is the inflection point and b3 > 1 the shape;
#'   \item normal: f(t) = b1 Phi((t - b2)/b3), the standard normal distribution function Phi of
#'     z = (t - b2)/b3; symmetric around b2, where the curve reaches half of its asymptote;
#'   \item gamma: f(t) = b1 F(t), with F the cumulative distribution function (CDF) of the gamma
#'     distribution with shape b3 > 1 and scale b2/(b3 - 1), whose density has its mode at b2: b2 is the inflection
#'     point and b3 the shape;
#'   \item linear: f(t) = b1 t + b2, a straight line with slope b1 and intercept b2, for comparison
#'     with the growth curves. It has no asymptote, does not use b3-b5 and is not truncated;
#'     \code{\link{timeto}} gives the times at which it reaches given responses
#'     (\code{type="absolute"}), and \code{\link{maxrate}} its constant rate b1.
#' }
#' The asymmetry b5, the truncation at the threshold b4, the maximum growth rates and inflection
#' points, the rate 1/b3 and the small-sample inference are described in the vignette,
#' \code{vignette("gcTMB")}; the links in the argument \code{link}.
#'
#' @return A gcTMB object, with methods for \code{print}, \code{summary}, \code{coef}, \code{vcov},
#'   \code{predict}, \code{confint}, \code{anova}, \code{plot} and others, see \code{\link{gcTMB-methods}}.
#' @seealso \code{\link{maxrate}} for the maximum growth rate and its timing,
#'   \code{\link{timeto}} for the times to reach a level of the asymptote,
#'   \code{\link{gcTMBstart}} for the starting values.
#' @examples
#' # logistic growth of the faba bean plants
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25, model="logistic")
#' summary(fit)
#'
#' # the other models, compared by AIC
#' fitg <- gcTMB(height ~ DAS, data=FabaHeights25, model="gompertz")
#' fith <- gcTMB(height ~ DAS, data=FabaHeights25, model="hill")
#' fitw <- gcTMB(height ~ DAS, data=FabaHeights25, model="weibull")
#' fitn <- gcTMB(height ~ DAS, data=FabaHeights25, model="normal")
#' fitgam <- gcTMB(height ~ DAS, data=FabaHeights25, model="gamma")
#' fitlin <- gcTMB(height ~ DAS, data=FabaHeights25, model="linear")
#' AIC(fit, fitg, fith, fitw, fitn, fitgam, fitlin)
#'
#' # the Richards curve: the logistic model with an estimated asymmetry b5
#' fitr5 <- gcTMB(height ~ DAS, data=FabaHeights25, fix=c(b5=FALSE))
#' coef(fitr5)
#'
#' # the rate 1/b3 instead of the scale b3
#' gcTMB(height ~ DAS, data=FabaHeights25, link=c(b3="inverse"))
#'
#' # treatment effects on the asymptote as ratios, with the log link
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitl <- gcTMB(height ~ DAS, fform=b1 ~ trt, data=FabaHeights25, link=c(b1="log"))
#' exp(coef(fitl)[2:5])
#'
#' # treatment effects on the asymptote and random effects of the plants
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitre <- gcTMB(height ~ DAS, fform=b1 ~ trt, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' summary(fitre)
#'
#' # pod counts with a poisson distribution
#' gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson")
#'
#' # no pods before a threshold b4, after which the count grows
#' fitt <- gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE))
#' coef(fitt)
#' @export

gcTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, weights=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", truncated=TRUE, link=NULL, gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  # control list
  control <- do.call(gcTMBcontrol, control)
  dd <- gcTMBdata(form, fform=fform, rform=rform, data=data, weights=weights, family=family, model=model, truncated=truncated, link=link)
  fixnamed <- !is.null(names(fix))
  start <- parvector(start, dd$npi, NA, "start")
  fix <- defaultfix(fix, dd$npi)
  if (fixnamed && any(fix & is.na(start) & dd$npi < 4)){
    stop("give the values of the fixed parameters in start, e.g. start=c(b1=10), fix=c(b1=TRUE).")
  }
  fix[dd$npi %in% unusedpars(dd$mod, dd$trunc)] <- TRUE

  # starting values for parameters with NA in start
  method <- if (is.null(gsl) & is.null(gsu) & is.null(gsn)) "lm" else "grid"
  start <- gcTMBstartvals(dd, start, fix, method=method, gsl=gsl, gsu=gsu, gsn=gsn)

  obj <- gcTMBobj(dd, start, fix, rform=rform)
  # random effects: if the Laplace approximation fails at the starting values, e.g. starting values
  # far from the population curve, start from the fit without random effects
  if (!is.null(rform) && !is.finite(suppressWarnings(obj$fn(obj$par)))){
    fit0 <- tryCatch(suppressWarnings(gcTMB(form, fform=fform, rform=NULL, data=data, start=start, fix=fix, weights=weights,
                                            lower=lower, upper=upper, family=family, model=model, truncated=truncated,
                                            link=link, control=control)),
                     error=function(e) NULL)
    if (!is.null(fit0) && is.finite(fit0$opt$objective)){
      start[!fix] <- fit0$estimates[seq_len(sum(!fix))]
      obj <- gcTMBobj(dd, start, fix, rform=rform)
    }
  }

  # bounds apply to the free coefficients only, not to the variance parameters
  bounds <- function(b, inf){
    if (length(b) == 1) b <- rep(b, length(start))
    if (length(b) != length(start)) stop("lower and upper must be a single value or have the same length as start.")
    c(b[!fix], rep(inf, length(obj$par) - sum(!fix)))
  }
  lower <- bounds(lower, -Inf)
  upper <- bounds(upper, Inf)

  ctrl <- list(iter.max=control$iter.max, eval.max=control$eval.max)
  if (is.null(rform)){
    # Newton steps with the Hessian; quasi-Newton steps if the Hessian is not finite on the way, e.g.
    # for parameters far outside the range of the data
    opt <- tryCatch(nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper, control=ctrl),
                    error=function(e) nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl))
  } else {
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl)
    opt <- retryvariances(obj, opt, dd$re$terms, lower, upper, ctrl)
  }
  # truncated curves: the likelihood has a kink in the threshold b4 at each observation time, where
  # the optimiser often stops before the other parameters are optimal
  i4 <- which(dd$npi[!fix] == 4)
  if (isTRUE(dd$trunc) && length(i4) > 0 && is.finite(opt$objective) && opt$convergence != 0){
    opt <- thresholdopt(obj, opt, i4, lower, upper, ctrl, dd, start, fix)
  }
  if (!is.finite(opt$objective)){
    warning("the negative log-likelihood is not finite at the estimates; some observations may be impossible under the model, e.g. positive counts where the curve is zero.", call.=FALSE)
  } else if (opt$convergence != 0){
    warning("the optimiser did not converge: ", opt$message, call.=FALSE)
  }
  res <- list()
  res$call <- match.call()
  res$data <- dd$data
  res$formula <- form
  res$family <- family
  res$model <- model
  res$truncated <- dd$trunc
  res$link <- dd$link
  res$start <- start
  res$fix <- fix
  res$Xs <- dd$Xs
  res$fform <- dd$fform
  res$fterms <- dd$fterms
  res$coefnames <- coefnames(dd$Xs)
  res$rform <- rform
  res$reterms <- dd$re$terms
  res$estimates <- opt$par
  res$lower <- lower
  res$upper <- upper
  res$opt <- opt
  res$control <- control
  res$obj <- obj
  # without random effects, the Hessian by automatic differentiation: the default finite differences of
  # sdreport use absolute steps, which fail for coefficients on a small scale
  res$sdrl <- if (isFALSE(control$sdreport)) nosdreport(obj, opt, dd$re$terms) else
    sdreport(obj, par.fixed=opt$par, hessian.fixed=if (is.null(rform)) obj$he(opt$par))
  # NaN standard errors of a Hessian that is not positive definite are reported by summary.gcTMB
  res$ssdrl <- if (inherits(res$sdrl, "nosdreport")) nosdrsummary(res$sdrl) else suppressWarnings(summary(res$sdrl))
  # residual degrees of freedom for small-sample inference of gaussian and lognormal models without
  # random effects (lognormal: gaussian for log(y)): observations with positive weights minus free coefficients
  if (isTRUE(control$smallsample) && family %in% c("gaussian", "lognormal") && length(dd$re$terms) == 0){
    rdf <- sum(dd$w > 0) - sum(!fix)
    if (rdf > 0) res$rdf <- rdf else warning("no residual degrees of freedom; the standard errors use the maximum likelihood estimate of the residual variance.", call.=FALSE)
  }
  class(res) <- "gcTMB"
  return(res)
}


# optimum of a truncated curve with the threshold b4 at a kink of the likelihood: alternately the other
# parameters with the coefficients i4 of b4 held, which is smooth, and all parameters, until the
# negative log-likelihood no longer decreases. With b4 at an observation time of its curve and the
# other parameters converged, the kink is the optimum in b4 and the fit is reported as converged.
thresholdopt <- function(obj, opt, i4, lower, upper, ctrl, dd, start, fix){
  best <- opt
  conditional <- FALSE
  for (k in 1:5){
    p <- best$par
    fn <- function(q){ p[-i4] <- q; obj$fn(p) }
    gr <- function(q){ p[-i4] <- q; as.vector(obj$gr(p))[-i4] }
    c1 <- tryCatch(nlminb(p[-i4], fn, gr, lower=lower[-i4], upper=upper[-i4], control=ctrl), error=function(e) NULL)
    improved <- FALSE
    if (!is.null(c1) && is.finite(c1$objective) && c1$objective < best$objective - 1e-10){
      p[-i4] <- c1$par
      best$par <- p
      best$objective <- c1$objective
      best$convergence <- c1$convergence
      best$message <- c1$message
      conditional <- TRUE
      improved <- TRUE
    }
    c2 <- tryCatch(nlminb(best$par, obj$fn, obj$gr, lower=lower, upper=upper, control=ctrl), error=function(e) NULL)
    if (!is.null(c2) && is.finite(c2$objective) && c2$objective < best$objective - 1e-10){
      best <- c2
      conditional <- FALSE
      improved <- TRUE
    }
    if (!improved || best$convergence == 0) break
  }
  obj$fn(best$par)
  # the threshold of each observation's curve at one of the observation times
  pv <- start
  pv[!fix] <- best$par[seq_len(sum(!fix))]
  b4 <- linkfun(dd$lnk[4], as.vector(dd$Xs[[4]] %*% pv[dd$npi == 4]), inverse=TRUE)
  corner <- any(abs(b4 - dd$x) <= 1e-6 * pmax(1, abs(dd$x)))
  if (conditional && best$convergence == 0 && corner){
    best$message <- "converged, with the threshold b4 at an observation time (a kink of the likelihood)"
  } else if (conditional){
    best$convergence <- 1L
  }
  best
}


# the parts of sdreport used by the methods, without the covariance matrix of the estimates
# (control sdreport=FALSE): the estimates, the random effects at their conditional modes, and the
# correlation matrices of the random effects (as the ADREPORT corr), with NA standard errors;
# nosdrsummary gives the same table as summary of an sdreport
nosdreport <- function(obj, opt, terms){
  env <- obj$env
  r <- env$random
  par.random <- if (!is.null(r)) env$last.par.best[r] else numeric(0)
  pl <- env$parList(opt$par)
  toff <- 0
  corr <- unlist(lapply(terms, function(tm){
    np <- length(tm$pars)
    cr <- corrmatrix(pl$theta[toff + np + seq_len(np*(np - 1)/2)], np)
    toff <<- toff + np + np*(np - 1)/2
    as.vector(cr)
  }))
  np <- length(opt$par)
  structure(list(par.fixed=opt$par, par.random=par.random, cov.fixed=matrix(NA_real_, np, np),
                 pdHess=NA, corr=corr), class="nosdreport")
}

nosdrsummary <- function(object){
  est <- c(object$corr, object$par.fixed, object$par.random)
  nm <- c(rep("corr", length(object$corr)), names(object$par.fixed), rep("u", length(object$par.random)))
  matrix(c(est, rep(NA_real_, length(est))), ncol=2, dimnames=list(nm, c("Estimate", "Std. Error")))
}
