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
#' @param fix Logical vector or NULL. One entry for each model coefficient. If TRUE, a model coefficient is fixed to its value in start. Alternatively, a vector named by parameter, e.g. c(b1=TRUE), fixes all coefficients of the named parameters; their values must then be given in start, e.g. start=c(b1=10). The reserved parameters b4 and b5, which the models do not use, are always fixed, at 0 unless given in start.
#' @param weights Numeric vector or NULL. Weights, one for each row of data.
#' @param lower Numeric, a single value or a vector. Lower bounds for the model coefficients b1-b5, a single value or a vector with the same length as start. Bounds of fixed coefficients are ignored; variance parameters are unbounded.
#' @param upper Numeric, a single value or a vector. Upper bounds for the model coefficients b1-b5, a single value or a vector with the same length as start.
#' @param family Character string. The distribution of the response, one of "gaussian", "lognormal", "binomial", "betabinomial", "beta", "poisson", "nbinom2", and "nbinom1". For the lognormal family, log(y) is normal with mean log f(t) and standard deviation sigma: the curve f is the median of a positive response, with a constant coefficient of variation; the mean is f(t) exp(sigma^2/2). The negative binomial families are for overdispersed counts, with variance mu (1 + mu/theta) for nbinom2 and mu (1 + phi) for nbinom1, and the beta-binomial family for overdispersed binomial counts, with variance n mu (1 - mu) (1 + (n - 1)/(phi + 1)) for n trials, as in glmmTMB.
#' @param model Character string. The growth model: "logistic" (default) or "gompertz", see Details.
#' @param link Named character vector or NULL. Links of the parameters b1-b5: NULL for the defaults (the identity link for all parameters), or links named by parameter. The scale b3 has the link "identity" (default) or "inverse", \code{link=c(b3="inverse")}, with which its coefficients are the growth rate 1/b3; the other parameters have the identity link. The coefficients, and start, fix, lower, upper and linfct of \code{\link{timeto}}, are on the scale of the links, see Details.
#' @param gsl Numeric vector or NULL. Lower bounds for a grid search, same length as start; NA for the default range. If gsl, gsu and gsn are all NULL, no grid search is performed. See \code{\link{gcTMBstart}}.
#' @param gsu Numeric vector or NULL. Upper bounds for a grid search, same length as start; NA for the default range.
#' @param gsn Integer, a single number or a vector, or NULL. Number of grid points for each parameter, a single number or a vector with the same length as start. \code{gsn=10} alone runs a grid search over the default ranges.
#' @param control List. Control parameters, the arguments of \code{\link{gcTMBcontrol}}, e.g. \code{list(iter.max=2000)}. See \code{\link{gcTMBcontrol}} for details.
#'
#' @details The growth models are of the form f(t) = b1 F((t - b2)/b3), with a function F increasing
#' from 0 to 1, which has its largest slope at 0:
#' \itemize{
#'   \item logistic: F(z) = 1/(1 + exp(-z)), f(t) = b1 / (1 + exp(-(t - b2)/b3)), symmetric around b2,
#'     where the curve reaches half of its asymptote;
#'   \item gompertz: F(z) = exp(-exp(-z)), f(t) = b1 exp(-exp(-(t - b2)/b3)), which reaches
#'     exp(-1), about 37\%, of its asymptote at b2 and approaches the asymptote more slowly than it
#'     leaves zero.
#' }
#' The curves start at zero and approach the asymptote b1 at late times. b2 is the time of the
#' maximum growth rate, the inflection point of the curve, and b3 the scale of the curve on the time
#' axis, positive for an increasing curve: the maximum growth rate is b1/(4 b3) for the logistic and
#' b1/(e b3) for the gompertz model. The parameters b4 and b5 are reserved for future models;
#' they are not used, and are always fixed.
#'
#' With the inverse link for b3, \code{link=c(b3="inverse")}, the linear predictor of b3, including
#' its random effects, is the rate 1/b3: the coefficients of b3 in \code{coef}, \code{start},
#' \code{fix}, \code{lower}, \code{upper} and \code{linfct} of \code{\link{timeto}} are rates per unit
#' of time, and covariates and random effects act additively on the rate. Without covariates and
#' random effects on b3, both links give the same fit.
#'
#' Gaussian models without random effects are estimated by maximum likelihood, but by default
#' (\code{smallsample=TRUE} in \code{\link{gcTMBcontrol}}) their inference follows lm and nls:
#' the residual standard deviation is the square root of the unbiased residual variance, the weighted
#' residual sum of squares divided by n - p residual degrees of freedom, for n observations with
#' positive weights and p free coefficients b1-b5; the covariance of the coefficients
#' (\code{vcov}, \code{summary}, \code{\link{timeto}}, \code{predict}) is scaled accordingly, the Wald
#' and prediction intervals use t quantiles with n - p degrees of freedom, and the profile intervals of
#' the coefficients and times an F-based cutoff. The log-likelihood, the estimate log_sigma in
#' \code{coef}, AIC, \code{anova}, \code{simulate} and \code{\link{bootstrap}} remain maximum likelihood.
#' The other families and models with random effects use maximum likelihood and normal quantiles.
#'
#' @return A gcTMB object, with methods for \code{print}, \code{summary}, \code{coef}, \code{vcov},
#'   \code{predict}, \code{confint}, \code{anova}, \code{plot} and others, see \code{\link{gcTMB-methods}}.
#' @seealso \code{\link{timeto}} for the times to reach a level of the asymptote,
#'   \code{\link{gcTMBstart}} for the starting values.
#' @examples
#' # logistic growth of the faba bean plants
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25, model="logistic")
#' summary(fit)
#'
#' # Gompertz model, compared by AIC
#' fitg <- gcTMB(height ~ DAS, data=FabaHeights25, model="gompertz")
#' AIC(fit, fitg)
#'
#' # the rate 1/b3 instead of the scale b3
#' gcTMB(height ~ DAS, data=FabaHeights25, link=c(b3="inverse"))
#'
#' # treatment effects on the asymptote and random effects of the plants
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitre <- gcTMB(height ~ DAS, fform=b1 ~ trt, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' summary(fitre)
#'
#' # pod counts with a poisson distribution
#' gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson")
#' @export

gcTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, weights=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", link=NULL, gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  # control list
  control <- do.call(gcTMBcontrol, control)
  dd <- gcTMBdata(form, fform=fform, rform=rform, data=data, weights=weights, family=family, model=model, link=link)
  fixnamed <- !is.null(names(fix))
  start <- parvector(start, dd$npi, NA, "start")
  fix <- parvector(fix, dd$npi, FALSE, "fix")
  if (fixnamed && any(fix & is.na(start))){
    stop("give the values of the fixed parameters in start, e.g. start=c(b1=10), fix=c(b1=TRUE).")
  }
  fix[dd$npi %in% unusedpars(dd$mod)] <- TRUE

  # starting values for parameters with NA in start
  method <- if (is.null(gsl) & is.null(gsu) & is.null(gsn)) "lm" else "grid"
  start <- gcTMBstartvals(dd, start, fix, method=method, gsl=gsl, gsu=gsu, gsn=gsn)

  obj <- gcTMBobj(dd, start, fix, rform=rform)

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
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper, control=ctrl)
  } else {
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl)
    opt <- retryvariances(obj, opt, dd$re$terms, lower, upper, ctrl)
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
  res$sdrl <- sdreport(obj, par.fixed=opt$par, hessian.fixed=if (is.null(rform)) obj$he(opt$par))
  # NaN standard errors of a Hessian that is not positive definite are reported by summary.gcTMB
  res$ssdrl <- suppressWarnings(summary(res$sdrl))
  # residual degrees of freedom for small-sample inference of gaussian and lognormal models without
  # random effects (lognormal: gaussian for log(y)): observations with positive weights minus free coefficients
  if (isTRUE(control$smallsample) && family %in% c("gaussian", "lognormal") && length(dd$re$terms) == 0){
    rdf <- sum(dd$w > 0) - sum(!fix)
    if (rdf > 0) res$rdf <- rdf else warning("no residual degrees of freedom; the standard errors use the maximum likelihood estimate of the residual variance.", call.=FALSE)
  }
  class(res) <- "gcTMB"
  return(res)
}
