#' Variance functions for heteroscedastic residuals
#'
#' Variance functions for the \code{variance} argument of \code{\link{gcTMB}}, with the arguments and
#' the parameterisation of \code{varIdent}, \code{varPower}, \code{varConstPower}, \code{varExp} and
#' \code{varComb} of the package nlme (whose objects \code{\link{gcTMB}} also accepts), and the
#' polynomial \code{varPoly}, for the gaussian and lognormal (on the log scale) families. The residual
#' variance of an observation is sigma^2 g^2, with the residual standard deviation sigma of the model and
#' the function g:
#' \itemize{
#'   \item \code{varIdent(form=~ 1 | s)}: g = delta_s, a ratio for each stratum s of the grouping
#'     variables, with delta = 1 for the reference stratum, so that each stratum has its own variance,
#'     e.g. a variance for each treatment (\code{~ 1 | Treatment}) or for each time point
#'     (\code{~ 1 | Time}); a numeric grouping variable is treated as a factor;
#'   \item \code{varPower(form=~ v)}: g = |v|^delta, variance sigma^2 |v|^(2 delta);
#'   \item \code{varConstPower(form=~ v)}: g = c + |v|^delta, c > 0, variance sigma^2 (c + |v|^delta)^2;
#'   \item \code{varExp(form=~ v)}: g = exp(t v), variance sigma^2 exp(2 t v);
#'   \item \code{varPoly(form=~ v, degree=k)}: g = 1 + d_1 (v - m) + d_2 (v - m)^2 + ... + d_k (v - m)^k,
#'     so that the standard deviation sigma g is a polynomial in v, e.g. one that is largest in the
#'     middle of the range of v (d_1 > 0, d_2 < 0, with its maximum at v = m - d_1/(2 d_2) for degree 2).
#'     The polynomial is centred at the mean m of v in the data (the mean response for \code{fitted(.)},
#'     m = 0 for \code{rate(.)}), shown by print, so that sigma is the standard deviation at v = m, and
#'     the standard deviation at v = 0 can take any value; an uncentred polynomial, with sigma the
#'     standard deviation at v = 0, could not describe a standard deviation that would be 0 or below
#'     at v = 0 but is positive on the data, e.g. a straight line in the time from the first day. A
#'     polynomial can be 0 or negative, where the standard deviation is not defined: the likelihood is zero
#'     (the log-likelihood minus infinity) if the polynomial is 0 or below at any observation, so that the
#'     optimiser stays where it is positive on the data. It starts from d = 0, a constant standard
#'     deviation. Outside the data, e.g. in \code{predict} and \code{simulate} at new times, the polynomial
#'     can still be 0 or below, which gives an error. A degree above 2 is rarely needed and tends to overfit;
#'   \item \code{varComb(...)}: the product of the multipliers g of two or more of the above, e.g.
#'     \code{varComb(varIdent(form=~ 1 | Treatment), varPower())} for a power of the fitted curve with a
#'     different factor for each treatment, or \code{varComb(varPower(), varExp())}, g = |v|^delta exp(t v),
#'     which with delta > 0 and t < 0 is largest at v = -delta/t, in the middle of growth. The variances are
#'     multiplied, not added: the standard deviations are sigma g_A g_B.
#' }
#' The covariate v of \code{varPower} and \code{varConstPower} is a numeric variable of the data,
#' e.g. \code{~ Time}, \code{fitted(.)} (the default), the fitted curve f(t) of each observation
#' (with its random effects), or \code{rate(.)}, the growth rate f'(t) of the fitted curve at the
#' time of the observation, its slope, which is largest at the inflection point. \code{rate(.)} is
#' a gcTMB addition, not available in nlme; it is computed by a central difference. With a grouping
#' factor, \code{form=~ v | g}, each level of g has its own power (and constant).
#'
#' The variance of \code{varPower} goes to zero where |v| does, e.g. the fitted values of a curve
#' that starts at zero or the growth rate before and after growth, which makes those observations
#' too influential and the fit unstable; \code{varConstPower} keeps the variance at least sigma^2
#' c^2 there. A negative power gives an infinite variance at |v| = 0; with \code{rate(.)}, this is the
#' case for curves whose growth rate passes through zero within the data, such as the beta curve after
#' its maximum, which the optimiser cannot handle, so \code{rate(.)} is meant for curves that keep
#' growing. Unlike nlme, which updates the fitted values between fits of the variance parameters,
#' gcTMB computes them inside the likelihood, so the fit is the joint maximum likelihood estimate.
#' The reference stratum of \code{varIdent} is the first level of the grouping variables (sorted
#' numerically for a numeric variable), unless \code{value} and \code{fixed} name all the others;
#' nlme takes the first in the order of the data.
#'
#' A variance function and a \code{correlation} structure can be combined: the correlation is then that
#' of the standardised residuals e/(sigma g). With \code{varIdent(~ 1 | Time)}, \code{simulate} and
#' prediction intervals need the time points of the fit.
#'
#' @param value Numeric. Starting values (or with \code{fixed}, see there): for \code{varIdent} the
#'   ratios of the standard deviations of the strata to the reference stratum, named by stratum (1 by
#'   default); for \code{varPower} the power and for \code{varExp} the exponent t (0 by default), named by
#'   level with a grouping factor; for \code{varPoly} the coefficients d_1, ..., d_k (0 by default), one for
#'   each power of v, or named \code{poly1}, ..., \code{polyk} (\code{poly1:<level>} for a single level
#'   of a grouping factor)
#' @param const,power Numeric. For \code{varConstPower}, the starting values of the constant c (0.1 by
#'   default, as nlme) and the power delta (0 by default)
#' @param degree Integer. The degree k of the polynomial of \code{varPoly}, 2 by default
#' @param form One-sided formula: \code{~ 1 | s} for \code{varIdent}, with s one or more grouping
#'   variables (\code{~ 1 | A * B}); \code{~ v} or \code{~ v | g} for \code{varPower},
#'   \code{varConstPower}, \code{varExp} and \code{varPoly}, with v a numeric variable, \code{fitted(.)}
#'   or \code{rate(.)}
#' @param fixed Values that are fixed and not estimated: for \code{varIdent} ratios named by stratum,
#'   for \code{varPower} the power and for \code{varExp} the exponent (named by level with a grouping
#'   factor), for \code{varConstPower} a list with \code{const} and/or \code{power}, e.g.
#'   \code{fixed=list(power=1)}, for \code{varPoly} coefficients named as in \code{value}, e.g.
#'   \code{fixed=c(poly2=0)}
#'
#' @details The estimated parameters are in \code{coef(fit)} as \code{var_par:<stratum>} (the log
#' ratios of \code{varIdent}), \code{var_par:power}, \code{var_par:const} (the log of c),
#' \code{var_par:expon} and \code{var_par:poly1}, \code{var_par:poly2}, ..., with \code{:<level>} for a
#' grouping factor, on the unconstrained scale of nlme; \code{\link{varpars}} gives them on the natural
#' scale. The parameters of the components of \code{varComb} are prefixed with \code{A.}, \code{B.}, ...,
#' as in nlme.
#' @return An object of class \code{c("varIdent", "varFunc")}, \code{c("varPower", "varFunc")},
#'   \code{c("varConstPower", "varFunc")}, \code{c("varExp", "varFunc")}, \code{c("varPoly", "varFunc")}
#'   or \code{c("varComb", "varFunc")}, with the structure of the nlme objects (a list of the components
#'   for \code{varComb}).
#' @references Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS.
#'   Springer, New York.
#' @seealso \code{\link{gcTMB}}, \code{\link{varpars}}, \code{\link{corStruct}}
#' @examples
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' # a residual variance for each treatment
#' fiti <- gcTMB(Height ~ DAS, fform=b1 ~ trt, data=FabaHeights25, variance=varIdent(form=~ 1 | trt))
#' varpars(fiti)
#' # residual standard deviation c + f'(t)^delta, largest where the curve grows fastest
#' fitr <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varConstPower(form=~ rate(.)))
#' varpars(fitr)
#' \donttest{
#' # a standard deviation that is a quadratic polynomial of the fitted curve
#' fitp <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varPoly(degree=2))
#' varpars(fitp)
#' # the product of a power and an exponential of the fitted curve
#' fitc <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varComb(varPower(), varExp()))
#' varpars(fitc)
#' }
#' @name varFunc
#' @aliases varIdent varPower varConstPower varExp varPoly varComb
NULL

#' @rdname varFunc
#' @export
varIdent <- function(value=numeric(0), form=~ 1, fixed=NULL){
  vf <- varform(form)
  if (is.null(vf$group)) stop("varIdent needs a grouping factor, e.g. form=~ 1 | Treatment.")
  value <- unlist(value)
  if (length(value) > 0){
    if (is.null(names(value)) && length(value) > 1) stop("initial values must have group names in 'varIdent'")
    if (any(value <= 0)) stop("initial values for 'varIdent' must be > 0")
    value <- log(value)
  }
  grp <- names(value)
  if (!is.null(fixed)){
    fixed <- unlist(fixed)
    if (is.null(names(fixed))) stop("fixed parameters must have names in 'varIdent'")
    if (any(fixed <= 0)) stop("fixed values for 'varIdent' must be > 0")
    attr(value, "fixed") <- log(fixed)
    grp <- c(grp, names(fixed))
  }
  attr(value, "groupNames") <- grp
  attr(value, "formula") <- form
  class(value) <- c("varIdent", "varFunc")
  value
}

#' @rdname varFunc
#' @export
varPower <- function(value=numeric(0), form=~ fitted(.), fixed=NULL){
  vf <- varform(form)
  if (is.null(vf$cov)) stop("'form' must have a covariate")
  value <- unlist(value)
  fixed <- unlist(fixed)
  attr(value, "fixed") <- fixed
  attr(value, "formula") <- form
  if (!is.null(vf$group)){
    if (is.null(names(value)) && length(value) > 1) stop("initial values must have group names in 'varPower'")
    if (!is.null(fixed) && is.null(names(fixed))) stop("fixed parameters must have group names in 'varPower'")
    attr(value, "groupNames") <- c(names(value), names(fixed))
  } else {
    attr(value, "whichFix") <- !is.null(fixed)
  }
  class(value) <- c("varPower", "varFunc")
  value
}

#' @rdname varFunc
#' @export
varConstPower <- function(const=numeric(0), power=numeric(0), form=~ fitted(.), fixed=NULL){
  vf <- varform(form)
  if (is.null(vf$cov)) stop("'form' must have a covariate")
  if (!is.null(fixed) && (!is.list(fixed) || is.null(names(fixed)) || any(!(names(fixed) %in% c("const", "power"))))){
    stop("fixed must be a list with the components const and/or power, e.g. fixed=list(power=1).")
  }
  logc <- function(v){
    v <- unlist(v)
    if (length(v) > 0 && any(v <= 0)) stop("constant in \"varConstPower\" structure must be > 0")
    if (length(v) > 0) log(v) else v
  }
  value <- list(const=logc(const), power=unlist(power))
  if (!is.null(fixed$const)) attr(value$const, "fixed") <- logc(fixed$const)
  if (!is.null(fixed$power)) attr(value$power, "fixed") <- unlist(fixed$power)
  attr(value, "formula") <- form
  if (is.null(vf$group)){
    wf <- array(FALSE, c(2, 1), list(c("const", "power"), NULL))
    wf[, 1] <- c(!is.null(fixed$const), !is.null(fixed$power))
    attr(value, "whichFix") <- wf
  }
  class(value) <- c("varConstPower", "varFunc")
  value
}

#' @rdname varFunc
#' @export
varExp <- function(value=numeric(0), form=~ fitted(.), fixed=NULL){
  vf <- varform(form)
  if (is.null(vf$cov)) stop("'form' must have a covariate")
  value <- unlist(value)
  fixed <- unlist(fixed)
  attr(value, "fixed") <- fixed
  attr(value, "formula") <- form
  if (!is.null(vf$group)){
    if (is.null(names(value)) && length(value) > 1) stop("initial values must have group names in 'varExp'")
    if (!is.null(fixed) && is.null(names(fixed))) stop("fixed parameters must have group names in 'varExp'")
    attr(value, "groupNames") <- c(names(value), names(fixed))
  } else {
    attr(value, "whichFix") <- !is.null(fixed)
  }
  class(value) <- c("varExp", "varFunc")
  value
}

#' @rdname varFunc
#' @export
varPoly <- function(value=numeric(0), degree=2, form=~ fitted(.), fixed=NULL){
  vf <- varform(form)
  if (is.null(vf$cov)) stop("'form' must have a covariate")
  if (length(degree) != 1 || !is.numeric(degree) || !is.finite(degree) || degree < 1 || degree != round(degree)) stop("degree must be a positive whole number.")
  degree <- as.integer(degree)
  check <- function(v, what){
    v <- unlist(v)
    if (length(v) == 0) return(numeric(0))
    if (!is.numeric(v) || any(!is.finite(v))) stop("the ", what, " of varPoly must be finite numbers.")
    if (is.null(names(v))){
      if (length(v) != degree) stop("the ", what, " of varPoly without names must have one value for each power of the covariate, ", degree, " here.")
      names(v) <- paste("poly", seq_len(degree), sep="")
    }
    v
  }
  value <- check(value, "initial values")
  fixed <- check(fixed, "fixed values")
  attr(value, "fixed") <- if (length(fixed) > 0) fixed
  attr(value, "formula") <- form
  attr(value, "degree") <- degree
  class(value) <- c("varPoly", "varFunc")
  value
}

#' @rdname varFunc
#' @param ... For \code{varComb}, two or more variance functions, whose standard deviation multipliers
#'   g are multiplied
#' @export
varComb <- function(...){
  vl <- list(...)
  vl <- unlist(lapply(vl, function(v) if (inherits(v, "varComb")) unname(unclass(v)) else list(v)), recursive=FALSE)
  if (length(vl) < 2) stop("varComb needs at least two variance functions.")
  for (v in vl){
    if (!inherits(v, "varFunc") || vartypeone(v) == 0) stop("the arguments of varComb must be varIdent(), varPower(), varConstPower(), varExp() or varPoly().")
  }
  names(vl) <- LETTERS[seq_along(vl)]
  class(vl) <- c("varComb", "varFunc")
  vl
}

# type of a single variance function: 1 varIdent, 2 varPower, 3 varConstPower, 4 varExp, 5 varPoly,
# 0 otherwise
vartypeone <- function(v){
  if (inherits(v, "varIdent")) return(1L)
  if (inherits(v, "varPower")) return(2L)
  if (inherits(v, "varConstPower")) return(3L)
  if (inherits(v, "varExp")) return(4L)
  if (inherits(v, "varPoly")) return(5L)
  0L
}

# type of a variance function: 0 none, 1-5 a single one (see vartypeone), 6 varComb
vartypeindex <- function(variance){
  if (is.null(variance)) return(0L)
  if (inherits(variance, "varComb")) return(6L)
  vt <- vartypeone(variance)
  if (vt == 0) stop("variance must be varIdent(), varPower(), varConstPower(), varExp(), varPoly() or varComb().")
  vt
}

# the single variance functions of a variance function: its components for varComb, itself otherwise,
# none for NULL
varcomps <- function(variance){
  if (is.null(variance)) return(list())
  if (inherits(variance, "varComb")) return(unname(unclass(variance)))
  list(variance)
}

# prefixes of the parameter names of the components: A., B., ... for varComb, none otherwise
varprefix <- function(variance){
  comps <- varcomps(variance)
  if (inherits(variance, "varComb")) paste(LETTERS[seq_along(comps)], ".", sep="") else rep("", length(comps))
}

# degree of the polynomial of varPoly, 0 for the other variance functions
vardegree <- function(v){
  if (inherits(v, "varPoly")) as.integer(attr(v, "degree")) else 0L
}

# number of parameters of a single variance function of type vt with G levels (the strata of varIdent,
# with the reference) and degree deg
varnpar <- function(vt, G, deg){
  switch(vt, G - 1L, G, 2L*G, G, deg*G)
}

# parts of the form ~ v | g of a variance function: the covariate expression (NULL for ~ 1) and the
# grouping expression (NULL without a grouping factor)
varform <- function(form){
  if (is.null(form)) form <- ~ 1
  if (!inherits(form, "formula") || length(form) != 2) stop("the form of a variance function must be a one-sided formula, e.g. ~ fitted(.) or ~ 1 | group.")
  rhs <- form[[2]]
  grp <- NULL
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))){
    grp <- rhs[[3]]
    rhs <- rhs[[2]]
  }
  cv <- if (identical(rhs, 1) || identical(rhs, 1L)) NULL else rhs
  list(cov=cv, group=grp)
}

# covariate type: 0 a variable of the data, 1 fitted(.), 2 rate(.)
varcovtype <- function(cv){
  if (is.null(cv)) return(0L)
  if (identical(cv, quote(fitted(.)))) return(1L)
  if (identical(cv, quote(rate(.)))) return(2L)
  0L
}

# variables of a variance function, for the missing values of the data
varvars <- function(variance){
  v <- unlist(lapply(varcomps(variance), function(vc){
    vf <- varform(attr(vc, "formula"))
    c(if (varcovtype(vf$cov) == 0) all.vars(vf$cov), all.vars(vf$group))
  }))
  setdiff(unique(c(character(0), v)), ".")
}

# strata of the grouping expression for the rows of data: the labels, and the levels in their order
# (factor levels, numeric values sorted numerically), interactions labelled a*b as in nlme
varstrata <- function(grp, data){
  vars <- all.vars(grp)
  miss <- setdiff(vars, names(data))
  if (length(miss) > 0) stop("the data need the grouping variables of the variance function: ", paste(miss, collapse=", "), ".")
  levs <- lapply(vars, function(v){
    x <- data[[v]]
    if (is.factor(x)) levels(droplevels(x)) else as.character(sort(unique(x)))
  })
  lab <- do.call(paste, c(lapply(vars, function(v) as.character(data[[v]])), sep="*"))
  all <- do.call(paste, c(rev(expand.grid(rev(levs), stringsAsFactors=FALSE)), sep="*"))
  list(label=lab, levels=all[all %in% lab])
}

# starting values of the parameters of each level: the named values (or an unnamed single value for
# all levels), and the default otherwise; the fixed ones, NA for free
varlevelvals <- function(val, fix, levels, default, what){
  st <- setNames(rep(default, length(levels)), levels)
  fx <- setNames(rep(NA_real_, length(levels)), levels)
  setvals <- function(target, v){
    if (length(v) == 0) return(target)
    if (is.null(names(v))){
      if (length(v) != 1) stop("the ", what, " values of the variance function must be named by the levels: ", paste(levels, collapse=", "), ".")
      target[] <- v
    } else {
      bad <- setdiff(names(v), levels)
      if (length(bad) > 0) stop("unknown level(s) ", paste(bad, collapse=", "), " in the ", what, " values of the variance function; the levels are ", paste(levels, collapse=", "), ".")
      target[names(v)] <- v
    }
    target
  }
  st <- setvals(st, setNames(as.numeric(val), names(val)))
  fx <- setvals(fx, if (is.null(fix)) numeric(0) else setNames(as.numeric(fix), names(fix)))
  st[!is.na(fx)] <- fx[!is.na(fx)]
  list(start=st, fixed=fx)
}

# the data of a variance function for the template, its components side by side (see the template):
# the types vartype, the covariate types varcovt, the covariates varcov (n x K), the strata (varIdent)
# or levels (the others) varstr (n x K, zero-based, the reference stratum 0), the numbers of levels
# varng, the degrees vardeg (varPoly) and the offsets varoff of the parameters in var_par; and the
# parameters var_par (on the unconstrained scale), their map (NA for fixed), names and the levels of
# each component
varsetup <- function(variance, data, y=NULL, x=NULL){
  n <- nrow(data)
  out <- list(vartype=integer(0), varcovt=integer(0), varcov=matrix(0, n, 0), varstr=matrix(0L, n, 0),
              varng=integer(0), vardeg=integer(0), varoff=integer(0), varscale=numeric(0), varcenter=numeric(0),
              var_par=numeric(0), map=integer(0), names=character(0), levels=list())
  comps <- varcomps(variance)
  if (length(comps) == 0) return(out)
  pre <- varprefix(variance)
  off <- 0L
  for (k in seq_along(comps)){
    s1 <- varsetupone(comps[[k]], data, y, x)
    out$vartype <- c(out$vartype, s1$vartype)
    out$varcovt <- c(out$varcovt, s1$varcovt)
    out$varcov <- cbind(out$varcov, s1$varcov)
    out$varstr <- cbind(out$varstr, s1$varstr)
    out$varng <- c(out$varng, s1$varng)
    out$vardeg <- c(out$vardeg, s1$vardeg)
    out$varoff <- c(out$varoff, off)
    out$varscale <- c(out$varscale, s1$varscale)
    out$varcenter <- c(out$varcenter, s1$varcenter)
    out$var_par <- c(out$var_par, s1$var_par)
    out$map <- c(out$map, s1$map)
    out$names <- c(out$names, paste("var_par:", pre[k], s1$names, sep=""))
    out$levels[[k]] <- s1$levels
    off <- off + length(s1$var_par)
  }
  storage.mode(out$varstr) <- "integer"
  storage.mode(out$vartype) <- storage.mode(out$varcovt) <- storage.mode(out$varng) <- storage.mode(out$vardeg) <- storage.mode(out$varoff) <- "integer"
  free <- !is.na(out$map)
  out$map[free] <- seq_len(sum(free))
  out
}

# the data of a single variance function (see varsetup), with the names of its parameters without the
# prefix var_par:. The polynomial of varPoly is in v - m, centred at the mean m of v in the data (the
# mean response for fitted(.), 0 for rate(.)), so that sigma is the standard deviation at v = m; inside
# the template it is in (v - m)/s, with the scale s of v - m in the data (the largest |y - m| for
# fitted(.), the largest |y| divided by the range of x for rate(.), the largest |v - m| for a variable;
# s = 1 and m = 0 without the response y), so that its coefficients are of a similar size. The centre
# and scale are the attributes "center" and "scale" of its levels
varsetupone <- function(variance, data, y=NULL, x=NULL){
  n <- nrow(data)
  vt <- vartypeone(variance)
  out <- list(vartype=vt, varcovt=0L, varcov=rep(0, n), varstr=rep(0L, n), varng=1L, vardeg=vardegree(variance), varscale=1, varcenter=0,
              var_par=numeric(0), map=integer(0), names=character(0), levels=NULL)
  form <- attr(variance, "formula")
  vf <- varform(form)
  if (vt == 1){
    if (is.null(vf$group)) stop("varIdent needs a grouping factor, e.g. form=~ 1 | Treatment.")
    st <- varstrata(vf$group, data)
    lev <- st$levels
    if (length(lev) < 2) stop("the grouping variables of varIdent have a single stratum in the data.")
    val <- as.numeric(variance)
    names(val) <- names(variance)
    fix <- attr(variance, "fixed")
    named <- unique(c(names(val), names(fix)))
    # the reference stratum: the one not named, if all others are, otherwise the first level
    ref <- if (length(named) == length(lev) - 1 && all(named %in% lev)) setdiff(lev, named) else lev[1]
    if (ref %in% named) stop("the reference stratum ", ref, " of varIdent cannot have a value; it is 1.")
    others <- setdiff(lev, ref)
    if (length(val) == 1 && is.null(names(val))){
      if (length(others) > 1) stop("initial values must have group names in 'varIdent'")
      names(val) <- others
    }
    lv <- varlevelvals(val, fix, others, 0, "varIdent")
    out$varstr <- as.integer(match(st$label, c(ref, others)) - 1L)
    out$varng <- length(lev)
    out$var_par <- unname(lv$start)
    out$map <- ifelse(is.na(lv$fixed), seq_along(others), NA)
    out$names <- others
    out$levels <- c(ref, others)
    return(out)
  }
  ct <- varcovtype(vf$cov)
  out$varcovt <- ct
  if (ct == 0){
    v <- eval(vf$cov, data, environment(form))
    if (!is.numeric(v) || length(v) != n) stop("the covariate of the variance function must be a numeric variable of the data, fitted(.) or rate(.).")
    out$varcov <- as.numeric(v)
  }
  lev <- "1"
  if (!is.null(vf$group)){
    st <- varstrata(vf$group, data)
    lev <- st$levels
    out$varstr <- as.integer(match(st$label, lev) - 1L)
  }
  G <- length(lev)
  out$varng <- G
  sfx <- if (is.null(vf$group)) "" else paste(":", lev, sep="")
  if (vt %in% c(2, 4)){
    nm <- if (vt == 2) "power" else "expon"
    lv <- varlevelvals(setNames(as.numeric(variance), names(variance)), attr(variance, "fixed"), lev, 0, nm)
    out$var_par <- unname(lv$start)
    out$map <- ifelse(is.na(lv$fixed), seq_len(G), NA)
    out$names <- paste(nm, sfx, sep="")
  } else if (vt == 3){
    vc <- varlevelvals(variance$const, attr(variance$const, "fixed"), lev, log(0.1), "const")
    vp <- varlevelvals(variance$power, attr(variance$power, "fixed"), lev, 0, "power")
    out$var_par <- unname(c(vc$start, vp$start))
    out$map <- ifelse(is.na(c(vc$fixed, vp$fixed)), seq_len(2*G), NA)
    out$names <- c(paste("const", sfx, sep=""), paste("power", sfx, sep=""))
  } else {
    deg <- out$vardeg
    ce <- if (is.null(y) || ct == 2) 0 else if (ct == 1) mean(y) else mean(out$varcov)
    sc <- if (is.null(y)) 1 else if (ct == 1) max(abs(y - ce)) else if (ct == 2) max(abs(y))/diff(range(x)) else max(abs(out$varcov - ce))
    if (!is.finite(ce)) ce <- 0
    if (!is.finite(sc) || sc <= 0) sc <- 1
    out$varscale <- sc
    out$varcenter <- ce
    st <- varpolyvals(setNames(as.numeric(variance), names(variance)), matrix(0, G, deg), lev, deg)
    fx <- varpolyvals(attr(variance, "fixed"), matrix(NA_real_, G, deg), lev, deg)
    st[!is.na(fx)] <- fx[!is.na(fx)]
    # the coefficient of (v - m)/s of each level, then of ((v - m)/s)^2, ...: d_j s^j
    out$var_par <- as.vector(st) * rep(sc^seq_len(deg), each=G)
    out$map <- ifelse(is.na(as.vector(fx)), seq_len(G*deg), NA)
    out$names <- paste(rep(paste("poly", seq_len(deg), sep=""), each=G), rep(sfx, deg), sep="")
    attr(lev, "scale") <- sc
    attr(lev, "center") <- ce
  }
  out$levels <- lev
  out
}

# the values v of varPoly in a matrix of the levels by the powers 1, ..., deg, starting from target:
# names polyj for all levels, or polyj:level for one
varpolyvals <- function(v, target, levels, deg){
  if (length(v) == 0) return(target)
  pn <- paste("poly", seq_len(deg), sep="")
  for (nm in names(v)){
    j <- match(sub(":.*$", "", nm), pn)
    if (is.na(j)) stop("unknown name ", nm, " in the values of varPoly; the names are ", paste(pn, collapse=", "),
                       ", or e.g. poly1:<level> with a grouping factor.", call.=FALSE)
    rows <- if (grepl(":", nm, fixed=TRUE)) match(sub("^[^:]*:", "", nm), levels) else seq_along(levels)
    if (anyNA(rows)) stop("unknown level in ", nm, " of varPoly; the levels are ", paste(levels, collapse=", "), ".", call.=FALSE)
    target[rows, j] <- v[[nm]]
  }
  target
}

# the factors s^j of the coefficients of varPoly (its parameters are d_j s^j, for the polynomial in
# (v - m)/s), in the order of its parameters; levels with the attribute "scale" s
varpolyfactor <- function(levels, deg){
  sc <- attr(levels, "scale")
  if (is.null(sc)) sc <- 1
  rep(sc^seq_len(deg), each=length(levels))
}

# parameters of a variance function on the natural scale from the unconstrained values v, named as
# varpars, with the levels of each component
varnatural <- function(variance, v, levels){
  comps <- varcomps(variance)
  pre <- varprefix(variance)
  out <- numeric(0)
  off <- 0
  for (k in seq_along(comps)){
    vt <- vartypeone(comps[[k]])
    deg <- vardegree(comps[[k]])
    np <- varnpar(vt, length(levels[[k]]), deg)
    nat <- varnaturalone(vt, v[off + seq_len(np)], levels[[k]], deg)
    off <- off + np
    out <- c(out, setNames(nat, paste(pre[k], names(nat), sep="")))
  }
  out
}

# parameters of a single variance function of type vt on the natural scale (varIdent: the ratios of
# all strata, 1 for the reference)
varnaturalone <- function(vt, v, levels, deg){
  if (vt == 1) return(setNames(c(1, exp(v)), levels))
  G <- length(levels)
  sfx <- if (identical(c(levels), "1")) "" else paste(":", levels, sep="")
  if (vt == 5) v <- v / varpolyfactor(levels, deg)
  switch(vt, NULL,
         setNames(v, paste("power", sfx, sep="")),
         setNames(c(exp(v[seq_len(G)]), v[G + seq_len(G)]), c(paste("const", sfx, sep=""), paste("power", sfx, sep=""))),
         setNames(v, paste("expon", sfx, sep="")),
         setNames(v, paste(rep(paste("poly", seq_len(deg), sep=""), each=G), rep(sfx, deg), sep="")))
}

# the parameters of a variance function, one row each, as in varpars: the name, the position in var_par
# (NA for the reference stratum of varIdent), whether the natural value is exp of the unconstrained one,
# the name in coef (var_par:...), the name of the profile interval and what it is, for messages
varparinfo <- function(variance, levels){
  comps <- varcomps(variance)
  pre <- varprefix(variance)
  out <- NULL
  off <- 0L
  for (k in seq_along(comps)){
    vt <- vartypeone(comps[[k]])
    deg <- vardegree(comps[[k]])
    G <- length(levels[[k]])
    np <- varnpar(vt, G, deg)
    nm <- names(varnaturalone(vt, rep(0, np), levels[[k]], deg))
    idx <- if (vt == 1) c(NA, off + seq_len(np)) else off + seq_len(np)
    lg <- if (vt == 1) rep(TRUE, G) else if (vt == 3) rep(c(TRUE, FALSE), each=G) else rep(FALSE, np)
    what <- if (vt == 1) rep("ratio of varIdent", G) else if (vt == 3) rep(c("constant of varConstPower", "power of varConstPower"), each=G) else
      rep(c(NA, "power of varPower", NA, "exponent of varExp", "coefficient of varPoly")[vt], np)
    prof <- paste(pre[k], if (vt == 1) "ratio:", nm, sep="")
    fac <- if (vt == 5) varpolyfactor(levels[[k]], deg) else rep(1, length(nm))
    out <- rbind(out, data.frame(name=paste(pre[k], nm, sep=""), index=idx, log=lg, coefname=paste("var_par:", pre[k], nm, sep=""),
                                 profname=prof, what=what, factor=fac, stringsAsFactors=FALSE))
    off <- off + np
  }
  out
}

# the parameters of the variance function of a fit (see varparinfo), NULL without
fitvarinfo <- function(object){
  if (is.null(object$variance)) NULL else varparinfo(object$variance, object$varlevels)
}

# the standard deviation multiplier g of each row of data for the parameters v (unconstrained) of a
# variance function, the product of its components: with the fitted curve mu and its growth rate rate
# (a function, called only for rate(.)), as in the template
varmult <- function(variance, v, levels, data, mu, rate){
  n <- length(mu)
  g <- rep(1, n)
  comps <- varcomps(variance)
  off <- 0
  rv <- NULL
  ratec <- function(){
    if (is.null(rv)) rv <<- rate()
    rv
  }
  for (k in seq_along(comps)){
    vt <- vartypeone(comps[[k]])
    deg <- vardegree(comps[[k]])
    np <- varnpar(vt, length(levels[[k]]), deg)
    g <- g * varmultone(comps[[k]], vt, v[off + seq_len(np)], levels[[k]], deg, data, mu, ratec)
    off <- off + np
  }
  g
}

# the multiplier g of a single variance function of type vt (see varmult)
varmultone <- function(variance, vt, v, levels, deg, data, mu, rate){
  n <- length(mu)
  vf <- varform(attr(variance, "formula"))
  if (vt == 1){
    st <- varstrata(vf$group, data)
    k <- match(st$label, levels)
    if (anyNA(k)) stop("the variance function varIdent has no variance for the strata ", paste(unique(st$label[is.na(k)]), collapse=", "),
                       " of newdata; only for those of the fit: ", paste(levels, collapse=", "), ".", call.=FALSE)
    return(c(1, exp(v))[k])
  }
  ct <- varcovtype(vf$cov)
  cv <- if (ct == 1) mu else if (ct == 2) rate() else {
    x <- eval(vf$cov, data, environment(attr(variance, "formula")))
    if (!is.numeric(x) || length(x) != n) stop("newdata needs the covariate of the variance function.")
    x
  }
  k <- rep(1L, n)
  if (!is.null(vf$group)){
    k <- match(varstrata(vf$group, data)$label, levels)
    if (anyNA(k)) stop("the variance function has no parameters for the levels ", paste(unique(varstrata(vf$group, data)$label[is.na(k)]), collapse=", "),
                       " of newdata.", call.=FALSE)
  }
  G <- length(levels)
  if (vt == 2) return(abs(cv)^v[k])
  if (vt == 3) return(exp(v[k]) + abs(cv)^v[G + k])
  if (vt == 4) return(exp(v[k]*cv))
  p <- rep(1, n)
  v <- v / varpolyfactor(levels, deg)
  ce <- attr(levels, "center")
  if (!is.null(ce)) cv <- cv - ce
  for (j in seq_len(deg)) p <- p + v[(j - 1)*G + k]*cv^j
  bad <- !is.na(p) & p <= 0
  if (any(bad)) stop("the polynomial of varPoly is 0 or below at ", sum(bad), " of the ", n, " observations (covariate ",
                     paste(signif(utils::head(cv[bad], 3), 4), collapse=", "), if (sum(bad) > 3) ", ...", "), where the standard deviation is not defined.", call.=FALSE)
  p
}

# step of the central difference of the growth rate at time x, as in the template
ratestep <- function(x) 1e-5*(1 + abs(x))

# the multiplier g of the rows of newdata (default: the data of the fit) for the fitted curves: with
# random=TRUE the curves of the fitted groups, pd and eta from preddata (computed if NULL)
fitvarmult <- function(object, newdata=object$data, random=TRUE, pd=NULL, eta=NULL, mu=NULL){
  if (is.null(object$variance)) return(rep(1, if (!is.null(mu)) length(mu) else nrow(newdata)))
  if (is.null(pd)) pd <- preddata(object, newdata, random)
  if (is.null(eta)) eta <- pd$eta
  if (is.null(mu)) mu <- predmean(pd, eta)
  rate <- function(){
    h <- ratestep(pd$x)
    p1 <- p2 <- pd
    p1$x <- pd$x + h
    p2$x <- pd$x - h
    (predmean(p1, eta) - predmean(p2, eta))/(2*h)
  }
  v <- object$obj$env$parList(object$opt$par)$var_par
  varmult(object$variance, v, object$varlevels, newdata, mu, rate)
}

# description of the variance function of a fit, for print and summary
vartext <- function(object){
  vs <- object$variance
  if (vartypeindex(vs) == 0) return(NULL)
  comps <- varcomps(vs)
  txt <- vapply(seq_along(comps), function(k) vartextone(comps[[k]], object$varlevels[[k]]), character(1))
  if (length(txt) == 1) return(txt)
  paste("the product of ", paste(LETTERS[seq_along(txt)], ". ", txt, sep="", collapse="; "), sep="")
}

# the centre of the covariate of varPoly in its description: " - m", or nothing for m = 0
varcentertext <- function(m){
  if (is.null(m) || m == 0) return("")
  paste(if (m > 0) " - " else " + ", format(signif(abs(m), 4)), sep="")
}

# description of a single variance function, with its levels (the centre of varPoly)
vartextone <- function(vs, levels=NULL){
  vt <- vartypeone(vs)
  vf <- varform(attr(vs, "formula"))
  cv <- if (!is.null(vf$cov)) deparse1(vf$cov)
  grp <- if (!is.null(vf$group)) paste(", by", deparse1(vf$group))
  switch(vt,
         paste("different standard deviations per stratum of", deparse1(vf$group)),
         paste("power of ", cv, grp, sep=""),
         paste("constant plus power of ", cv, grp, sep=""),
         paste("exponential of ", cv, grp, sep=""),
         paste("polynomial of degree ", vardegree(vs), " in ", cv, varcentertext(attr(levels, "center")), grp, sep=""))
}

printvar <- function(object, digits){
  txt <- vartext(object)
  if (is.null(txt)) return(invisible())
  cat("\nVariance function: ", txt, "\n", sep="")
  vp <- if (!is.null(object$varpars)) object$varpars else varpars(object)
  cat(paste(names(vp), dec2text(vp), sep=": ", collapse=", "), "\n", sep="")
}
