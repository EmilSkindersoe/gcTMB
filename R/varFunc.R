#' Variance functions for heteroscedastic residuals
#'
#' Variance functions for the \code{variance} argument of \code{\link{gcTMB}}, with the arguments and
#' the parameterisation of \code{varIdent}, \code{varPower} and \code{varConstPower} of the package
#' nlme (whose objects \code{\link{gcTMB}} also accepts), for the gaussian and lognormal (on the log
#' scale) families. The residual variance of an observation is sigma^2 g^2, with the residual standard
#' deviation sigma of the model and the function g:
#' \itemize{
#'   \item \code{varIdent(form=~ 1 | s)}: g = delta_s, a ratio for each stratum s of the grouping
#'     variables, with delta = 1 for the reference stratum, so that each stratum has its own variance,
#'     e.g. a variance for each treatment (\code{~ 1 | Treatment}) or for each time point
#'     (\code{~ 1 | Time}); a numeric grouping variable is treated as a factor;
#'   \item \code{varPower(form=~ v)}: g = |v|^delta, variance sigma^2 |v|^(2 delta);
#'   \item \code{varConstPower(form=~ v)}: g = c + |v|^delta, c > 0, variance sigma^2 (c + |v|^delta)^2.
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
#'   default); for \code{varPower} the power (0 by default), named by level with a grouping factor
#' @param const,power Numeric. For \code{varConstPower}, the starting values of the constant c (0.1 by
#'   default, as nlme) and the power delta (0 by default)
#' @param form One-sided formula: \code{~ 1 | s} for \code{varIdent}, with s one or more grouping
#'   variables (\code{~ 1 | A * B}); \code{~ v} or \code{~ v | g} for \code{varPower} and
#'   \code{varConstPower}, with v a numeric variable, \code{fitted(.)} or \code{rate(.)}
#' @param fixed Values that are fixed and not estimated: for \code{varIdent} ratios named by stratum,
#'   for \code{varPower} the power (named by level with a grouping factor), for \code{varConstPower}
#'   a list with \code{const} and/or \code{power}, e.g. \code{fixed=list(power=1)}
#'
#' @details The estimated parameters are in \code{coef(fit)} as \code{var_par:<stratum>} (the log
#' ratios of \code{varIdent}), \code{var_par:power} and \code{var_par:const} (the log of c), with
#' \code{:<level>} for a grouping factor, on the unconstrained scale of nlme; \code{\link{varpars}}
#' gives them on the natural scale.
#' @return An object of class \code{c("varIdent", "varFunc")}, \code{c("varPower", "varFunc")} or
#'   \code{c("varConstPower", "varFunc")}, with the structure of the nlme objects.
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
#' @name varFunc
#' @aliases varIdent varPower varConstPower
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

# type of a variance function: 0 none, 1 varIdent, 2 varPower, 3 varConstPower
vartypeindex <- function(variance){
  if (is.null(variance)) return(0L)
  if (inherits(variance, "varIdent")) return(1L)
  if (inherits(variance, "varPower")) return(2L)
  if (inherits(variance, "varConstPower")) return(3L)
  stop("variance must be varIdent(), varPower() or varConstPower().")
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
  if (is.null(variance)) return(character(0))
  vf <- varform(attr(variance, "formula"))
  v <- c(if (varcovtype(vf$cov) == 0) all.vars(vf$cov), all.vars(vf$group))
  setdiff(v, ".")
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

# the data of a variance function for the template: vartype, the covariate type varcovt, the covariate
# varcov, the stratum (varIdent) or level (varPower, varConstPower) varstr of each row (zero-based,
# the reference stratum 0), the number of levels varng; and the parameters var_par (on the
# unconstrained scale), their map (NA for fixed) and names
varsetup <- function(variance, data){
  n <- nrow(data)
  vt <- vartypeindex(variance)
  out <- list(vartype=vt, varcovt=0L, varcov=rep(0, n), varstr=rep(0L, n), varng=1L,
              var_par=numeric(0), map=integer(0), names=character(0), levels=NULL)
  if (vt == 0) return(out)
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
    out$names <- paste("var_par:", others, sep="")
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
  if (vt == 2){
    lv <- varlevelvals(setNames(as.numeric(variance), names(variance)), attr(variance, "fixed"), lev, 0, "power")
    if (is.null(vf$group)) names(lv$start) <- names(lv$fixed) <- NULL
    out$var_par <- unname(lv$start)
    out$map <- ifelse(is.na(lv$fixed), seq_len(G), NA)
    out$names <- paste("var_par:power", sfx, sep="")
  } else {
    vc <- varlevelvals(variance$const, attr(variance$const, "fixed"), lev, log(0.1), "const")
    vp <- varlevelvals(variance$power, attr(variance$power, "fixed"), lev, 0, "power")
    out$var_par <- unname(c(vc$start, vp$start))
    out$map <- ifelse(is.na(c(vc$fixed, vp$fixed)), seq_len(2*G), NA)
    out$names <- c(paste("var_par:const", sfx, sep=""), paste("var_par:power", sfx, sep=""))
  }
  out$levels <- lev
  out
}

# parameters on the natural scale from the unconstrained values v
varnatural <- function(vt, v, levels, G){
  if (vt == 1) return(setNames(c(1, exp(v)), levels))
  sfx <- if (identical(levels, "1")) "" else paste(":", levels, sep="")
  if (vt == 2) return(setNames(v, paste("power", sfx, sep="")))
  setNames(c(exp(v[seq_len(G)]), v[G + seq_len(G)]), c(paste("const", sfx, sep=""), paste("power", sfx, sep="")))
}

# the standard deviation multiplier g of each row of data for the parameters v (unconstrained) of a
# variance function: with the fitted curve mu and its growth rate rate (a function, called only for
# rate(.)), as in the template
varmult <- function(variance, v, levels, data, mu, rate){
  vt <- vartypeindex(variance)
  n <- length(mu)
  if (vt == 0) return(rep(1, n))
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
  exp(v[k]) + abs(cv)^v[G + k]
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
  vt <- vartypeindex(vs)
  if (vt == 0) return(NULL)
  vf <- varform(attr(vs, "formula"))
  cv <- if (!is.null(vf$cov)) deparse1(vf$cov)
  grp <- if (!is.null(vf$group)) deparse1(vf$group)
  switch(vt,
         paste("different standard deviations per stratum of", grp),
         paste("power of ", cv, if (!is.null(grp)) paste(", by", grp), sep=""),
         paste("constant plus power of ", cv, if (!is.null(grp)) paste(", by", grp), sep=""))
}

printvar <- function(object, digits){
  txt <- vartext(object)
  if (is.null(txt)) return(invisible())
  cat("\nVariance function: ", txt, "\n", sep="")
  vp <- if (!is.null(object$varpars)) object$varpars else varpars(object)
  cat(paste(names(vp), format(vp, digits=digits), sep=": ", collapse=", "), "\n", sep="")
}
