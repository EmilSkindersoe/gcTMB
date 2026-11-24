#' Autocorrelation of the residuals within groups
#'
#' Correlation structures for the residuals of the observations of a curve over time, for the
#' \code{correlation} argument of \code{\link{gcTMB}}, with the same arguments and the same
#' parameterisation as \code{corAR1} and \code{corCAR1} of the package nlme (whose objects
#' \code{\link{gcTMB}} also accepts).
#'
#' The residuals of different groups are independent. Within a group, the residuals e are stationary,
#' with the residual variance sigma^2 of the model (sigma^2/w for weights w) and the correlation
#' \deqn{cor(e_s, e_t) = rho^{|t - s|},}{cor(e_s, e_t) = rho^|t - s|,}
#' with rho the correlation of residuals one unit of time apart (Zuur et al., 2009, Section 6.1;
#' Pinheiro and Bates, 2000, Section 5.3):
#' \itemize{
#'   \item \code{corAR1}: the autoregressive process of order 1 at integer times (positions), with
#'     -1 < rho < 1; gaps in the times are allowed, e.g. measurements every second day with
#'     \code{form = ~ day | ID} and a day without measurement, with the correlation rho^2 across a
#'     gap of two days. Without a time covariate (\code{form = ~ 1 | ID}), the time is the position
#'     of the observation within its group, in the order of the data;
#'   \item \code{corCAR1}: the continuous-time autoregressive process of order 1, for times that need
#'     not be integers and unequally spaced observations, with 0 < rho < 1.
#' }
#' The times must be unique within a group. As the process is Markov, the likelihood is the product of
#' the density of the first residual of each group and the conditional densities of each residual given
#' the previous one, normal with mean rho^dt e_prev and variance sigma^2 (1 - rho^(2 dt)) for the time
#' dt between them, which is exact and linear in the number of observations.
#'
#' @param value Numeric. The starting value of rho, or its value with \code{fixed=TRUE}: between -1
#'   and 1 (default 0) for \code{corAR1}, between 0 and 1 (default 0.2) for \code{corCAR1}
#' @param form One-sided formula. The time covariate and the grouping factor, \code{~ time | group}.
#'   Without a grouping factor, the residuals of all observations form a single series when the model has
#'   no random effects, and are grouped by the innermost grouping factor of the random effects (the one
#'   with most levels) otherwise, as in \code{nlme::lme}. Without a time covariate (\code{~ 1} or
#'   \code{~ 1 | group}), the times are the positions of the observations within their group.
#' @param fixed Logical. If TRUE, rho is fixed at \code{value} and not estimated
#'
#' @return An object of class \code{c("corAR1", "corStruct")} or \code{c("corCAR1", "corStruct")}:
#'   \code{value} on the unconstrained scale of nlme, log((1 + rho)/(1 - rho)) for AR(1) and
#'   log(rho/(1 - rho)) for CAR(1), with the attributes \code{formula} and \code{fixed}.
#' @references Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS.
#'   Springer, New York.
#'
#'   Zuur, A.F., Ieno, E.N., Walker, N.J., Saveliev, A.A. and Smith, G.M. (2009). Mixed Effects Models
#'   and Extensions in Ecology with R. Springer, New York.
#' @details The estimated parameter is in \code{coef(fit)} as \code{cor_par:AR1} or \code{cor_par:CAR1}, on
#' the unconstrained scale of nlme (\code{coef} of an nlme \code{corStruct}): the logit link
#' log(rho/(1 - rho)) for CAR(1), and log((1 + rho)/(1 - rho)), the logit of (1 + rho)/2, for AR(1),
#' whose rho lies between -1 and 1. Wald intervals are computed on this scale and transformed back.
#' \code{\link{corpars}} gives
#' rho, and \code{confint(fit, parm="cor_par:CAR1")} an interval for rho. If nlme is attached after
#' gcTMB, its constructors mask these; gcTMB accepts their objects, and gcTMB's correlation structures can
#' be passed to \code{nlme::gls}, \code{nlme::gnls}, \code{nlme::lme} and \code{nlme::nlme}.
#' @seealso \code{\link{gcTMB}}, \code{\link{corpars}} for the estimate of rho, and
#'   \code{residuals(fit, type="normalized")} for residuals that are independent under the model.
#' @examples
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fit <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25,
#'              correlation=corCAR1(form=~ DAS | plant_ID))
#' corpars(fit)
#' @name corStruct
#' @aliases corAR1 corCAR1
NULL

#' @rdname corStruct
#' @export
corAR1 <- function(value=0, form=~ 1, fixed=FALSE){
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || abs(value) >= 1) stop("parameter in AR(1) structure must be between -1 and 1")
  value <- log((1 + value)/(1 - value))
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  class(value) <- c("corAR1", "corStruct")
  value
}

#' @rdname corStruct
#' @export
corCAR1 <- function(value=0.2, form=~ 1, fixed=FALSE){
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value <= 0 || value >= 1) stop("parameter in CAR(1) structure must be between 0 and 1")
  value <- log(value/(1 - value))
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  class(value) <- c("corCAR1", "corStruct")
  value
}

# type of a correlation structure: 0 none, 1 AR(1), 2 CAR(1)
cortypeindex <- function(correlation){
  if (is.null(correlation)) return(0L)
  if (inherits(correlation, "corAR1")) return(1L)
  if (inherits(correlation, "corCAR1")) return(2L)
  stop("correlation must be corAR1() or corCAR1().")
}

# rho from its unconstrained value v, for correlation type ct (1 AR(1), 2 CAR(1))
corrho <- function(ct, v){
  if (ct == 1) (exp(v) - 1)/(exp(v) + 1) else stats::plogis(v)
}

# parts of the form ~ time | group of a correlation structure: the time expression (NULL for ~ 1)
# and the grouping expression (NULL without a grouping factor)
corform <- function(correlation){
  form <- attr(correlation, "formula")
  if (is.null(form)) form <- ~ 1
  if (!inherits(form, "formula") || length(form) != 2) stop("the form of the correlation structure must be a one-sided formula, e.g. ~ time | group.")
  rhs <- form[[2]]
  grp <- NULL
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))){
    grp <- rhs[[3]]
    rhs <- rhs[[2]]
  }
  tm <- if (identical(rhs, 1) || identical(rhs, 1L)) NULL else rhs
  list(time=tm, group=grp)
}

# variables of a correlation structure, for the missing values of the data
corvars <- function(correlation){
  if (is.null(correlation)) return(character(0))
  cf <- corform(correlation)
  c(if (!is.null(cf$time)) all.vars(cf$time), if (!is.null(cf$group)) all.vars(cf$group))
}

# group of each row of data: the levels of the grouping expression of the form (with / or : for
# nested factors, the innermost groups); without one, the innermost grouping factor of the random effect
# terms (the term with most levels), or a single group
corgroups <- function(cf, data, reterms, termgroup){
  n <- nrow(data)
  if (!is.null(cf$group)){
    vars <- all.vars(cf$group)
    miss <- setdiff(vars, names(data))
    if (length(miss) > 0) stop("the data need the grouping variables of the correlation structure: ", paste(miss, collapse=", "), ".")
    return(as.character(interaction(lapply(vars, function(v) data[[v]]), drop=TRUE, lex.order=TRUE)))
  }
  if (length(reterms) > 0){
    t <- which.max(vapply(reterms, function(tm) length(tm$levels), integer(1)))
    return(termgroup(t, data))
  }
  rep("1", n)
}

# the chains of the correlated residuals: for each observation of data its group, time (or position),
# the previous observation of its group in time (corprev, zero-based, -1 for the first), the time or
# lag since then, and an order with each observation after its previous one. Observations with weight
# zero are left out of the chains.
corsetup <- function(correlation, data, w=rep(1, nrow(data)), reterms=list(), termgroup=NULL){
  n <- nrow(data)
  ct <- cortypeindex(correlation)
  empty <- list(cortype=0L, corprev=rep(-1L, n), cordt=rep(0, n), corlag=rep(0L, n), corord=seq_len(n) - 1L)
  if (ct == 0) return(empty)
  cf <- corform(correlation)
  g <- corgroups(cf, data, reterms, termgroup)
  use <- w > 0
  if (is.null(cf$time)){
    # positions within the group, in the order of the data
    tt <- numeric(n)
    tt[use] <- stats::ave(seq_len(sum(use)), g[use], FUN=seq_along)
  } else {
    tt <- eval(cf$time, data, environment(attr(correlation, "formula")))
    if (!is.numeric(tt) || length(tt) != n) stop("the time covariate of the correlation structure must be numeric.")
  }
  cls <- class(correlation)[1]
  if (any(duplicated(data.frame(g, tt)[use, ]))) stop("the time covariate must have unique values within groups for \"", cls, "\" objects.")
  if (ct == 1 && any(abs(tt[use] - round(tt[use])) > 1e-8)) stop("the time covariate of \"corAR1\" must be integer valued; use corCAR1() for continuous times.")
  prev <- rep(-1L, n)
  dt <- rep(0, n)
  ord <- order(g, tt)
  ord <- ord[use[ord]]
  for (grp in split(ord, g[ord])){
    if (length(grp) > 1){
      prev[grp[-1]] <- grp[-length(grp)] - 1L
      dt[grp[-1]] <- diff(tt[grp])
    }
  }
  list(cortype=ct, corprev=as.integer(prev), cordt=as.numeric(dt), corlag=as.integer(round(dt)),
       corord=as.integer(c(ord, which(!use)) - 1L), group=g, time=tt)
}

# description of the correlation structure of a fit, for print and summary
cortext <- function(object){
  cs <- object$correlation
  if (is.null(cs)) return(NULL)
  cf <- corform(cs)
  ct <- cortypeindex(cs)
  grp <- if (!is.null(cf$group)) deparse1(cf$group) else if (length(object$reterms) > 0) object$corgroup else NULL
  paste(c("AR(1)", "continuous-time AR(1)")[ct],
        if (!is.null(cf$time)) paste(" in", deparse1(cf$time)) else " in the order of the observations",
        if (!is.null(grp)) paste(" within", grp), sep="")
}

printcor <- function(object, digits){
  txt <- cortext(object)
  if (is.null(txt)) return(invisible())
  cat("\nCorrelation structure: ", txt, "\n", sep="")
  cp <- if (!is.null(object$corrho)) object$corrho else corpars(object)
  cat("rho: ", format(cp[["rho"]], digits=digits),
      if (isTRUE(attr(object$correlation, "fixed"))) " (fixed)", "\n", sep="")
}

# residuals that are independent under the model: the standardised residuals r (correlated within the
# chains of the fit) minus their conditional mean given the previous residual, divided by the
# conditional standard deviation
normalizeresid <- function(object, r){
  dat <- object$obj$env$data
  ct <- dat$cortype
  if (is.null(ct) || ct == 0 || !(object$family %in% c("gaussian", "lognormal"))) return(r)
  rho <- corpars(object)[["rho"]]
  prev <- dat$corprev
  has <- prev >= 0
  cr <- rep(0, length(r))
  cr[has] <- if (ct == 1) rho^dat$corlag[has] else rho^dat$cordt[has]
  rp <- rep(0, length(r))
  rp[has] <- r[prev[has] + 1]
  out <- (r - cr*rp)/sqrt(1 - cr^2)
  out[dat$w <= 0] <- 0
  out
}

# standard normal residuals with the correlation of the fit for the rows of newdata (simulate with
# newdata): chains by group and time as in the fit, with weight 1
corsimresid <- function(object, newdata){
  n <- nrow(newdata)
  if (is.null(object$correlation)) return(stats::rnorm(n))
  cs <- corsetup(object$correlation, newdata, reterms=object$reterms,
                 termgroup=function(t, data) obsgroups(object, t, data))
  rho <- corpars(object)[["rho"]]
  r <- numeric(n)
  for (i in cs$corord + 1L){
    p <- cs$corprev[i]
    if (p < 0){
      r[i] <- stats::rnorm(1)
    } else {
      cr <- if (cs$cortype == 1) rho^cs$corlag[i] else rho^cs$cordt[i]
      r[i] <- cr*r[p + 1] + sqrt(1 - cr^2)*stats::rnorm(1)
    }
  }
  r
}
