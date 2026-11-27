#' Autocorrelation of the residuals within groups
#'
#' Correlation structures for the residuals of the observations of a curve over time (or in space), for
#' the \code{correlation} argument of \code{\link{gcTMB}}, with the same arguments and the same
#' parameterisation as \code{corAR1}, \code{corCAR1}, \code{corExp} and \code{corGaus} of the package
#' nlme (whose objects \code{\link{gcTMB}} also accepts).
#'
#' The residuals of different groups are independent. Within a group, the residuals e are stationary,
#' with the residual variance sigma^2 of the model (sigma^2/w for weights w). For \code{corAR1} and
#' \code{corCAR1}, the correlation is
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
#' The spatial structures \code{corExp} and \code{corGaus} give the correlation as a function of the
#' distance d between the positions of two observations of a group, with a range r > 0 and, with
#' \code{nugget=TRUE}, a nugget 0 < n < 1 (n = 0 otherwise):
#' \itemize{
#'   \item \code{corExp}: \eqn{cor = (1 - n) \exp(-d/r)}{cor = (1 - n) exp(-d/r)};
#'   \item \code{corGaus}: \eqn{cor = (1 - n) \exp(-(d/r)^2)}{cor = (1 - n) exp(-(d/r)^2)}.
#' }
#' The position is the time covariate (\code{~ time | group}), with d the time between the observations,
#' or several covariates (\code{~ x + y | group}), with d the distance of \code{metric}. Without a nugget,
#' \code{corExp} in time is \code{corCAR1} with rho = exp(-1/r). The nugget is the part of the residual
#' variance that is independent between observations, such as measurement error, besides the correlated
#' deviations of the process (Diggle, 1988); with it, observations at the same position are allowed. The
#' residuals of a group are multivariate normal with this correlation matrix, so the time of a fit grows
#' with the cube of the number of observations per group. As in nlme, the default starting value of the
#' range is 90\% of the smallest positive distance between observations of a group, and of the nugget 0.1.
#'
#' @param value Numeric. The starting value of the parameters, or their values with \code{fixed=TRUE}:
#'   rho, between -1 and 1 (default 0) for \code{corAR1} and between 0 and 1 for \code{corCAR1}, by
#'   default the rho with a correlation of 0.2 at the median time d between successive observations of
#'   a group, \eqn{0.2^{1/d}}{0.2^(1/d)}, so that the start does not depend on the unit of time (0.2 if
#'   \code{fixed=TRUE}, as in nlme); for \code{corExp} and \code{corGaus} the range, or with a nugget the range and the
#'   nugget (default: see Details)
#' @param form One-sided formula. The time covariate (or the position variables, \code{~ x + y | group}, for
#'   \code{corExp} and \code{corGaus}) and the grouping factor, \code{~ time | group}.
#'   Without a grouping factor, the residuals of all observations form a single series when the model has
#'   no random effects, and are grouped by the innermost grouping factor of the random effects (the one
#'   with most levels) otherwise, as in \code{nlme::lme}. Without a time covariate (\code{~ 1} or
#'   \code{~ 1 | group}), the times are the positions of the observations within their group.
#' @param nugget Logical. If TRUE, the spatial structure has a nugget
#' @param metric Character. The distance between positions of several variables: "euclidean" (default),
#'   "maximum" or "manhattan", as in \code{\link[stats]{dist}}; all give |t - s| for a single time covariate
#' @param fixed Logical. If TRUE, the parameters are fixed at \code{value} and not estimated
#'
#' @return An object of class \code{c("corAR1", "corStruct")} or \code{c("corCAR1", "corStruct")}:
#'   \code{value} on the unconstrained scale of nlme, log((1 + rho)/(1 - rho)) for AR(1) and
#'   log(rho/(1 - rho)) for CAR(1), with the attributes \code{formula} and \code{fixed}; or of class
#'   \code{c("corExp", "corSpatial", "corStruct")} or \code{c("corGaus", "corSpatial", "corStruct")}:
#'   \code{value} on the natural scale, as in nlme, with the attributes \code{formula}, \code{nugget},
#'   \code{metric} and \code{fixed}.
#' @references Diggle, P.J. (1988). An approach to the analysis of repeated measurements. Biometrics
#'   44, 959-971.
#'
#'   Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS.
#'   Springer, New York.
#'
#'   Zuur, A.F., Ieno, E.N., Walker, N.J., Saveliev, A.A. and Smith, G.M. (2009). Mixed Effects Models
#'   and Extensions in Ecology with R. Springer, New York.
#' @details The estimated parameters are in \code{coef(fit)} as \code{cor_par:AR1} or \code{cor_par:CAR1},
#' or \code{cor_par:range} and \code{cor_par:nugget}, on
#' the unconstrained scale of nlme (\code{coef} of an nlme \code{corStruct}): the logit link
#' log(rho/(1 - rho)) for CAR(1), and log((1 + rho)/(1 - rho)), the logit of (1 + rho)/2, for AR(1),
#' whose rho lies between -1 and 1; log(r) for the range and log(n/(1 - n)) for the nugget. Wald
#' intervals are computed on this scale and transformed back. \code{\link{corpars}} gives
#' rho (or the range and nugget), and \code{confint(fit, parm="cor_par:CAR1")} an interval for rho. If nlme is attached after
#' gcTMB, its constructors mask these; gcTMB accepts their objects, and gcTMB's correlation structures can
#' be passed to \code{nlme::gls}, \code{nlme::gnls}, \code{nlme::lme} and \code{nlme::nlme}.
#' @seealso \code{\link{gcTMB}}, \code{\link{corpars}} for the estimate of rho, and
#'   \code{residuals(fit, type="normalized")} for residuals that are independent under the model.
#' @examples
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fit <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25,
#'              correlation=corCAR1(form=~ DAS | Plant_ID))
#' corpars(fit)
#' # spatial correlation between neighbouring trees of a plot, with a nugget (independent
#' # measurement error), on five of the plots
#' pines <- droplevels(subset(PineHeightDiameter, Plot %in% 1:5))
#' fitn <- gcTMB(Height ~ DBH, rform=b1 ~ 1 | Plot, data=pines, model="chapmanrichards",
#'               correlation=corExp(form=~ X + Y | Plot, nugget=TRUE))
#' corpars(fitn)
#' @name corStruct
#' @aliases corAR1 corCAR1 corExp corGaus
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
  # the default start is replaced in the fit by the rho that gives a correlation of 0.2 at the median
  # spacing of the times (carstart), unless fixed
  bydefault <- missing(value)
  value <- log(value/(1 - value))
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  if (bydefault) attr(value, "spacingstart") <- TRUE
  class(value) <- c("corCAR1", "corStruct")
  value
}

#' @rdname corStruct
#' @export
corExp <- function(value=numeric(0), form=~ 1, nugget=FALSE, metric=c("euclidean", "maximum", "manhattan"), fixed=FALSE){
  corspatial(value, form, nugget, match.arg(metric), fixed, "corExp")
}

#' @rdname corStruct
#' @export
corGaus <- function(value=numeric(0), form=~ 1, nugget=FALSE, metric=c("euclidean", "maximum", "manhattan"), fixed=FALSE){
  corspatial(value, form, nugget, match.arg(metric), fixed, "corGaus")
}

# spatial correlation structure as in nlme: value on the natural scale (range, and nugget), transformed
# when the model is set up
corspatial <- function(value, form, nugget, metric, fixed, cls){
  spatialcheck(value, nugget)
  attr(value, "formula") <- form
  attr(value, "nugget") <- nugget
  attr(value, "metric") <- metric
  attr(value, "fixed") <- fixed
  class(value) <- c(cls, "corSpatial", "corStruct")
  value
}

spatialcheck <- function(value, nugget){
  if (!is.numeric(value) || any(!is.finite(value))) stop("the value of a \"corSpatial\" structure must be numeric.")
  if (length(value) == 0) return(invisible())
  if (value[1] <= 0) stop("'range' must be > 0 in \"corSpatial\" initial value")
  if (isTRUE(nugget)){
    if (length(value) > 2) stop("initial value for \"corSpatial\" parameters of wrong dimension")
    if (length(value) == 2 && (value[2] <= 0 || value[2] >= 1)) stop("initial value of nugget ratio must be between 0 and 1")
  } else if (length(value) != 1) stop("initial value for \"corSpatial\" parameters of wrong dimension")
  invisible()
}

# type of a correlation structure: 0 none, 1 AR(1), 2 CAR(1), 3 exponential, 4 Gaussian
cortypeindex <- function(correlation){
  if (is.null(correlation)) return(0L)
  if (inherits(correlation, "corAR1")) return(1L)
  if (inherits(correlation, "corCAR1")) return(2L)
  if (inherits(correlation, "corExp")) return(3L)
  if (inherits(correlation, "corGaus")) return(4L)
  stop("correlation must be corAR1(), corCAR1(), corExp() or corGaus().")
}

corhasnugget <- function(correlation) cortypeindex(correlation) >= 3 && isTRUE(attr(correlation, "nugget"))

# parameters on the natural scale from the unconstrained values v, for correlation type ct: rho for
# AR(1) and CAR(1), the range (and nugget) for the spatial structures
corrho <- function(ct, v){
  if (ct == 1) return(c(rho=(exp(v[1]) - 1)/(exp(v[1]) + 1)))
  if (ct == 2) return(c(rho=stats::plogis(v[1])))
  c(range=exp(v[1]), if (length(v) > 1) c(nugget=stats::plogis(v[2])))
}

# starting (or fixed) value of a spatial structure on the unconstrained scale, log(range) and
# logit(nugget): the value given, or as nlme, 90% of the smallest distance and a nugget of 0.1
spatialstart <- function(correlation, mind){
  nug <- isTRUE(attr(correlation, "nugget"))
  val <- as.numeric(correlation)
  if (!is.null(attr(correlation, "minD"))) return(val)   # an initialised nlme object: unconstrained already
  spatialcheck(val, nug)
  if (length(val) == 0) val <- 0.9*mind
  if (nug && length(val) == 1) val <- c(val, 0.1)
  c(log(val[1]), if (nug) stats::qlogis(val[2]))
}

# starting (or fixed) value of CAR(1) on the unconstrained scale, logit(rho): the value given, or, for
# the default of corCAR1() (not fixed), the rho with the correlation 0.2 at the median time between
# successive observations of a group, rho = 0.2^(1/d), so that the start does not depend on the unit of
# time (0.2 per unit of time is practically no correlation between observations many units apart)
carstart <- function(correlation, cor){
  val <- as.numeric(correlation)
  if (!isTRUE(attr(correlation, "spacingstart")) || isTRUE(attr(correlation, "fixed"))) return(val)
  dt <- cor$cordt[cor$corprev >= 0]
  dt <- dt[dt > 0]
  if (length(dt) == 0) return(val)
  stats::qlogis(min(0.2^(1/stats::median(dt)), 0.999))
}

# correlation matrix of a spatial structure (type ct) for the distances D and the parameters cp
# (range, nugget) on the natural scale
spatialcor <- function(ct, D, cp){
  nug <- if ("nugget" %in% names(cp)) cp[["nugget"]] else 0
  d <- D/cp[["range"]]
  R <- (1 - nug)*(if (ct == 3) exp(-d) else exp(-d^2))
  diag(R) <- 1
  R
}

# the expressions of the position variables of a spatial form, e.g. x and y of ~ x + y | g
spatialterms <- function(e){
  if (is.call(e) && identical(e[[1]], as.name("+")) && length(e) == 3) c(spatialterms(e[[2]]), spatialterms(e[[3]])) else list(e)
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
  empty <- list(cortype=0L, corprev=rep(-1L, n), cordt=rep(0, n), corlag=rep(0L, n), corord=seq_len(n) - 1L,
                corgsize=integer(0), cordist=numeric(0))
  if (ct == 0) return(empty)
  cf <- corform(correlation)
  g <- corgroups(cf, data, reterms, termgroup)
  use <- w > 0
  if (ct >= 3) return(spatialsetup(correlation, cf, data, g, use, empty))
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
       corord=as.integer(c(ord, which(!use)) - 1L), corgsize=integer(0), cordist=numeric(0), group=g, time=tt)
}

# the groups of a spatial structure: the observations of each group consecutive in corord (with the
# sizes corgsize), and the matrix of distances between the positions of each group, by column,
# concatenated (cordist); mind is the smallest positive distance, for the starting value of the range
spatialsetup <- function(correlation, cf, data, g, use, empty){
  n <- nrow(data)
  if (is.null(cf$time)){
    pos <- matrix(0, n, 1)
    pos[use, 1] <- stats::ave(seq_len(sum(use)), g[use], FUN=seq_along)
  } else {
    env <- environment(attr(correlation, "formula"))
    pos <- sapply(spatialterms(cf$time), function(e){
      v <- eval(e, data, env)
      if (!is.numeric(v) || length(v) != n) stop("the position variables of the correlation structure must be numeric.")
      v
    })
    pos <- matrix(pos, nrow=n)
  }
  metric <- attr(correlation, "metric")
  if (is.null(metric)) metric <- "euclidean"
  ord <- order(g, pos[, 1])
  ord <- ord[use[ord]]
  blocks <- split(ord, factor(g[ord], levels=unique(g[ord])))
  dl <- lapply(blocks, function(b) as.matrix(stats::dist(pos[b, , drop=FALSE], method=metric)))
  dd <- unlist(lapply(dl, function(D) D[lower.tri(D)]))
  if (length(dd) > 0 && any(dd == 0) && !isTRUE(attr(correlation, "nugget"))){
    stop("cannot have zero distances in \"corSpatial\" without a nugget; give the observations of a group distinct positions, or use nugget=TRUE.")
  }
  pd <- dd[dd > 0]
  mind <- if (length(pd) > 0) min(pd) else 1
  out <- empty
  out$cortype <- cortypeindex(correlation)
  out$corord <- as.integer(c(unlist(blocks, use.names=FALSE), which(!use)) - 1L)
  out$corgsize <- as.integer(lengths(blocks))
  out$cordist <- as.numeric(unlist(lapply(dl, as.vector), use.names=FALSE))
  out$cormind <- mind
  out$group <- g
  out$time <- pos
  out
}

# description of the correlation structure of a fit, for print and summary
cortext <- function(object){
  cs <- object$correlation
  if (is.null(cs)) return(NULL)
  cf <- corform(cs)
  ct <- cortypeindex(cs)
  grp <- if (!is.null(cf$group)) deparse1(cf$group) else if (length(object$reterms) > 0) object$corgroup else NULL
  metric <- attr(cs, "metric")
  paste(c("AR(1)", "continuous-time AR(1)", "exponential", "Gaussian")[ct],
        if (!is.null(cf$time)) paste(" in", deparse1(cf$time)) else " in the order of the observations",
        if (ct >= 3 && !is.null(metric) && metric != "euclidean") paste(" (", metric, " distance)", sep=""),
        if (!is.null(grp)) paste(" within", grp),
        if (corhasnugget(cs)) ", with nugget", sep="")
}

printcor <- function(object, digits){
  txt <- cortext(object)
  if (is.null(txt)) return(invisible())
  cat("\nCorrelation structure: ", txt, "\n", sep="")
  cp <- if (!is.null(object$corrho)) object$corrho else corpars(object)
  cat(paste(names(cp), dec2text(cp), sep=": ", collapse=", "),
      if (isTRUE(attr(object$correlation, "fixed"))) " (fixed)", "\n", sep="")
}

# applies fun(idx, R) to each group of a spatial structure: the (one-based) observations idx of the group
# and their correlation matrix R, from the TMB data dat (or corsetup) and the parameters cp
spatialblocks <- function(dat, cp, fun){
  off <- 0
  doff <- 0
  for (m in dat$corgsize){
    idx <- dat$corord[off + seq_len(m)] + 1L
    R <- spatialcor(dat$cortype, matrix(dat$cordist[doff + seq_len(m*m)], m, m), cp)
    fun(idx, R)
    off <- off + m
    doff <- doff + m*m
  }
}

# residuals that are independent under the model: the standardised residuals r (correlated within the
# chains of the fit) minus their conditional mean given the previous residual, divided by the
# conditional standard deviation
normalizeresid <- function(object, r){
  dat <- object$obj$env$data
  ct <- dat$cortype
  if (is.null(ct) || ct == 0 || !(object$family %in% c("gaussian", "lognormal"))) return(r)
  if (ct >= 3){
    # L^-1 r per group, with R = L L' the Cholesky factorisation of the correlation matrix
    out <- r
    spatialblocks(dat, corpars(object), function(idx, R) out[idx] <<- forwardsolve(t(chol(R)), r[idx]))
    out[dat$w <= 0] <- 0
    return(out)
  }
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
  if (cs$cortype >= 3){
    r <- stats::rnorm(n)
    spatialblocks(cs, corpars(object), function(idx, R) r[idx] <<- as.vector(t(chol(R)) %*% r[idx]))
    return(r)
  }
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
