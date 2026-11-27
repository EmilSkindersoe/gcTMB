#' Refit a gcTMB model with other components
#'
#' \code{refit} fits a model again with some of its components changed, e.g. a new \code{fform},
#' \code{rform}, \code{correlation}, \code{variance}, \code{model} or \code{family}: the call of the
#' original fit with the arguments in \code{...} changed and all its other arguments reused (the
#' formula, the data, the other structures, \code{start}, \code{fix}, \code{control}, ...). The new fit is
#' a fit like any other, e.g. to compare with the original by \code{\link{anova.gcTMB}} or AIC. To change
#' the parameters of a fit without fitting it again, e.g. for the true model of a simulation, use
#' \code{\link{update.gcTMB}}.
#'
#' The starting values of the new fit are chosen by \code{startfrom}:
#' \itemize{
#'   \item \code{"data"} (default): from the data, as for a new call of \code{\link{gcTMB}} (see
#'     \code{\link{gcTMBstart}});
#'   \item \code{"fit"}: from the estimates of the fit being refitted, where they carry over, as a
#'     \code{start} list (see \code{\link{gcTMB}}): the coefficients with the same name, growth model and
#'     link (e.g. the intercepts when covariates are added to \code{fform}, or all coefficients when
#'     only the correlation, the variance function or the random effects change); the variances and
#'     covariances of the random effects of the parameters on the same grouping factor (e.g. that of b1
#'     when b2 is added to \code{b1 ~ 1 | Plant_ID}); the residual standard deviation or dispersion with
#'     the same family and variance function; and the parameters of the correlation structure and of the
#'     variance function when the refit keeps them (and they are not fixed). The rest starts from the
#'     data, e.g. the new coefficients of covariates, all coefficients of a new growth model, and the
#'     correlations of a random effect with a new one. With a new \code{rform}, the residual standard
#'     deviation of the fit may not suit the new random effects, but the optimiser moves it with them.
#'     A \code{start}, \code{correlation} or \code{variance} given in \code{...} is used as given. The
#'     call of the new fit shows the starting values used.
#' }
#' The data are found again by their name in the call, so a data frame changed since the original fit
#' gives a fit to the changed data.
#'
#' The package lme4 has a \code{refit} function of its own, for a fit to a new response. For objects
#' that are not gcTMB fits, \code{refit} calls that of lme4; if lme4 is attached after gcTMB, its
#' \code{refit} masks this one, and \code{gcTMB::refit} refits a gcTMB model.
#'
#' @param object A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param ... The arguments of \code{\link{gcTMB}} to change, named, e.g. \code{fform=b1 + b2 ~ trt},
#'   \code{correlation=corCAR1(form=~ DAS | Plant_ID)} or \code{model="gompertz"}; NULL removes a
#'   component, e.g. \code{correlation=NULL}. An unnamed first argument changes the formula, as
#'   \code{\link[stats]{update.formula}}. For other objects, the arguments of the refit method of lme4.
#' @param startfrom Character string. The starting values: "data" (default) for those of
#'   \code{\link{gcTMBstart}}, or "fit" for the estimates of \code{object} where they carry over, see
#'   Details.
#'
#' @return A gcTMB object.
#' @seealso \code{\link{update.gcTMB}} to set the parameters of a fit without refitting it.
#'
#' @examples
#' fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
#' # another growth model, and covariates on b1 and b2, the other arguments reused
#' fitg <- refit(fit, model="gompertz")
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitt <- refit(fit, fform=b1 + b2 ~ trt)
#' anova(fit, fitt)
#' # an autocorrelation structure with a random asymptote only
#' fitc <- refit(fit, rform=b1 ~ 1 | Plant_ID, correlation=corCAR1(form=~ DAS | Plant_ID))
#' # a refit that starts from the estimates of the fit it refits
#' fitc2 <- refit(fitc, fform=b1 + b2 ~ trt, startfrom="fit")
#' fitc2$call
#' @export
refit <- function(object, ...) UseMethod("refit")

#' @rdname refit
#' @export
refit.default <- function(object, ...){
  if (!requireNamespace("lme4", quietly=TRUE)) stop("no refit method for an object of class ", class(object)[1], ".", call.=FALSE)
  getExportedValue("lme4", "refit")(object, ...)
}

#' @rdname refit
#' @export
refit.gcTMB <- function(object, ..., startfrom=c("data", "fit")){
  startfrom <- match.arg(startfrom)
  # the call of the fit with the arguments in ... changed, as stats::update (an unnamed first argument
  # updates the formula)
  extras <- match.call(expand.dots=FALSE)$...
  cl <- stats::getCall(object)
  if (length(extras) > 0){
    nm <- names(extras)
    if (is.null(nm)) nm <- rep("", length(extras))
    if (nm[1] == ""){
      cl$form <- stats::update.formula(stats::formula(eval(cl$form, parent.frame())), eval(extras[[1]], parent.frame()))
      extras <- extras[-1]
      nm <- nm[-1]
    }
    if (any(nm == "")) stop("the arguments to change must be named, e.g. refit(fit, fform=b1 ~ trt).")
    bad <- intersect(nm, c("fixef", "VarCorr", "sigma", "corpar", "varpar"))
    if (length(bad) > 0) stop("refit() fits the model again with other components; to set ", paste(bad, collapse=", "),
                              " without refitting, use update().")
    for (a in nm) cl[a] <- list(extras[[a]])
  }
  if (startfrom == "fit") cl <- fitstartcall(object, cl, names(extras), parent.frame())
  eval(cl, parent.frame())
}

# the refit call cl of refit(startfrom="fit"): starting values from the fit object for what the call
# does not set itself (changed, the names of the arguments given in refit): the coefficients with the
# same name, growth model and link as in the fit (the others NA, from gcTMBstart), unless start is
# given; the parameters of the correlation structure and the variance function, unless they are given
# or fixed. The random effect variances and the residual standard deviation start as in gcTMB.
fitstartcall <- function(object, cl, changed, env){
  if (!("start" %in% changed)){
    dargs <- as.list(cl)[-1]
    dargs <- dargs[intersect(names(dargs), names(formals(gcTMBdata)))]
    dd <- do.call(gcTMBdata, lapply(dargs, eval, envir=env))
    newnames <- coefnames(dd$Xs)
    pv <- object$start
    pv[!object$fix] <- object$sdrl$par.fixed[seq_len(sum(!object$fix))]
    st <- rep(NA_real_, length(newnames))
    if (identical(dd$model, object$model)){
      m <- match(newnames, object$coefnames)
      par <- sub(":.*", "", newnames)
      samelink <- vapply(par, function(p) identical(dd$link[[p]], object$link[[p]]), logical(1))
      ok <- !is.na(m) & samelink
      st[ok] <- pv[m[ok]]
    }
    # the covariance matrices of the random effects: the entries of the parameters on the same grouping
    # factor; sigma with the same family and variance function
    vc <- refitvarcorr(object, dd$re$terms)
    sg <- NULL
    if (!(dd$family %in% c("poisson", "binomial")) && identical(dd$family, object$family) && !("variance" %in% changed)){
      sg <- exp(object$obj$env$parList(object$opt$par)$log_sigma)
    }
    fe <- if (any(!is.na(st))) unname(st)
    if (length(vc) > 0 || !is.null(sg)){
      cl$start <- c(if (!is.null(fe)) list(fixef=fe), if (length(vc) > 0) list(VarCorr=vc), if (!is.null(sg)) list(sigma=sg))
    } else if (!is.null(fe)) cl$start <- fe
  }
  if (!("correlation" %in% changed)){
    cc <- corstartcall(object)
    if (!is.null(cc)) cl$correlation <- cc
  }
  if (!("variance" %in% changed)){
    vc <- varstartcall(object)
    if (!is.null(vc)) cl$variance <- vc
  }
  cl
}

# starting covariance matrices of the random effect terms of a refit from those of the fit object: for
# each new term, the variances and covariances of the parameters (and covariates of random slopes) that
# a term of the fit has on the same grouping factor, NA otherwise; terms without any are left out
refitvarcorr <- function(object, terms){
  if (length(terms) == 0 || length(object$reterms) == 0) return(list())
  old <- VarCorr(object)
  oldgroup <- vapply(object$reterms, function(tm) tm$group, character(1))
  tn <- make.unique(vapply(terms, function(tm) tm$group, character(1)))
  out <- list()
  for (t in seq_along(terms)){
    lab <- relabels(terms[[t]])
    np <- length(lab)
    S <- matrix(NA_real_, np, np, dimnames=list(lab, lab))
    for (o in which(oldgroup == terms[[t]]$group)){
      So <- old[[o]]
      idx <- match(lab, rownames(So))
      have <- which(!is.na(idx))
      if (length(have) > 0) S[have, have] <- So[idx[have], idx[have]]
    }
    if (any(is.finite(S))) out[[tn[t]]] <- unname(S)
  }
  out
}

# a call of the constructor of the correlation structure of a fit with its estimates as the starting
# values; NULL without an estimated structure
corstartcall <- function(object){
  cs <- object$correlation
  ct <- cortypeindex(cs)
  if (ct == 0 || isTRUE(attr(cs, "fixed"))) return(NULL)
  cp <- unname(corpars(object))
  form <- attr(cs, "formula")
  if (ct <= 2) return(as.call(list(as.name(c("corAR1", "corCAR1")[ct]), value=cp, form=form)))
  as.call(list(as.name(c("corExp", "corGaus")[ct - 2]), value=cp, form=form, nugget=isTRUE(attr(cs, "nugget")),
               metric=if (is.null(attr(cs, "metric"))) "euclidean" else attr(cs, "metric")))
}

# a call of the constructor of the variance function of a fit with its estimates as the starting
# values; NULL without a variance function or with fixed parameters
varstartcall <- function(object){
  vf <- object$variance
  vt <- vartypeindex(vf)
  if (vt == 0) return(NULL)
  fx <- if (vt == 3) c(attr(vf$const, "fixed"), attr(vf$power, "fixed")) else attr(vf, "fixed")
  if (length(fx) > 0) return(NULL)
  vp <- varpars(object)
  form <- attr(vf, "formula")
  lev <- object$varlevels
  if (vt == 1) return(as.call(list(as.name("varIdent"), value=vp[-1], form=form)))
  G <- length(lev)
  named <- function(v) if (G == 1) unname(v) else setNames(unname(v), lev)
  if (vt == 2) return(as.call(list(as.name("varPower"), value=named(vp), form=form)))
  as.call(list(as.name("varConstPower"), const=named(vp[seq_len(G)]), power=named(vp[G + seq_len(G)]), form=form))
}
