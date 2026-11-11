#' Time to reach a level of the asymptote
#'
#' Estimates the times at which growth curves reach a percentage of their asymptote, or a given
#' response, with standard errors by the delta method and automatic differentiation.
#'
#' @param x A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param respLev Numeric vector. Levels in percent of the asymptote b1, between 0 and 100 (default 50): the time t p at which the curve reaches p percent of b1, b2 + b3 Finv(p/100) with the inverse Finv of F of the model (see \code{\link{gcTMB}}). For the logistic model, t 50 is b2. For \code{type="absolute"}, values of the response.
#' @param linfct Numeric matrix, list of matrices, or NULL. Linear functions to identify the parameters of the growth curve, with 5 rows (b1-b5) and one column for each model coefficient, or a list of such matrices; default is the identity matrix, for models without covariates in fform. Each matrix is paired with an entry of respLev.
#' @param newdata Data frame or NULL. The values of the covariates in fform, one row for each growth curve, as an alternative to linfct; the time is not needed. With grouping variables of the random effects, the curves of fitted groups, see Details. The times are computed at each level of respLev for each row, named by the covariate and group values (or by the row names of newdata, if set).
#' @param random Logical. If TRUE, the times of rows of newdata with fitted groups are conditional on their random effects, and for the other random effect terms, the variance-covariance matrix of the random effects is included in the standard error estimation. If FALSE, the times are for the population curves.
#' @param residual Logical. If TRUE, the residual variance is included in the standard error estimation (gaussian family only)
#' @param type Character string. "relative" (default) for levels relative to the asymptote, or "absolute" for the times at which the curves reach the responses in respLev, on the scale of the response (proportions for binomial, beta-binomial and beta). Responses outside the range of a curve, between zero and the asymptote, give NA, with a warning.
#'
#' @details Without random effects, or with \code{random=FALSE}, the times are those of the
#' population curves, with all random effects at zero, and their standard errors include the
#' uncertainty of the fixed effects only. With \code{random=TRUE}, the standard errors also include the
#' variance of the random effects, the variation of the time between groups.
#'
#' If a row of newdata gives a level of a grouping factor that is in the fitted data, the time of
#' that row is for the curve of the group, with its predicted random effects (the conditional
#' modes), as in \code{\link{predict.gcTMB}}. These are shrunken towards the population curve, more so
#' for groups with few observations. Their standard errors include the uncertainty of the fixed
#' effects and the prediction error of the random effects, from the joint precision matrix of
#' \code{\link[TMB]{sdreport}}. Random effect terms without a grouping variable in newdata, or with a
#' level that is not in the fitted data, are treated as above. Profile intervals are not available
#' for the times of fitted groups, see \code{\link{confint.timetoTMB}}.
#'
#' @return A timetoTMB object: a list with the estimated times (\code{coef}), their covariance
#'   matrix (\code{vcov}), and what \code{confint} and \code{profile} need.
#' @seealso \code{\link{confint.timetoTMB}} for confidence intervals, \code{\link{profile.timetoTMB}}
#' @examples
#' fit <- gcTMB(height ~ DAS, data=FabaHeights25, model="logistic")
#' # days after sowing at which the plants reach 10%, 50% and 90% of their final height
#' tt <- timeto(fit, respLev=c(10, 50, 90))
#' tt
#' confint(tt)
#' # the day at which they reach a height of 50
#' timeto(fit, respLev=50, type="absolute")
#'
#' # a curve for each treatment level
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fit2 <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25)
#' timeto(fit2, respLev=90, newdata=data.frame(trt=levels(FabaHeights25$trt)))
#' @export
#'
timeto <- function(x, respLev=50, linfct=NULL, newdata=NULL, random=TRUE, residual=FALSE, type=c("relative", "absolute")){
  type <- match.arg(type)
  absolute <- type == "absolute"
  if (!is.null(newdata)){
    if (!is.null(linfct)) stop("give either linfct or newdata, not both.")
    # every level of respLev for every curve
    nl <- curvelinfct(x, newdata)
    rows <- rep(seq_along(nl), length(respLev))
    linfct <- rep(nl, length(respLev))
    respLev <- rep(respLev, each=length(nl))
  }
  if (is.null(linfct)) linfct <- diag(5)
  if (!is.list(linfct)) linfct <- list(linfct)
  if (length(respLev) != length(linfct)){
    if (length(respLev) == 1){
      respLev <- rep(respLev, length(linfct))
    }
    if (length(linfct) == 1){
      linfct <- rep(linfct, length(respLev))
    }
    if (length(respLev) != length(linfct)){
      stop("respLev and linfct must have the same length or one of them must be of length 1.")
    }
  }
  if (is.null(names(linfct))) names(linfct) <- paste("lf", 1:length(linfct), sep="")
  ncf <- length(x$start)
  for (L in linfct){
    if (!is.matrix(L) || nrow(L) != 5 || ncol(L) != ncf){
      stop("each matrix in linfct must have 5 rows (b1-b5) and one column for each of the ", ncf,
           " model coefficients: ", paste(x$coefnames, collapse=", "), ".",
           if (is.null(newdata) && ncf > 5) " Alternatively, give the covariate values of the curves in newdata.")
    }
  }
  if (!absolute && any(respLev <= 0 | respLev >= 100)) stop("respLev must be between 0 and 100 percent of the asymptote.")

  if (residual == TRUE && x$family != "gaussian"){
    stop("residual=TRUE is only available for the gaussian family.")
  }
  rvar <- if (residual) residpar(x)^2 else 0

  pv <- x$start
  pv[!x$fix] <- x$sdrl$par.fixed[seq_len(sum(!x$fix))]
  nb <- sum(!x$fix)
  terms <- x$reterms
  uhat <- if (length(terms) > 0) x$sdrl$par.random else numeric(0)
  nu <- length(uhat)

  # random effects of each time: the weights W of the fitted random effects on b1-b5, for the rows of
  # newdata with a fitted group (random=TRUE), and the terms whose variance is added for a new group
  # (random=TRUE, no fitted level)
  nt <- length(linfct)
  W <- rep(list(matrix(0, 5, nu)), nt)
  newterm <- matrix(random, nt, length(terms))
  # columns of each term for a new group (intercept, covariates of slopes), for each time
  loads <- lapply(terms, function(tm) reloads(tm, matrix(as.numeric(tm$cnms == "(Intercept)"), nt, length(tm$cnms), byrow=TRUE)))
  if (!is.null(newdata) && random && length(terms) > 0){
    rl <- relevels(x, newdata, skip=TRUE)
    W <- lapply(rows, function(r) reweights(x, rl, r))
    newterm <- rl$newgroup[rows, , drop=FALSE]
    loads <- lapply(seq_along(terms), function(t) reloads(terms[[t]], rl$X[[t]][rows, , drop=FALSE]))
    if (any(is.na(unlist(loads)))) stop("newdata has missing values in the covariates of the random slopes.")
  } else if (random && any(unlist(lapply(terms, function(tm) tm$cnms)) != "(Intercept)")){
    stop("the variance of the times depends on the covariates of the random slopes: give their values in newdata, or use random=FALSE.")
  }
  cond <- vapply(W, function(w) any(w != 0), logical(1))

  mod <- modelindex(x$model)
  lnk <- linkindex(x$link)

  # gradient of each time w.r.t. the free coefficients, the fitted random effects, the random effects
  # of a new group for each term, and the residual
  tl <- numeric(nt)
  ql <- numeric(nt)
  grad <- matrix(0, nt, nb + nu + sum(vapply(terms, function(tm) length(tm$pars), integer(1))) + 1)

  for (i in seq_len(nt)){
    lpv <- as.vector(linfct[[i]] %*% pv)
    u <- as.vector(W[[i]] %*% uhat)
    # fraction of the asymptote, or the response for absolute times
    q <- if (absolute) respLev[i] else respLev[i]/100
    ql[i] <- q
    qf <- timeq(x, q, lpv + u, absolute)
    if (!is.finite(qf) || qf <= 0 || qf >= 1){
      tl[i] <- NA
      grad[i, ] <- NA
      next
    }

    # time and its gradient w.r.t. (b1-b5, u1-u5, r) from the gcTime TMB template
    tobj <- MakeADFun(list(model="gcTime", mod=mod, lnk=as.integer(lnk), absolute=as.integer(absolute), q=q),
                      list(b=lpv, u=u, r=0),
                      ADreport=TRUE, DLL="gcTMB_TMBExports", silent=TRUE)
    tl[i] <- tobj$fn()
    tgrad <- tobj$gr()
    gu <- tgrad[6:10]
    gnew <- unlist(lapply(seq_along(terms), function(t) gu[terms[[t]]$pars] * loads[[t]][i, ] * newterm[i, t]))
    grad[i, ] <- c((tgrad[1:5] %*% linfct[[i]])[!x$fix], gu %*% W[[i]], gnew, tgrad[11])
  }

  # covariance of the free coefficients and the fitted random effects: jointly for times of fitted
  # groups, otherwise of the coefficients only
  cbu <- if (any(cond)) jointcov(x) else as.matrix(bdiag(fixedcov(x)[seq_len(nb), seq_len(nb), drop=FALSE], matrix(0, nu, nu)))
  cvc <- as.matrix(bdiag(c(list(cbu), lapply(recov(x), function(vc) vc$cov), list(matrix(rvar)))))
  tvar <- grad %*% cvc %*% t(grad)
  tn <- paste(if (absolute) paste("y=", respLev, sep="") else respLev, names(linfct), sep=":")
  if (anyNA(tl)){
    warning("the curves do not reach the response for ", paste(tn[is.na(tl)], collapse=", "), "; the times are NA.", call.=FALSE)
  }
  names(tl) <- tn
  colnames(tvar) <- rownames(tvar) <- tn
  # the fit, linear functions, random effect weights and levels are kept for profile and confint
  names(ql) <- names(tl)
  res <- list(estimates=tl, variance=tvar, fit=x, linfct=linfct, q=ql, type=type, random=random, residual=residual,
              uweights=W, conditional=setNames(cond, names(tl)),
              newvar=setNames(apply(newterm, 1, any) & length(terms) > 0, names(tl)))
  class(res) <- "timetoTMB"
  return(res)
}

# fraction of the asymptote at a time with level q (absolute=FALSE) or response q (absolute=TRUE), for
# the parameters lpv (b1-b5 on the scale of the linear predictors) of the curve
timeq <- function(fit, q, lpv, absolute){
  if (!absolute) return(q)
  b <- linkpars(linkindex(fit$link), lpv)
  q / b[1]
}

# weights of the fitted random effects u on b1-b5 for row r of the random effect design rl of newdata
reweights <- function(object, rl, r){
  w <- matrix(0, 5, length(object$sdrl$par.random))
  uoff <- 0
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    for (k in seq_along(tm$pars)) w[tm$pars[k], uoff + (k - 1)*nl + seq_len(nl)] <- as.vector(rl$Zl[[t]][[tm$zcol[k]]][r, ])
    uoff <- uoff + nl*length(tm$pars)
  }
  w
}

# joint covariance of the free coefficients b1-b5 (all fixed parameters with all=TRUE) and the
# prediction errors of the random effects, from the joint precision matrix of sdreport; the likelihood
# is unchanged, so the Hessian of the fit is reused
jointcov <- function(object, all=FALSE){
  env <- object$obj$env
  H <- if (isTRUE(object$sdrl$pdHess)) solve(object$sdrl$cov.fixed)
  sdr <- sdreport(object$obj, par.fixed=object$opt$par, hessian.fixed=H, getJointPrecision=TRUE)
  id <- c(if (all) seq_along(env$par)[-env$random] else seq_len(sum(!object$fix)), env$random)
  C <- tryCatch(solve(sdr$jointPrecision), error=function(e) NULL)
  if (is.null(C)) stop("the joint precision matrix of the coefficients and random effects is singular; check the convergence of the fit.")
  as.matrix(C)[id, id, drop=FALSE]
}

# linear functions of the curves in the rows of newdata: the rows of the design matrices of b1-b5,
# named by the values of the covariates in fform and the grouping variables in newdata, or by the
# row names of newdata if set
curvelinfct <- function(x, newdata){
  if (!is.data.frame(newdata) || nrow(newdata) == 0) stop("newdata must be a data frame with at least one row.")
  vars <- unique(unlist(lapply(x$fterms, function(ft) all.vars(ft$terms))))
  miss <- setdiff(vars, names(newdata))
  if (length(miss) > 0) stop("newdata needs the covariates of fform: ", paste(miss, collapse=", "), ".")
  Xs <- lapply(x$fterms, function(ft) designmatrix(ft, newdata))
  if (any(sapply(Xs, anyNA))) stop("newdata has missing values in the covariates of fform.")
  npi <- rep(1:5, sapply(Xs, ncol))
  L <- lapply(seq_len(nrow(newdata)), function(r){
    L <- matrix(0, 5, length(npi))
    for (i in 1:5) L[i, npi == i] <- Xs[[i]][r, ]
    L
  })
  # grouping variables of the random effects in newdata
  gvars <- unique(unlist(lapply(x$reterms, function(tm) all.vars(tm$formula))))
  vars <- c(vars, intersect(setdiff(gvars, vars), names(newdata)))
  if (.row_names_info(newdata) > 0){
    names(L) <- rownames(newdata)
  } else if (length(vars) > 0){
    names(L) <- do.call(paste, c(lapply(vars, function(v) paste(v, newdata[[v]], sep="=")), sep=","))
  } else {
    names(L) <- seq_len(nrow(newdata))
  }
  L
}
