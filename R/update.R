#' Update a gcTMB fit, or set its parameters
#'
#' With \code{fixef}, \code{VarCorr}, \code{sigma} or \code{rho}, \code{update} returns the model with some of
#' its parameters set to new values, e.g. as the true model of a simulation study (see
#' \code{\link{simulate.gcTMB}}); the model is not refitted. Without them, it refits the model with
#' the arguments in \code{...} changed, as \code{\link[stats]{update}}.
#'
#' Every parameter of the model can be set: the coefficients of b1-b5 (\code{fixef}, including
#' coefficients fixed in the fit), the covariance matrices of the random effects (\code{VarCorr}), the
#' residual standard deviation or dispersion (\code{sigma}) and the autocorrelation (\code{rho}),
#' i.e. all the entries of \code{coef}. The parameters that are not given keep their values, and the random effects are the conditional
#' modes at the new parameters. The new values are not estimates: the standard errors, the covariance
#' matrix of \code{vcov} and the confidence intervals are NA, and \code{\link{converged}} returns NA.
#' The log-likelihood is that of the data at the new parameters. To fit the data as well as possible
#' with a chosen value of a parameter, refit the model with that parameter fixed instead, e.g.
#' \code{correlation=corCAR1(0.8, form=~ day | plant, fixed=TRUE)} for the autocorrelation, or
#' \code{start=c(b5=2)} for the asymmetry, which is fixed by default; the other parameters are then estimated given it, whereas
#' \code{update} keeps their values. For gaussian models without random
#' effects, the residual standard deviation is the value given, without the small-sample correction
#' of \code{\link{gcTMB}}.
#'
#' @param object A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param fixef Named numeric vector or NULL. New values of coefficients of b1-b5, named as in
#'   \code{\link{fixef}} (e.g. "b1:(Intercept)"), or by the parameter (e.g. "b1") for a parameter with
#'   a single coefficient. Coefficients fixed in the fit, such as the asymmetry b5 (1 by default), can
#'   be set too, e.g. \code{fixef=c(b5=4)}; they stay fixed.
#' @param VarCorr List or NULL. New covariance matrices of the random effects, named by the terms as
#'   in \code{\link{VarCorr}}, with the rows and columns in the same order; a single matrix for a model
#'   with one random effect term. The terms that are not given keep their covariance matrices.
#' @param sigma Numeric or NULL. The new residual standard deviation (gaussian, of log(y) for
#'   lognormal), or the dispersion parameter of the other families, as returned by
#'   \code{\link{sigma.gcTMB}}.
#' @param rho Numeric or NULL. The new autocorrelation rho of the residuals, as returned by
#'   \code{\link{corpars}}, for a model fitted with a \code{correlation} structure: between 0 and 1 for
#'   \code{\link{corCAR1}}, between -1 and 1 for \code{\link{corAR1}}.
#' @param ... Without fixef, VarCorr, sigma and rho: the arguments of \code{\link{gcTMB}} to change, as in
#'   \code{\link[stats]{update}}.
#'
#' @return A gcTMB object.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' # larger asymptote, an asymmetric curve, and a larger variance of b1
#' vc <- VarCorr(fit)$plant_ID
#' vc[1, 1] <- 2*vc[1, 1]
#' fit2 <- update(fit, fixef=c(b1=150, b5=4), VarCorr=list(plant_ID=vc))
#' fixef(fit2)
#' VarCorr(fit2)
#' # refit with another growth model
#' fitg <- update(fit, model="gompertz")
#' # a weaker autocorrelation of the residuals, the other parameters unchanged
#' fitc <- gcTMB(height ~ DAS, rform=b1 ~ 1 | plant_ID, data=FabaHeights25,
#'               correlation=corCAR1(form=~ DAS | plant_ID))
#' corpars(update(fitc, rho=0.7))
#'
#' @importFrom stats update
#' @export
update.gcTMB <- function(object, ..., fixef=NULL, VarCorr=NULL, sigma=NULL, rho=NULL){
  if (is.null(fixef) && is.null(VarCorr) && is.null(sigma) && is.null(rho)) return(NextMethod())
  if (length(list(...)) > 0) stop("give either new parameter values (fixef, VarCorr, sigma, rho) or arguments to refit the model, not both.")
  env <- object$obj$env
  pl <- env$parList(object$opt$par)
  start <- object$start

  if (!is.null(fixef)){
    j <- matchcoefs(object, fixef)
    npi <- as.integer(sub("^b([1-5]):.*", "\\1", object$coefnames))
    pos <- ave(npi, npi, FUN=seq_along)
    for (k in seq_along(j)){
      pl[[paste("b", npi[j[k]], sep="")]][pos[j[k]]] <- unname(fixef[k])
      start[j[k]] <- unname(fixef[k])
    }
  }

  if (!is.null(VarCorr)){
    tn <- retermnames(object)
    if (is.matrix(VarCorr)){
      if (length(tn) != 1) stop("VarCorr must be a list of covariance matrices named by the random effect terms: ", paste(tn, collapse=", "), ".")
      VarCorr <- setNames(list(VarCorr), tn)
    }
    if (!is.list(VarCorr) || is.null(names(VarCorr)) || any(!(names(VarCorr) %in% tn))){
      stop("VarCorr must be a list of covariance matrices named by the random effect terms: ", paste(tn, collapse=", "), ".")
    }
    toff <- 0
    for (t in seq_along(object$reterms)){
      np <- length(object$reterms[[t]]$pars)
      nt <- np + np*(np - 1)/2
      S <- VarCorr[[tn[t]]]
      if (!is.null(S)) pl$theta[toff + seq_len(nt)] <- covtheta(as.matrix(S), np, tn[t])
      toff <- toff + nt
    }
  }

  if (!is.null(sigma)){
    if (object$family %in% c("poisson", "binomial")) stop("the ", object$family, " family has no residual standard deviation or dispersion parameter.")
    if (length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) stop("sigma must be a single positive number.")
    pl$log_sigma <- log(sigma)
  }

  if (!is.null(rho)){
    ct <- cortypeindex(object$correlation)
    if (ct == 0) stop("the model has no autocorrelation; fit it with a correlation structure, e.g. correlation=corCAR1(form=~ time | group).")
    lo <- if (ct == 1) -1 else 0
    if (length(rho) != 1 || !is.finite(rho) || rho <= lo || rho >= 1) stop("rho must be a single number between ", lo, " and 1 for ", class(object$correlation)[1], ".")
    # on the unconstrained scale of nlme
    pl$cor_par <- if (ct == 1) log((1 + rho)/(1 - rho)) else stats::qlogis(rho)
  }

  # a new objective function: the one of the fit is shared with the original object
  rnd <- !is.null(env$random)
  obj <- MakeADFun(env$data, pl, map=env$map, random=if (rnd) "u", DLL="gcTMB_TMBExports", silent=TRUE)
  par <- obj$par
  nll <- obj$fn(par)
  obj$env$last.par.best <- obj$env$last.par
  # the values at the new parameters, e.g. the conditional modes of the random effects and the
  # correlations of the random effects; no standard errors
  sdr <- suppressWarnings(sdreport(obj, par.fixed=par, hessian.fixed=diag(length(par)), getReportCovariance=FALSE))
  sdr$cov.fixed[] <- NA
  sdr$pdHess <- FALSE
  ssdr <- suppressWarnings(summary(sdr))
  ssdr[, 2] <- NA

  object$obj <- obj
  object$start <- start
  object$estimates <- par
  object$opt$par <- par
  object$opt$objective <- nll
  object$opt$message <- "parameters set by update()"
  object$sdrl <- sdr
  object$ssdrl <- ssdr
  object$rdf <- NULL
  object$modified <- TRUE
  object$call <- match.call()
  object
}

# indices in coefnames of the coefficients named in x: full names, or the parameter (e.g. "b1") for
# a parameter with a single coefficient
matchcoefs <- function(object, x){
  if (is.null(names(x)) || any(names(x) == "")) stop("fixef must be a named vector, e.g. c(b1=150).")
  cn <- object$coefnames
  par <- sub(":.*", "", cn)
  vapply(names(x), function(nm){
    j <- which(cn == nm)
    if (length(j) == 0) j <- which(par == nm)
    if (length(j) == 0) stop("unknown coefficient ", nm, "; the coefficients are ", paste(cn, collapse=", "), ".", call.=FALSE)
    if (length(j) > 1) stop("parameter ", nm, " has several coefficients; give their full names: ", paste(cn[j], collapse=", "), ".", call.=FALSE)
    j
  }, integer(1))
}

# variance parameters theta of a term with np random coefficients from its covariance matrix S: the
# log standard deviations and the unconstrained parameters of the correlation matrix (see corrmatrix)
covtheta <- function(S, np, name){
  if (!all(dim(S) == np)) stop("the covariance matrix of ", name, " must be ", np, " x ", np, ".", call.=FALSE)
  if (!isSymmetric(unname(S)) || any(diag(S) <= 0)) stop("the covariance matrix of ", name, " must be symmetric with positive variances.", call.=FALSE)
  sd <- sqrt(diag(S))
  if (np == 1) return(log(sd))
  Lc <- tryCatch(t(chol(cov2cor(S))), error=function(e) stop("the covariance matrix of ", name, " must be positive definite.", call.=FALSE))
  L <- Lc / diag(Lc)
  c(log(sd), t(L)[upper.tri(L)])
}

#' Convergence of a gcTMB fit
#'
#' Whether the fit converged to a maximum of the likelihood: the negative log-likelihood is finite,
#' the optimiser reports convergence, the Hessian of the parameters is positive definite (not checked
#' for fits with \code{control=list(sdreport=FALSE)}, see \code{\link{gcTMBcontrol}}), and the
#' largest absolute gradient is below \code{gradtol}. By default, fits with standard deviations of
#' the random effects close to zero count as converged; with \code{singular=TRUE} such singular fits
#' do not (a standard deviation below \code{sdtol}, or with a standard error of its log that is not
#' finite or above 100, as in \code{summary}). Note that the Hessian of a singular fit is often not
#' positive definite.
#'
#' A fit with \code{control=list(restart=TRUE)} (see \code{\link{gcTMBcontrol}}) that was optimised a
#' second time from its estimates, and whose log-likelihood the restart improved by less than
#' \code{lltol}, is at the maximum: the convergence message of the optimiser and the gradient are
#' then not checked. With few groups and several correlated random effects, the likelihood is often
#' too flat in some directions for nlminb to report convergence or reach a gradient below gradtol,
#' although the maximum is reached.
#'
#' @param object A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param gradtol Numeric. The largest acceptable absolute gradient of the negative log-likelihood,
#'   projected on the bounds of the coefficients.
#' @param singular Logical. If TRUE, singular fits, with a standard deviation of the random effects
#'   close to zero, count as not converged; FALSE (default) accepts them.
#' @param sdtol Numeric. With singular=TRUE, standard deviations of the random effects below sdtol
#'   make the fit singular.
#' @param lltol Numeric. For a fit restarted by \code{control=list(restart=TRUE)}: if the restart
#'   improved the log-likelihood by less than lltol, the optimiser message and the gradient are not
#'   checked.
#'
#' @return TRUE or FALSE, with the attribute \code{reasons}, a character vector of the failed checks
#'   (empty for a converged fit). NA for a model with parameters set by \code{\link{update.gcTMB}}.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' converged(fit)
#'
#' @export
converged <- function(object, gradtol=1e-3, singular=FALSE, sdtol=1e-8, lltol=1e-6){
  if (isTRUE(object$modified)) return(structure(NA, reasons="parameters set by update()"))
  s <- summary(object)
  sds <- unlist(lapply(s$varcor, function(vc) setNames(vc$sd, paste(vc$group, vc$labels))))
  small <- names(sds)[sds < sdtol]
  # a restart (control restart=TRUE) that hardly changed the log-likelihood: the fit is at the maximum
  stable <- !is.null(object$restart) && is.finite(object$restart$improvement) && object$restart$improvement < lltol
  reasons <- c(
    if (!is.finite(as.numeric(s$logLik))) "the log-likelihood is not finite",
    if (!stable && s$convergence != 0) paste("the optimiser did not converge:", s$message),
    if (!isTRUE(s$nosdreport) && !isTRUE(s$pdHess)) "the Hessian is not positive definite",
    if (!stable && is.finite(s$maxgrad) && s$maxgrad > gradtol) paste("the largest absolute gradient is", format(s$maxgrad, digits=2)),
    if (!is.finite(s$maxgrad)) "the gradient is not finite",
    if (singular && length(union(s$singular, small)) > 0) paste("singular fit:", paste(union(s$singular, small), collapse=", ")))
  structure(length(reasons) == 0, reasons=reasons)
}
