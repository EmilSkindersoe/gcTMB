#' Set the parameters of a gcTMB fit
#'
#' \code{update} returns the model with some of its parameters replaced by new values, without fitting
#' it again: the model keeps its structure (fform, rform, model, family, correlation and variance), and
#' \itemize{
#'   \item \code{fixef} sets the coefficients of b1-b5 of \code{fform} (including coefficients fixed in the fit);
#'   \item \code{VarCorr} the covariance matrices of the random effects of \code{rform};
#'   \item \code{sigma} the residual standard deviation, or the dispersion of the other families;
#'   \item \code{corpar} the parameters of the \code{correlation} structure (rho, or the range and nugget);
#'   \item \code{varpar} the parameters of the \code{variance} function.
#' }
#' The parameters that are not given keep their estimates, and the random effects are the conditional
#' modes at the new parameters. The result is not a fit: it is moved away from the maximum of the
#' likelihood, e.g. to be the true model of a simulation study (see \code{\link{simulate.gcTMB}}). The new
#' values are not estimates: the standard errors, the covariance matrix of \code{vcov} and the
#' confidence intervals are NA, and \code{\link{converged}} returns NA. The log-likelihood is that of the
#' data at the new parameters.
#'
#' To fit the model again with other components, e.g. a new \code{fform}, \code{rform},
#' \code{correlation}, \code{variance}, \code{model} or \code{family}, use \code{\link{refit}}; the other
#' arguments of \code{\link{gcTMB}} give an error here. To set parameters of a model with a new
#' component, refit it first and then set the parameters of the new fit, e.g.
#' \code{update(refit(fit, fform=b1 ~ trt), fixef=...)}, as the coefficients of the new fform are not
#' known before the refit.
#'
#' To fit the data as well as possible with a chosen value of a parameter, refit with that parameter
#' fixed instead, e.g. \code{refit(fit, correlation=corCAR1(0.8, form=~ time | ID, fixed=TRUE))} for
#' the autocorrelation, or \code{refit(fit, start=c(b5=2))} for the asymmetry, which is fixed by
#' default; the other parameters are then estimated given it, whereas setting it keeps their values.
#' For gaussian models without random effects, the residual standard deviation set by \code{sigma} is
#' the value given, without the small-sample correction of \code{\link{gcTMB}}.
#'
#' @param object A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param fixef Named numeric vector or NULL. New values of coefficients of b1-b5, named as in
#'   \code{\link{fixef}} (e.g. "b1:(Intercept)"), or by the parameter (e.g. "b1") for a parameter with
#'   a single coefficient. Coefficients fixed in the fit, such as the asymmetry b5 (1 by default), can
#'   be set too, e.g. \code{fixef=c(b5=4)}; they stay fixed. NULL (default) keeps all the coefficients.
#' @param VarCorr List or NULL. New covariance matrices of the random effects, named by the terms as
#'   in \code{\link{VarCorr}}, with the rows and columns in the same order; a single matrix for a model
#'   with one random effect term. The terms that are not given keep their covariance matrices, and NULL (default) keeps all of them.
#' @param sigma Numeric or NULL. The new residual standard deviation (gaussian, of log(y) for
#'   lognormal), or the dispersion parameter of the other families, as returned by
#'   \code{\link{sigma.gcTMB}}. NULL (default) keeps it.
#' @param corpar Numeric or NULL. New parameters of the correlation structure, as returned by
#'   \code{\link{corpars}}, for a model fitted with a \code{correlation} structure: the autocorrelation
#'   rho, between 0 and 1 for \code{\link{corCAR1}} and between -1 and 1 for \code{\link{corAR1}}; for
#'   \code{\link{corExp}} and \code{\link{corGaus}}, a named vector with the new range and/or nugget,
#'   e.g. \code{c(range=10, nugget=0.2)}. NULL (default) keeps them.
#' @param varpar Named numeric vector or NULL. New parameters of the variance function (see
#'   \code{\link{varFunc}}) on their natural scale, named as returned by \code{\link{varpars}}, e.g.
#'   \code{c(power=1)} or \code{c(const=2, power=0.5)} for \code{varConstPower}, or the ratio of a
#'   stratum for \code{varIdent}; the parameters not given keep their values. NULL (default) keeps them.
#' @param ... Not used; other arguments give an error that points to \code{\link{refit}}.
#'
#' @return A gcTMB object.
#' @seealso \code{\link{refit}} to fit the model again with other components.
#'
#' @examples
#' fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
#' # a larger asymptote, an asymmetric curve, and a larger variance of b1
#' vc <- VarCorr(fit)$Plant_ID
#' vc[1, 1] <- 2*vc[1, 1]
#' fit2 <- update(fit, fixef=c(b1=150, b5=4), VarCorr=list(Plant_ID=vc))
#' fixef(fit2)
#' VarCorr(fit2)
#' # a weaker autocorrelation of the residuals, the other parameters unchanged
#' fitc <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25,
#'               correlation=corCAR1(form=~ DAS | Plant_ID))
#' corpars(update(fitc, corpar=0.7))
#' # a residual standard deviation proportional to the fitted curve
#' fitv <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varPower())
#' varpars(update(fitv, varpar=c(power=1)))
#'
#' @importFrom stats update
#' @export
update.gcTMB <- function(object, fixef=NULL, VarCorr=NULL, sigma=NULL, corpar=NULL, varpar=NULL, ...){
  other <- names(match.call(expand.dots=FALSE)$...)
  if (length(match.call(expand.dots=FALSE)$...) > 0){
    stop("update() sets parameters of a fit (fixef, VarCorr, sigma, corpar, varpar)",
         if (length(other) > 0 && any(other != "")) paste(" and has no argument ", paste(other[other != ""], collapse=", "), sep=""),
         "; to fit the model again with other components, use refit(), e.g. refit(fit, fform=b1 ~ trt).", call.=FALSE)
  }
  if (is.null(fixef) && is.null(VarCorr) && is.null(sigma) && is.null(corpar) && is.null(varpar)){
    stop("give new parameter values in fixef, VarCorr, sigma, corpar or varpar; to fit the model again, use refit().", call.=FALSE)
  }
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

  if (!is.null(corpar)){
    rho <- corpar
    ct <- cortypeindex(object$correlation)
    if (ct == 0) stop("the model has no autocorrelation; fit it with a correlation structure, e.g. correlation=corCAR1(form=~ time | group).")
    if (ct >= 3){
      # spatial structures: the range and nugget, as named by corpars
      pn <- names(corpars(object))
      if (!is.numeric(rho) || is.null(names(rho)) || any(!(names(rho) %in% pn)) || any(!is.finite(rho))){
        stop("for ", class(object$correlation)[1], ", corpar must be a named vector with ", paste(pn, collapse=" and/or "), ", as returned by corpars().")
      }
      if ("range" %in% names(rho)){
        if (rho[["range"]] <= 0) stop("the range must be positive.")
        pl$cor_par[1] <- log(rho[["range"]])
      }
      if ("nugget" %in% names(rho)){
        if (rho[["nugget"]] <= 0 || rho[["nugget"]] >= 1) stop("the nugget must be between 0 and 1.")
        pl$cor_par[2] <- stats::qlogis(rho[["nugget"]])
      }
    } else {
      lo <- if (ct == 1) -1 else 0
      if (length(rho) != 1 || !is.finite(rho) || rho <= lo || rho >= 1) stop("corpar must be a single number between ", lo, " and 1 for ", class(object$correlation)[1], ".")
      # on the unconstrained scale of nlme
      pl$cor_par <- if (ct == 1) log((1 + rho)/(1 - rho)) else stats::qlogis(rho)
    }
  }

  if (!is.null(varpar)){
    vt <- vartypeindex(object$variance)
    if (vt == 0) stop("the model has no variance function; fit it with one, e.g. variance=varPower().")
    vn <- names(varpars(object))
    if (!is.numeric(varpar) || is.null(names(varpar)) || any(!(names(varpar) %in% vn)) || any(!is.finite(varpar))){
      stop("varpar must be a named vector with some of ", paste(vn, collapse=", "), ", as returned by varpars().")
    }
    G <- env$data$varng
    for (nm in names(varpar)){
      k <- match(nm, vn)
      v <- varpar[[nm]]
      if (vt == 1){
        if (k == 1) stop("the reference stratum ", nm, " of varIdent has the ratio 1.")
        if (v <= 0) stop("the ratios of varIdent must be positive.")
        pl$var_par[k - 1] <- log(v)
      } else if (vt == 3 && k <= G){
        if (v <= 0) stop("the constant of varConstPower must be positive.")
        pl$var_par[k] <- log(v)
      } else {
        pl$var_par[k] <- v
      }
    }
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
#' Whether the fit converged to a maximum of the likelihood. A fit that failed with an error and was
#' returned by \code{control=list(onerror="warning")} (see \code{\link{gcTMBcontrol}}) has not
#' converged, with its error as the reason. Otherwise: the negative log-likelihood is finite,
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
#' fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
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
  # a fit that failed with an error (control onerror="warning"): its error is the reason
  if (!is.null(object$fiterror)) return(structure(FALSE, reasons=paste("the fit failed:", object$fiterror)))
  reasons <- c(
    if (!is.finite(as.numeric(s$logLik))) "the log-likelihood is not finite",
    if (!stable && s$convergence != 0) paste("the optimiser did not converge:", s$message),
    if (!isTRUE(s$nosdreport) && !isTRUE(s$pdHess)) "the Hessian is not positive definite",
    if (!stable && is.finite(s$maxgrad) && s$maxgrad > gradtol) paste("the largest absolute gradient is", format(s$maxgrad, digits=2)),
    if (!is.finite(s$maxgrad)) "the gradient is not finite",
    if (singular && length(union(s$singular, small)) > 0) paste("singular fit:", paste(union(s$singular, small), collapse=", ")))
  structure(length(reasons) == 0, reasons=reasons)
}
