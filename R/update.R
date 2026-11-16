#' Update a gcTMB fit, or set its parameters
#'
#' With \code{fixef}, \code{VarCorr} or \code{sigma}, \code{update} returns the model with some of
#' its parameters set to new values, e.g. as the true model of a simulation study (see
#' \code{\link{simulate.gcTMB}}); the model is not refitted. Without them, it refits the model with
#' the arguments in \code{...} changed, as \code{\link[stats]{update}}.
#'
#' The parameters that are not given keep their values, and the random effects are the conditional
#' modes at the new parameters. The new values are not estimates: the standard errors, the covariance
#' matrix of \code{vcov} and the confidence intervals are NA, and \code{\link{converged}} returns NA.
#' The log-likelihood is that of the data at the new parameters. For gaussian models without random
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
#' @param ... Without fixef, VarCorr and sigma: the arguments of \code{\link{gcTMB}} to change, as in
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
#'
#' @importFrom stats update
#' @export
update.gcTMB <- function(object, ..., fixef=NULL, VarCorr=NULL, sigma=NULL){
  if (is.null(fixef) && is.null(VarCorr) && is.null(sigma)) return(NextMethod())
  if (length(list(...)) > 0) stop("give either new parameter values (fixef, VarCorr, sigma) or arguments to refit the model, not both.")
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
#' Whether the fit converged to a proper maximum of the likelihood: the negative log-likelihood is
#' finite, the optimiser reports convergence, the Hessian of the parameters is positive definite, the
#' largest absolute gradient is below \code{gradtol}, and no standard deviation of the random effects
#' is close to zero (a singular fit: below \code{sdtol}, or with a standard error of its log that is
#' not finite or above 100, as in \code{summary}).
#'
#' @param object A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param gradtol Numeric. The largest acceptable absolute gradient of the negative log-likelihood,
#'   projected on the bounds of the coefficients.
#' @param sdtol Numeric. Standard deviations of the random effects below sdtol make the fit singular.
#'
#' @return TRUE or FALSE, with the attribute \code{reasons}, a character vector of the failed checks
#'   (empty for a converged fit). NA for a model with parameters set by \code{\link{update.gcTMB}}.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' converged(fit)
#'
#' @export
converged <- function(object, gradtol=1e-3, sdtol=1e-8){
  if (isTRUE(object$modified)) return(structure(NA, reasons="parameters set by update()"))
  s <- summary(object)
  sds <- unlist(lapply(s$varcor, function(vc) setNames(vc$sd, paste(vc$group, vc$labels))))
  small <- names(sds)[sds < sdtol]
  reasons <- c(
    if (!is.finite(as.numeric(s$logLik))) "the log-likelihood is not finite",
    if (s$convergence != 0) paste("the optimiser did not converge:", s$message),
    if (!isTRUE(s$pdHess)) "the Hessian is not positive definite",
    if (is.finite(s$maxgrad) && s$maxgrad > gradtol) paste("the largest absolute gradient is", format(s$maxgrad, digits=2)),
    if (!is.finite(s$maxgrad)) "the gradient is not finite",
    if (length(union(s$singular, small)) > 0) paste("singular fit:", paste(union(s$singular, small), collapse=", ")))
  structure(length(reasons) == 0, reasons=reasons)
}
