#' Set control parameters for gcTMB
#'
#' @param optimiser Character string. The optimiser: only "nlminb" (\code{\link{nlminb}}, with the
#'   gradients of TMB) is available, for the fit, the profiles and the bootstrap.
#' @param iter.max Integer. Maximum number of iterations of the optimiser, see \code{\link{nlminb}}.
#' @param eval.max Integer. Maximum number of evaluations of the objective function, see \code{\link{nlminb}}.
#' @param smallsample Logical. If TRUE (default), gaussian and lognormal models without random effects and without autocorrelation use the unbiased
#'   residual variance with n - p residual degrees of freedom and t quantiles, as lm and nls; if
#'   FALSE, the maximum likelihood estimate and normal quantiles. See \code{\link{gcTMB}}.
#' @param sdreport Logical. If TRUE (default), the fit computes the covariance matrix of the
#'   estimates (TMB's sdreport), for standard errors, Wald intervals and the standard errors of
#'   predictions. FALSE skips it, which saves about a quarter of the time of a fit with random
#'   effects, e.g. in simulation studies that only need the estimates and the log-likelihood: the
#'   standard errors, \code{vcov} and the Wald intervals are then NA, and \code{\link{converged}}
#'   does not check the Hessian. Profile and bootstrap intervals are still available.
#' @param restart Logical. If TRUE, a fit whose optimiser does not report convergence, or whose
#'   largest absolute gradient is above 1e-3, is optimised once more from its estimates, and keeps the
#'   better of the two optima. With few groups and several correlated random effects, the likelihood
#'   is often so flat in some directions that nlminb reports "singular convergence" or "false
#'   convergence", or stops with a gradient slightly above the tolerance, at the maximum: the restart
#'   then hardly changes the log-likelihood. The improvement of the log-likelihood by the restart is
#'   kept in the fit (\code{fit$restart$improvement}), and \code{\link{converged}} accepts a fit that
#'   the restart improved by less than its argument lltol. FALSE (default) for no restart.
#' @param onerror Character string. What a fit does when the fitting itself fails with an error, e.g.
#'   the optimiser stops at a point where the likelihood cannot be evaluated, or the covariance of the
#'   estimates (sdreport) cannot be computed: "stop" (default) stops with the error, as before;
#'   "warning" turns the error into a warning and returns the gcTMB object, so that a loop over many
#'   fits, e.g. a simulation study, keeps the failed ones. If the optimiser fails, the returned
#'   object has the starting values as its estimates (or the optimum before a failed refinement);
#'   if only sdreport fails, it has the estimates without standard errors, as with
#'   \code{sdreport=FALSE}. The error message is kept in \code{fit$fiterror}, printed by print and
#'   summary, and \code{\link{converged}} reports the fit as not converged. Errors in the
#'   specification of the model (formulas, data, links, starting values) stop in either case.
#'
#' @return A list with the control parameters.
#' @examples
#' gcTMBcontrol(iter.max=2000)
#' # small-sample inference switched off: maximum likelihood estimate of the residual sd
#' gcTMB(Height ~ DAS, data=FabaHeights25, control=list(smallsample=FALSE))
#' # estimates only, without standard errors
#' gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25, control=list(sdreport=FALSE))
#' # a second optimisation from the estimates when the first one does not converge
#' fit <- gcTMB(Height ~ DAS, rform=b1 + b2 + b3 ~ 1 | Plant_ID, data=FabaHeights25,
#'              control=list(restart=TRUE))
#' converged(fit)
#' # a loop over models: a fit that fails is returned with a warning instead of stopping the loop
#' fits <- lapply(c("logistic", "gompertz"), function(m)
#'   gcTMB(Height ~ DAS, data=FabaHeights25, model=m, control=list(onerror="warning")))
#' # the error of each fit, NULL if it did not fail
#' lapply(fits, function(f) f$fiterror)
#' @export
#'

gcTMBcontrol <- function(optimiser="nlminb", iter.max=1000, eval.max=2000, smallsample=TRUE, sdreport=TRUE,
                         restart=FALSE, onerror=c("stop", "warning")){
    if (!identical(optimiser, "nlminb")) stop("the only available optimiser is \"nlminb\".")
    onerror <- match.arg(onerror)
    list(optimiser=optimiser, iter.max=iter.max, eval.max=eval.max, smallsample=smallsample, sdreport=sdreport,
         restart=restart, onerror=onerror)
}
