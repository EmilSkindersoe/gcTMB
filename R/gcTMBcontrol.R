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
#'
#' @return A list with the control parameters.
#' @examples
#' gcTMBcontrol(iter.max=2000)
#' # small-sample inference switched off: maximum likelihood estimate of the residual sd
#' gcTMB(height ~ DAS, data=FabaHeights25, control=list(smallsample=FALSE))
#' # estimates only, without standard errors
#' gcTMB(height ~ DAS, rform=b1 ~ 1 | plant_ID, data=FabaHeights25, control=list(sdreport=FALSE))
#' @export
#'

gcTMBcontrol <- function(optimiser="nlminb", iter.max=1000, eval.max=2000, smallsample=TRUE, sdreport=TRUE){
    if (!identical(optimiser, "nlminb")) stop("the only available optimiser is \"nlminb\".")
    list(optimiser=optimiser, iter.max=iter.max, eval.max=eval.max, smallsample=smallsample, sdreport=sdreport)
}
