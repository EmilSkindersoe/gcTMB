#' Set control parameters for gcTMB
#'
#' @param optimiser Character string. Specifies the optimiser to use
#' @param iter.max Integer. Maximum number of iterations of the optimiser, see \code{\link{nlminb}}.
#' @param eval.max Integer. Maximum number of evaluations of the objective function, see \code{\link{nlminb}}.
#' @param smallsample Logical. If TRUE (default), gaussian and lognormal models without random effects use the unbiased
#'   residual variance with n - p residual degrees of freedom and t quantiles, as lm and nls; if
#'   FALSE, the maximum likelihood estimate and normal quantiles. See \code{\link{gcTMB}}.
#'
#' @return A list with the control parameters.
#' @examples
#' gcTMBcontrol(iter.max=2000)
#' # small-sample inference switched off: maximum likelihood estimate of the residual sd
#' gcTMB(height ~ DAS, data=FabaHeights25, control=list(smallsample=FALSE))
#' @export
#'

gcTMBcontrol <- function(optimiser="nlminb", iter.max=1000, eval.max=2000, smallsample=TRUE){
    list(optimiser=optimiser, iter.max=iter.max, eval.max=eval.max, smallsample=smallsample)
}
