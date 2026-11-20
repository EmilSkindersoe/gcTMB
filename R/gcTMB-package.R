#' @section Growth models and fitting:
#' \code{\link{gcTMB}} fits the growth curves: logistic (Richards), Gompertz, Hill, Weibull, normal and
#' gamma, each with an asymmetry b5 and an additive b6, and a straight line. The parameters are the asymptote b1 (or the
#' maximum growth rate, \code{b1type="rate"}), the time b2 of the maximum growth rate or of half the
#' asymptote, the scale, slope or shape b3, the threshold b4 of a truncated curve, the asymmetry
#' b5 and the additive b6 (the response from which the curve grows), each with a link and optionally fixed. \code{\link{gcTMBstart}} gives the starting values and
#' \code{\link{gcTMBcontrol}} the settings of the optimisation and the inference.
#'
#' @section Effects, distributions and autocorrelation:
#' Fixed effects of covariates on the parameters (\code{fform}) and random effects (\code{rform}),
#' correlated across parameters, with several terms, nested or crossed factors and random slopes; the
#' families gaussian, lognormal, poisson, binomial, betabinomial, beta, nbinom1 and nbinom2; and
#' autocorrelated residuals within individuals, \code{\link{corAR1}} and \code{\link{corCAR1}}.
#'
#' @section Inference:
#' \itemize{
#'   \item \code{\link{summary.gcTMB}}, \code{coef}, \code{vcov}, \code{logLik}, \code{\link{converged}};
#'   \item confidence intervals: \code{\link{confint.gcTMB}} (profile likelihood, Wald, bootstrap),
#'     \code{\link{profile.gcTMB}}, \code{\link{bootstrap}};
#'   \item tests and model choice: \code{\link{anova.gcTMB}} (likelihood ratio, F and parametric
#'     bootstrap tests), \code{\link{criteria}} (information and goodness-of-fit criteria),
#'     \code{\link{pairwise}} (pairwise comparisons with single-step adjustment).
#' }
#'
#' @section Timing and growth rates:
#' \code{\link{timeto}}, the times to reach a percentage of the asymptote or a given response, and
#' \code{\link{maxrate}}, the maximum growth rate, its timing and the fraction of the asymptote reached
#' then, with standard errors and \code{confint} and \code{profile} methods.
#'
#' @section Predictions, diagnostics and simulation:
#' \code{\link{predict.gcTMB}} (curves and parameters, with confidence and prediction intervals),
#' \code{fitted}, \code{\link{residuals.gcTMB}} (Pearson, randomised quantile, normalized),
#' \code{\link{plot.gcTMB}} (validation plots and curves), \code{\link{simulate.gcTMB}} (at the observed
#' or new times and individuals), and \code{\link{update.gcTMB}} (refit, or set the parameters of a true
#' model for simulation).
#'
#' @section Extractors:
#' \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{VarCorr}}, \code{sigma}, \code{\link{corpars}}
#' and \code{\link{varpars}}, compatible with nlme and lme4 (see \code{\link{gcTMB}}, Relation to nlme).
#'
#' @section Data and vignettes:
#' The data sets \code{\link{FabaHeights25}} and \code{\link{FabaCounts23}}. The vignettes:
#' \code{vignette("gcTMB")}, an overview of most features; \code{vignette("growth-models")},
#' \code{vignette("fixed-random-effects")}, \code{vignette("timing-rates")} and
#' \code{vignette("autocorrelation-simulation")}.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib gcTMB, .registration = TRUE
#' @useDynLib gcTMB_TMBExports
#' @importFrom stats AIC BIC approx optim coef coefficients confint cor fitted lm lm.fit lm.wfit logLik
#' @importFrom stats model.frame model.matrix model.response nlminb nobs pbeta pbinom
#' @importFrom stats plogis pnorm ppois predict printCoefmat profile qchisq qlogis qnorm rlnorm
#' @importFrom stats quantile residuals runif sd setNames simulate splinefun vcov
#' @importFrom stats complete.cases terms .getXlevels na.pass as.formula
#' @importFrom stats rnorm rpois rbinom rbeta rnbinom pnbinom anova pchisq pf lowess ppoints qt qf df.residual uniroot optimize dlogis
#' @importFrom stats cov2cor dgamma dnorm na.omit p.adjust.methods pgamma pt
#' @importFrom graphics plot
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix t bdiag
#' @importFrom reformulas subbars mkReTrms findbars
#' @importFrom TMB MakeADFun sdreport
#' @importFrom methods as
## usethis namespace: end
NULL

# release the DLLs when the namespace is unloaded, e.g. by unloadNamespace("gcTMB"): on Windows a
# loaded DLL cannot be overwritten, so this lets a new version be installed from the same R session
.onUnload <- function(libpath){
  for (dll in c("gcTMB_TMBExports", "gcTMB")){
    try(library.dynam.unload(dll, libpath), silent=TRUE)
  }
}
