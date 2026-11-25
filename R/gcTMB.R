#' Fit a growth curve
#'
#' Fits a nonlinear growth curve over time by maximum likelihood with the Template Model Builder
#' (TMB), optionally with covariates on the parameters and random effects, which are integrated out by
#' the Laplace approximation.
#'
#' @param form Formula, with the names of the response and the explanatory variable x, e.g. \code{weight ~ day}. x is usually time, but can be any variable along which the response grows towards an asymptote, e.g. the diameter of trees in a height-diameter model \code{height ~ dbh}.
#' @param fform Formula, list of formulas, or NULL. Fixed effects for the model parameters: a formula of the form b1 + b2 ~ Treatment, or a list of such formulas, e.g. list(b1 ~ Treatment, b2 ~ Treatment + Site). Parameters not listed have a single intercept (~ 1); NULL (default) gives every parameter a single intercept, one curve for all the data. A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~ 1), is also accepted.
#' @param rform Formula, list of formulas, or NULL. Random effects: a formula of the form b1 + b2 ~ 1 | ID for random intercepts, or b1 ~ x | ID for random intercepts and slopes in a covariate x, or a list of such formulas, e.g. list(b1 + b2 ~ 1 | ID, b1 ~ 1 | site). The random effects within a formula are correlated across the listed parameters and their intercepts and slopes; different formulas and the levels of nested factors (1 | site/ID) are independent terms. With a double bar, all the random coefficients of the formula are independent, across the parameters as well as their intercepts and slopes: b1 + b2 ~ 1 || ID is list(b1 ~ 1 | ID, b2 ~ 1 | ID), a diagonal covariance matrix as random=pdDiag(b1 + b2 ~ 1) in nlme. A list of formulas combines the two into a block-diagonal covariance matrix: list(b1 + b2 ~ 1 | ID, b3 ~ 1 | ID) gives correlated random effects of b1 and b2, and a random effect of b3 independent of both, as random=pdBlocked(list(pdSymm(b1 + b2 ~ 1), pdSymm(b3 ~ 1))) in nlme. The random slopes of a parameter usually go with a fixed effect of the covariate in fform. NULL (default) for no random effects.
#' @param data Data frame. Observations with missing values in the variables of the model or in the weights are left out.
#' @param start Numeric vector or NULL. Starting values for each model coefficient. NA entries are replaced by the starting values of \code{\link{gcTMBstart}}, and NULL (default) uses them for all coefficients. Alternatively, a vector named by parameter, e.g. c(b1=10), gives the same starting value to all coefficients of each named parameter.
#' @param fix Logical vector or NULL. One entry for each model coefficient. If TRUE, a model coefficient is fixed to its value in start. Alternatively, a vector named by parameter, e.g. c(b1=TRUE), fixes all coefficients of the named parameters; their values must then be given in start, e.g. start=c(b1=10). NULL (default) fixes only the parameters that are fixed by default. The threshold (start of growth) b4 and the asymmetry b5 are fixed by default, at b4 = 0 (the curve starts at zero at t = 0) and b5 = 1 (the symmetric curve) unless given in start; so is the intercept b4 of the chapmanrichards and korf models, at the breast height of 1.3, while the rates b4 and b5 of the preecebaines model are estimated; fix=c(b4=FALSE) or fix=c(b5=FALSE), or an unnamed fix vector, estimates them. Without truncation, b4 is not used and always fixed, at 0 unless given in start.
#' @param correlation A correlation structure or NULL. The autocorrelation of the residuals within
#'   groups over time, \code{\link{corAR1}}, \code{\link{corCAR1}}, \code{\link{corExp}} or
#'   \code{\link{corGaus}} (or the same structures of
#'   the package nlme), e.g. \code{corCAR1(form=~ time | ID)}, for the gaussian and lognormal (on the
#'   log scale) families. NULL (default) for independent residuals. See \code{\link{corStruct}}; with
#'   a correlation structure, the weights are precision weights, with residual variance sigma^2/w.
#' @param variance A variance function or NULL. Heteroscedastic residuals with the standard deviation
#'   sigma g, \code{\link{varIdent}}, \code{\link{varPower}} or \code{\link{varConstPower}} (or the
#'   same functions of the package nlme), e.g. \code{varIdent(form=~ 1 | Treatment)} for a variance
#'   for each treatment, or \code{varConstPower(form=~ rate(.))} for a standard deviation that grows
#'   with the growth rate of the curve, for the gaussian and lognormal (on the log scale) families.
#'   NULL (default) for a constant variance. See \code{\link{varFunc}}.
#' @param lower Numeric, a single value or a vector. Lower bounds for the model coefficients b1-b5, a single value or a vector with the same length as start. Bounds of fixed coefficients are ignored; variance parameters are unbounded.
#' @param upper Numeric, a single value or a vector. Upper bounds for the model coefficients b1-b5, a single value or a vector with the same length as start.
#' @param family Character string. The distribution of the response, one of "gaussian", "lognormal", "binomial", "betabinomial", "beta", "poisson", "nbinom2", and "nbinom1". For the lognormal family, log(y) is normal with mean log f(t) and standard deviation sigma: the curve f is the median of a positive response, with a constant coefficient of variation; the mean is f(t) exp(sigma^2/2). The negative binomial families are for overdispersed counts, with variance mu (1 + mu/theta) for nbinom2 and mu (1 + phi) for nbinom1, and the beta-binomial family for overdispersed binomial counts, with variance n mu (1 - mu) (1 + (n - 1)/(phi + 1)) for n trials, as in glmmTMB.
#' @param model Character string. The growth model: "logistic" (default), "gompertz", "hill", "weibull", "normal" or "gamma", each with an asymmetry parameter b5, "beta" or "beta_capped" for determinate growth, "chapmanrichards" or "korf" for height-diameter models in forestry, "preecebaines" for human growth in height, or "linear" for a straight line, see Details.
#' @param truncated Logical. If TRUE (default), the curve is truncated at a threshold time b4, by default fixed at 0: the response is zero up to b4 and then grows, see Details. FALSE for the curves without truncation. The beta, beta_capped, chapmanrichards, korf and preecebaines curves are never truncated.
#' @param link Named character vector or NULL. Links of the parameters b1-b5: NULL for the defaults, or links named by parameter, e.g. \code{link=c(b1="log", b3="inverse")}. The default is the identity link for all parameters, except the shape b3 of the weibull, gamma and korf models, which must be larger than 1 and has the link "logm1", log(b3 - 1), by default (\code{link=c(b3="identity")} gives the shape itself). The links available for each parameter depend on its role in the model:
#'   \itemize{
#'     \item the asymptote b1: "identity", "log" or "logit" (for proportions); the slope b1 of the linear model: "identity" or "log";
#'     \item the time b2 (inflection point, half the asymptote, middle of puberty), the rate b2 of chapmanrichards and the coefficient b2 of korf: "identity" or "log"; the intercept b2 of the linear model: "identity";
#'     \item a scale b3 (logistic, gompertz, normal): "identity", "inverse", with which its coefficients are the growth rate 1/b3, or "log"; the shape b3 of weibull, gamma and korf: "logm1", "identity" or "log"; any other b3 (the slope of hill, the shape of chapmanrichards, the end of growth of beta, the height of preecebaines): "identity" or "log";
#'     \item b4 (the threshold or start of growth, the intercept of chapmanrichards and korf, the rate of preecebaines) and the asymmetry b5: "identity" or "log".
#'   }
#'   A link that a model does not allow gives an error naming the allowed ones. The coefficients, and start, fix, lower, upper and curves of \code{\link{timeto}}, are on the scale of the links, see Details.
#' @param b1type Character string. The meaning of b1: "asymptote" (default), or "rate" for the maximum growth rate of the curve, see Details. With "rate", covariates and random effects on b1 are effects on the maximum growth rate, which \code{\link{anova.gcTMB}} and \code{\link{partest}} can test.
#' @param weights Numeric vector or NULL. Weights, one for each row of data; NULL (default) for a weight of 1 for every row. Unlike the
#'   \code{weights} of nlme, which are variance functions (\code{varFunc}), these are numeric case
#'   weights of the log-likelihood; with a correlation structure, they are precision weights, as
#'   \code{weights=varFixed(~ 1/w)} in nlme.
#' @param gsl Numeric vector or NULL. Lower bounds for a grid search, same length as start; NA for the default range. If gsl, gsu and gsn are all NULL (default), no grid search is performed and the starting values come from linear models (method "lm" of gcTMBstart); if any of them is given, NULL for the others means the default ranges and 10 grid points. See \code{\link{gcTMBstart}}.
#' @param gsu Numeric vector or NULL. Upper bounds for a grid search, same length as start; NA for the default range, see gsl for NULL.
#' @param gsn Integer, a single number or a vector, or NULL. Number of grid points for each parameter, a single number or a vector with the same length as start. \code{gsn=10} alone runs a grid search over the default ranges; see gsl for NULL.
#' @param control List. Control parameters, the arguments of \code{\link{gcTMBcontrol}}, e.g. \code{list(iter.max=2000)}. See \code{\link{gcTMBcontrol}} for details.
#'
#' @details The logistic, gompertz, hill, weibull, normal and gamma curves are sigmoid growth curves,
#' often used in agronomy and ecology; they increase from 0 towards the asymptote b1:
#' \itemize{
#'   \item logistic: f(t) = b1 / (1 + exp(-(t - b2)/b3)), symmetric around b2, where the curve
#'     reaches half of its asymptote;
#'   \item gompertz: f(t) = b1 exp(-exp(-(t - b2)/b3)); it reaches exp(-1), about 37\%, of its
#'     asymptote at b2 and approaches the asymptote more slowly than it leaves zero;
#'   \item hill: f(t) = b1 / (1 + (b2/t)^b3), a logistic curve in log t (log-logistic),
#'     b1 / (1 + exp(-b3 (log t - log b2))); b2 is the time at which the curve reaches half of its
#'     asymptote, and b3 the slope on the log time axis;
#'   \item weibull: f(t) = b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3)); b2 is the inflection point and b3 > 1
#'     the shape;
#'   \item normal: f(t) = b1 Phi((t - b2)/b3), with the standard normal distribution function Phi;
#'     symmetric around b2, where the curve reaches half of its asymptote;
#'   \item gamma: f(t) = b1 F(t), with F the cumulative distribution function (CDF) of the gamma
#'     distribution with shape b3 > 1 and scale b2/(b3 - 1), whose density has its mode at b2: b2 is the inflection
#'     point and b3 the shape;
#'   \item linear: f(t) = b1 t + b2, a straight line with slope b1 and intercept b2, for comparison
#'     with the growth curves. It has no asymptote, does not use b3-b5 and is not truncated;
#'     \code{\link{timeto}} gives the times at which it reaches given responses
#'     (\code{type="absolute"}), and \code{\link{maxrate}} its constant rate b1.
#' }
#' Five more curves are never truncated: the beta curves of determinate growth, the height-diameter
#' models chapmanrichards and korf, common in forestry, and the preecebaines model of human growth:
#' \itemize{
#'   \item chapmanrichards: f(t) = b4 + b1 (1 - exp(-b2 t))^b3, the Chapman-Richards curve of
#'     Richards (1959) in its usual form for height-diameter models,
#'     zero above the intercept b4 up to t = 0: b1 is the asymptote, b2 the rate, b3 the shape and b4
#'     the intercept, fixed by default at the breast height of 1.3 m (\code{start=c(b4=0)} for a curve
#'     from zero, \code{fix=c(b4=FALSE)} to estimate it).
#'     For b3 > 1 the inflection point is at log(b3)/b2, where the curve reaches (1 - 1/b3)^b3 of its
#'     asymptote; for b3 <= 1 it grows fastest at t = 0;
#'   \item beta: the beta growth function of Yin et al. (2003), eqn 8 (with b4 = 0 there),
#'     f(t) = b1 (1 + (b3 - t)/(b3 - b2)) ((t - b4)/(b3 - b4))^((b3 - b4)/(b3 - b2)) for
#'     b4 < t < 2 b3 - b2, and zero otherwise: b1 is the maximum w_max, b2 the time t_m of the maximum
#'     growth rate, b3 the end of growth t_e, where the curve reaches b1, and b4 the start of growth
#'     t_b, b4 < b2 < b3. After b3 the curve falls back to zero at 2 b3 - b2, for responses that
#'     decline after their maximum, e.g. biomass lost after maturity;
#'   \item beta_capped: the same curve up to b3, and b1 from b3 on, the beta sigmoid growth function
#'     of Yin et al. for determinate growth that stops at its maximum;
#'   \item korf: f(t) = b4 + b1 exp(b2/((1 - b3) t^(b3 - 1))) for t > 0, and b4 up to t = 0, the Korf
#'     curve: b1 is the asymptote above the intercept b4, b2 > 0 a coefficient, b3 > 1 the shape, and
#'     b4 the intercept as for chapmanrichards. The inflection point is at t = (b2/b3)^(1/(b3 - 1)),
#'     where the curve reaches exp(-b3/(b3 - 1)) of its asymptote;
#'   \item preecebaines: model 1 of Preece and Baines (1978) for human growth in height,
#'     f(t) = b1 - 2 (b1 - b3)/(exp(b4 (t - b2)) + exp(b5 (t - b2))): the adult height b1, the age b2 in
#'     the middle of puberty, the height b3 at b2, and the rates b4 < b5 of growth before and during
#'     puberty, which are estimated by default. \code{\link{maxrate}} gives the peak of the pubertal
#'     growth spurt (the peak height velocity); \code{b1type="rate"} is not available.
#' }
#' \strong{Relation to nlme.} Parts of the interface follow the package nlme (Pinheiro and Bates,
#' 2000), so that models can be moved between \code{nlme::nlme}/\code{nlme::gnls} and gcTMB:
#' \itemize{
#'   \item the fixed and random effects are formulas for the parameters, as the arguments \code{fixed}
#'     and \code{random} of \code{nlme::nlme}: \code{fform=b1 + b2 ~ trt} or
#'     \code{fform=list(b1 ~ trt, b2 ~ 1)}, and \code{rform=b1 + b2 ~ 1 | ID}. The random effects of a
#'     formula have an unstructured covariance matrix, as \code{pdSymm}; uncorrelated random effects,
#'     as \code{pdDiag}, are separate terms, e.g. \code{rform=list(b1 ~ 1 | ID, b2 ~ 1 | ID)} or
#'     \code{b1 ~ x || ID}. Nested factors (\code{1 | site/ID}) and the double bar use the syntax
#'     of lme4;
#'   \item the correlation structures \code{\link{corAR1}}, \code{\link{corCAR1}}, \code{\link{corExp}}
#'     and \code{\link{corGaus}} have the
#'     arguments and the parameterisation of nlme, whose objects are accepted as well, and
#'     \code{residuals(fit, type="normalized")} gives the normalized residuals of nlme;
#'   \item the variance functions \code{\link{varIdent}}, \code{\link{varPower}} and
#'     \code{\link{varConstPower}} go in the argument \code{variance} (the argument \code{weights} of
#'     nlme), with the arguments and parameterisation of nlme, plus the growth rate \code{rate(.)} as a
#'     covariate;
#'   \item \code{\link{fixef}}, \code{\link{ranef}} and \code{\link{VarCorr}} are generics with the names
#'     of those of nlme, and work for gcTMB, nlme and lme4 fits whichever package is attached last;
#'   \item \code{predict} and \code{fitted} with \code{random=FALSE} correspond to \code{level=0} of
#'     nlme (the population curves), and with \code{random=TRUE} to the innermost level.
#' }
#' The differences: gcTMB fits by maximum likelihood (the Laplace approximation with random
#' effects), so its log-likelihoods, AIC and likelihood ratio tests compare with nlme fits with
#' \code{method="ML"}, not with the default REML of \code{nlme::lme} and \code{nlme::gls}; the
#' integration over the random effects is exact for gaussian responses with random effects that enter
#' the curve linearly, e.g. random asymptotes, where nlme's linearisation can stop slightly short of
#' the maximum. \code{coef}
#' returns the estimated coefficients and variance parameters, not the coefficients of each group as
#' \code{coef} of nlme (use \code{predict(fit, type="parameters")} for those); \code{weights} are
#' numeric weights, not variance functions; and the results of \code{ranef} and \code{VarCorr} have
#' the format of lme4 (a list of data frames, and of covariance matrices). \code{\link{confint.gcTMB}}
#' gives the intervals of \code{intervals} of nlme.
#'
#' The asymmetry b5, the truncation at the threshold b4, the maximum growth rates and inflection
#' points, the rate 1/b3 and the small-sample inference are described in the vignette,
#' \code{vignette("gcTMB")}; the links in the argument \code{link}.
#'
#' With \code{b1type="rate"}, b1 is the maximum growth rate of the curve, the slope at its inflection
#' point (see \code{\link{maxrate}}), instead of the asymptote; the asymptote is b1 divided by the
#' maximum growth rate of the curve with asymptote 1, which depends on b2-b5. The model is the same as
#' with the asymptote, unless b1 has random effects, which are then random maximum growth rates.
#' \code{predict(type="parameters")} gives the asymptote as well. The maximum growth rate has a closed
#' form for all the models with b5 = 1, and for the logistic (Richards) and hill curves also with an
#' asymmetry b5; the gompertz, weibull, normal and gamma models need b5 fixed at 1. The hill curve has
#' a maximum growth rate at a positive time only for b3 b5 > 1. For the linear model, b1 is the rate
#' anyway.
#'
#' @return A gcTMB object, with methods for \code{print}, \code{summary}, \code{coef}, \code{vcov},
#'   \code{predict}, \code{confint}, \code{anova}, \code{plot} and others, see \code{\link{gcTMB-methods}}.
#' @references Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS.
#'   Springer, New York.
#'
#'   Preece, M.A. and Baines, M.J. (1978). A new family of mathematical models describing the human
#'   growth curve. Annals of Human Biology 5, 1-24.
#'
#'   Richards, F.J. (1959). A flexible growth function for empirical use. Journal of Experimental
#'   Botany 10, 290-300.
#'
#'   Yin, X., Goudriaan, J., Lantinga, E.A., Vos, J. and Spiertz, H.J. (2003). A flexible sigmoid
#'   function of determinate growth. Annals of Botany 91, 361-371.
#' @seealso \code{\link{maxrate}} for the maximum growth rate and its timing,
#'   \code{\link{timeto}} for the times to reach a level of the asymptote,
#'   \code{\link{gcTMBstart}} for the starting values.
#' @examples
#' # logistic growth of the faba bean plants
#' fit <- gcTMB(Height ~ DAS, data=FabaHeights25, model="logistic")
#' summary(fit)
#'
#' # the other models, compared by AIC
#' fitg <- gcTMB(Height ~ DAS, data=FabaHeights25, model="gompertz")
#' fith <- gcTMB(Height ~ DAS, data=FabaHeights25, model="hill")
#' fitw <- gcTMB(Height ~ DAS, data=FabaHeights25, model="weibull")
#' fitn <- gcTMB(Height ~ DAS, data=FabaHeights25, model="normal")
#' fitgam <- gcTMB(Height ~ DAS, data=FabaHeights25, model="gamma")
#' fitlin <- gcTMB(Height ~ DAS, data=FabaHeights25, model="linear")
#' # the forestry models from zero instead of their default intercept, the breast height of 1.3 m
#' fitcr <- gcTMB(Height ~ DAS, data=FabaHeights25, model="chapmanrichards", start=c(b4=0))
#' fitb <- gcTMB(Height ~ DAS, data=FabaHeights25, model="beta_capped")
#' fitk <- gcTMB(Height ~ DAS, data=FabaHeights25, model="korf", start=c(b4=0))
#' AIC(fit, fitg, fith, fitw, fitn, fitgam, fitlin, fitcr, fitb, fitk)
#'
#' # a height-diameter model from the breast height of 1.3 m, the default intercept b4
#' gcTMB(Height ~ DBH, data=PineHeightDiameter, model="chapmanrichards")
#'
#' # the Richards curve: the logistic model with an estimated asymmetry b5
#' fitr5 <- gcTMB(Height ~ DAS, data=FabaHeights25, fix=c(b5=FALSE))
#' coef(fitr5)
#'
#' # the rate 1/b3 instead of the scale b3
#' gcTMB(Height ~ DAS, data=FabaHeights25, link=c(b3="inverse"))
#'
#' # treatment effects on the asymptote as ratios, with the log link
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitl <- gcTMB(Height ~ DAS, fform=b1 ~ trt, data=FabaHeights25, link=c(b1="log"))
#' exp(coef(fitl)[2:5])
#'
#' # treatment effects on the asymptote and random effects of the plants
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitre <- gcTMB(Height ~ DAS, fform=b1 ~ trt, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
#' summary(fitre)
#'
#' # treatment effects on the maximum growth rate, b1 with b1type = "rate", tested by a likelihood
#' # ratio test and compared with the first treatment
#' fitr <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25, b1type="rate")
#' fitr0 <- gcTMB(Height ~ DAS, fform=b2 ~ trt, data=FabaHeights25, b1type="rate")
#' anova(fitr0, fitr)
#' partest(fitr, newdata=data.frame(trt=levels(FabaHeights25$trt)))
#'
#' # residuals of each plant correlated over time, continuous-time AR(1)
#' fitc <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25,
#'               correlation=corCAR1(form=~ DAS | Plant_ID))
#' corpars(fitc)
#'
#' # pod counts with a poisson distribution
#' gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson")
#'
#' # no pods before a threshold b4, after which the count grows
#' fitt <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE))
#' coef(fitt)
#' @export

gcTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, correlation=NULL, variance=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", truncated=TRUE, link=NULL, b1type=c("asymptote", "rate"), weights=NULL, gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  b1type <- match.arg(b1type)
  # control list
  control <- do.call(gcTMBcontrol, control)
  dd <- gcTMBdata(form, fform=fform, rform=rform, data=data, weights=weights, family=family, model=model, truncated=truncated, link=link, b1type=b1type,
                  correlation=correlation, variance=variance)
  fixnamed <- !is.null(names(fix))
  start <- defaultstart(parvector(start, dd$npi, NA, "start"), dd)
  fix <- defaultfix(fix, dd$npi, dd$mod)
  if (fixnamed && any(fix & is.na(start) & dd$npi < 4)){
    stop("give the values of the fixed parameters in start, e.g. start=c(b1=10), fix=c(b1=TRUE).")
  }
  fix[dd$npi %in% unusedpars(dd$mod, dd$trunc)] <- TRUE
  checkb1rate(dd, start, fix)

  # starting values for parameters with NA in start
  method <- if (is.null(gsl) & is.null(gsu) & is.null(gsn)) "lm" else "grid"
  start <- gcTMBstartvals(dd, start, fix, method=method, gsl=gsl, gsu=gsu, gsn=gsn)

  obj <- gcTMBobj(dd, start, fix, rform=rform)
  # random effects: if the Laplace approximation fails at the starting values, e.g. starting values
  # far from the population curve, start from the fit without random effects
  if (!is.null(rform) && !is.finite(suppressWarnings(obj$fn(obj$par)))){
    fit0 <- tryCatch(suppressWarnings(gcTMB(form, fform=fform, rform=NULL, data=data, start=start, fix=fix, weights=weights,
                                            lower=lower, upper=upper, family=family, model=model, truncated=truncated,
                                            link=link, b1type=b1type, control=control)),
                     error=function(e) NULL)
    if (!is.null(fit0) && is.finite(fit0$opt$objective)){
      start[!fix] <- fit0$estimates[seq_len(sum(!fix))]
      obj <- gcTMBobj(dd, start, fix, rform=rform)
    }
  }

  # bounds apply to the free coefficients only, not to the variance parameters
  bounds <- function(b, inf){
    if (length(b) == 1) b <- rep(b, length(start))
    if (length(b) != length(start)) stop("lower and upper must be a single value or have the same length as start.")
    c(b[!fix], rep(inf, length(obj$par) - sum(!fix)))
  }
  lower <- bounds(lower, -Inf)
  upper <- bounds(upper, Inf)

  ctrl <- list(iter.max=control$iter.max, eval.max=control$eval.max)
  # control onerror="warning": an error of the fitting is kept in fiterror, and the fit continues
  # from the fallback, e.g. the starting values
  fiterror <- NULL
  catchfit <- function(expr, fallback){
    if (control$onerror == "stop") return(expr)
    tryCatch(expr, error=function(e){
      fiterror <<- c(fiterror, conditionMessage(e))
      fallback()
    })
  }
  failedopt <- function(par) function(){
    objective <- tryCatch(as.numeric(obj$fn(par)), error=function(e) NaN)
    list(par=par, objective=objective, convergence=1L, message=paste("error:", fiterror[length(fiterror)]),
         iterations=0L, evaluations=c("function"=0L, gradient=0L))
  }
  opt <- catchfit({
    if (is.null(rform)){
      # Newton steps with the Hessian; quasi-Newton steps if the Hessian is not finite on the way, e.g.
      # for parameters far outside the range of the data
      tryCatch(nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper, control=ctrl),
               error=function(e) nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl))
    } else {
      opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl)
      retryvariances(obj, opt, dd$re$terms, lower, upper, ctrl)
    }
  }, failedopt(obj$par))
  # truncated curves: the likelihood has a kink in the threshold b4 at each observation time, where
  # the optimiser often stops before the other parameters are optimal
  i4 <- which(dd$npi[!fix] == 4)
  if (isTRUE(dd$trunc) && length(i4) > 0 && is.finite(opt$objective) && opt$convergence != 0){
    opt0 <- opt
    opt <- catchfit(thresholdopt(obj, opt, i4, lower, upper, ctrl, dd, start, fix), function() opt0)
  }
  # restart: once more from the estimates, if the optimiser did not report convergence or the gradient
  # is large; the better of the two optima is kept, and the improvement shows whether the first one
  # was already at the maximum (see converged)
  restart <- NULL
  if (isTRUE(control$restart) && is.finite(opt$objective)){
    gr0 <- tryCatch(max(abs(obj$gr(opt$par))), error=function(e) Inf)
    if (opt$convergence != 0 || !is.finite(gr0) || gr0 > 1e-3){
      opt2 <- tryCatch(nlminb(start=opt$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl),
                       error=function(e) NULL)
      improvement <- if (!is.null(opt2) && is.finite(opt2$objective)) opt$objective - opt2$objective else NA_real_
      restart <- list(improvement=improvement, convergence=opt$convergence, message=opt$message, maxgrad=gr0)
      if (!is.na(improvement) && improvement >= 0) opt <- opt2
      obj$fn(opt$par)
    }
  }
  if (!is.null(fiterror)){
    warning("the fit failed: ", fiterror[1], if (is.null(opt$iterations) || opt$iterations == 0)
      "; the object is returned at the starting values" else "; the object is returned at the last optimum", call.=FALSE)
  } else if (!is.finite(opt$objective)){
    warning("the negative log-likelihood is not finite at the estimates; some observations may be impossible under the model, e.g. positive counts where the curve is zero.", call.=FALSE)
  } else if (opt$convergence != 0){
    warning("the optimiser did not converge: ", opt$message, call.=FALSE)
  }
  res <- list()
  res$call <- match.call()
  res$data <- dd$data
  res$formula <- form
  res$family <- family
  res$model <- model
  res$truncated <- dd$trunc
  res$b1type <- b1type
  res$link <- dd$link
  res$start <- start
  res$fix <- fix
  res$Xs <- dd$Xs
  res$fform <- dd$fform
  res$fterms <- dd$fterms
  res$coefnames <- coefnames(dd$Xs)
  res$rform <- rform
  res$reterms <- dd$re$terms
  res$correlation <- correlation
  res$variance <- variance
  res$varlevels <- dd$var$levels
  # the grouping factor of the random effects that groups the residuals, for a form without a group
  if (!is.null(correlation) && is.null(corform(correlation)$group) && length(dd$re$terms) > 0){
    res$corgroup <- dd$re$terms[[which.max(vapply(dd$re$terms, function(tm) length(tm$levels), integer(1)))]]$group
  }
  res$estimates <- opt$par
  res$lower <- lower
  res$upper <- upper
  res$opt <- opt
  res$restart <- restart
  res$control <- control
  res$obj <- obj
  # without random effects, the Hessian by automatic differentiation: the default finite differences of
  # sdreport use absolute steps, which fail for coefficients on a small scale
  # a failed optimisation has no covariance of the estimates; a failed sdreport gives the estimates
  # without standard errors
  res$sdrl <- if (isFALSE(control$sdreport) || !is.null(fiterror)) nosdreport(obj, opt, dd$re$terms) else
    catchfit(sdreport(obj, par.fixed=opt$par, hessian.fixed=if (is.null(rform)) obj$he(opt$par)),
             function(){
               warning("the fit failed: ", fiterror[length(fiterror)], "; the object is returned without standard errors", call.=FALSE)
               nosdreport(obj, opt, dd$re$terms)
             })
  res$fiterror <- fiterror
  # NaN standard errors of a Hessian that is not positive definite are reported by summary.gcTMB
  res$ssdrl <- if (inherits(res$sdrl, "nosdreport")) nosdrsummary(res$sdrl) else suppressWarnings(summary(res$sdrl))
  # residual degrees of freedom for small-sample inference of gaussian and lognormal models without
  # random effects (lognormal: gaussian for log(y)): observations with positive weights minus free coefficients
  if (isTRUE(control$smallsample) && family %in% c("gaussian", "lognormal") && length(dd$re$terms) == 0 && is.null(correlation) && is.null(variance)){
    rdf <- sum(dd$w > 0) - sum(!fix)
    if (rdf > 0) res$rdf <- rdf else warning("no residual degrees of freedom; the standard errors use the maximum likelihood estimate of the residual variance.", call.=FALSE)
  }
  class(res) <- "gcTMB"
  return(res)
}


# optimum of a truncated curve with the threshold b4 at a kink of the likelihood: alternately the other
# parameters with the coefficients i4 of b4 held, which is smooth, and all parameters, until the
# negative log-likelihood no longer decreases. With b4 at an observation time of its curve and the
# other parameters converged, the kink is the optimum in b4 and the fit is reported as converged.
thresholdopt <- function(obj, opt, i4, lower, upper, ctrl, dd, start, fix){
  best <- opt
  conditional <- FALSE
  for (k in 1:5){
    p <- best$par
    fn <- function(q){ p[-i4] <- q; obj$fn(p) }
    gr <- function(q){ p[-i4] <- q; as.vector(obj$gr(p))[-i4] }
    c1 <- tryCatch(nlminb(p[-i4], fn, gr, lower=lower[-i4], upper=upper[-i4], control=ctrl), error=function(e) NULL)
    improved <- FALSE
    if (!is.null(c1) && is.finite(c1$objective) && c1$objective < best$objective - 1e-10){
      p[-i4] <- c1$par
      best$par <- p
      best$objective <- c1$objective
      best$convergence <- c1$convergence
      best$message <- c1$message
      conditional <- TRUE
      improved <- TRUE
    }
    c2 <- tryCatch(nlminb(best$par, obj$fn, obj$gr, lower=lower, upper=upper, control=ctrl), error=function(e) NULL)
    if (!is.null(c2) && is.finite(c2$objective) && c2$objective < best$objective - 1e-10){
      best <- c2
      conditional <- FALSE
      improved <- TRUE
    }
    if (!improved || best$convergence == 0) break
  }
  obj$fn(best$par)
  # the threshold of each observation's curve at one of the observation times
  pv <- start
  pv[!fix] <- best$par[seq_len(sum(!fix))]
  b4 <- linkfun(dd$lnk[4], as.vector(dd$Xs[[4]] %*% pv[dd$npi == 4]), inverse=TRUE)
  corner <- any(abs(b4 - dd$x) <= 1e-6 * pmax(1, abs(dd$x)))
  if (conditional && best$convergence == 0 && corner){
    best$message <- "converged, with the threshold b4 at an observation time (a kink of the likelihood)"
  } else if (conditional){
    best$convergence <- 1L
  }
  best
}


# the parts of sdreport used by the methods, without the covariance matrix of the estimates
# (control sdreport=FALSE): the estimates, the random effects at their conditional modes, and the
# correlation matrices of the random effects (as the ADREPORT corr), with NA standard errors;
# nosdrsummary gives the same table as summary of an sdreport
nosdreport <- function(obj, opt, terms){
  env <- obj$env
  r <- env$random
  # the conditional modes of the random effects at the estimates (not last.par.best, the best point
  # evaluated so far, which need not be the estimates)
  par.random <- if (!is.null(r)){
    # a failed fit (control onerror="warning") may not be evaluable at its estimates
    try(obj$fn(opt$par), silent=TRUE)
    env$last.par[r]
  } else numeric(0)
  pl <- env$parList(opt$par)
  toff <- 0
  corr <- unlist(lapply(terms, function(tm){
    np <- length(tm$pars)
    cr <- corrmatrix(pl$theta[toff + np + seq_len(np*(np - 1)/2)], np)
    toff <<- toff + np + np*(np - 1)/2
    as.vector(cr)
  }))
  np <- length(opt$par)
  structure(list(par.fixed=opt$par, par.random=par.random, cov.fixed=matrix(NA_real_, np, np),
                 pdHess=NA, corr=corr), class="nosdreport")
}

nosdrsummary <- function(object){
  est <- c(object$corr, object$par.fixed, object$par.random)
  nm <- c(rep("corr", length(object$corr)), names(object$par.fixed), rep("u", length(object$par.random)))
  matrix(c(est, rep(NA_real_, length(est))), ncol=2, dimnames=list(nm, c("Estimate", "Std. Error")))
}
