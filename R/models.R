# Properties of the growth models. The position of a model in the list is its index (mod) in the TMB
# templates (gcfunctions.h), which define the same functions; a new model needs an entry here, at the
# same index as in the templates.
#
# All models but the linear one are of the form f(t) = b1 F(z(t)): a distribution function F,
# increasing from 0 to 1, of a base distribution (basedists), and a transformed time z(t), increasing
# in t:
#   logistic: F logistic,               z = (t - b2) / b3
#   gompertz: F gompertz,               z = (t - b2) / b3
#   hill:     F logistic,               z = b3 (log t - log b2), f(t) = b1 / (1 + (b2 / t)^b3)
#   weibull:  F complementary log-log,  z = b3 (log t - log b2) + log((b3 - 1) / b3),
#             f(t) = b1 (1 - exp(-(b3 - 1) / b3 (t / b2)^b3))
#   normal:   F standard normal,        z = (t - b2) / b3
#   gamma:    F gamma with shape b3,    z = t (b3 - 1) / b2: the gamma distribution function with shape
#             b3 and scale b2 / (b3 - 1), whose density has its mode at b2
#   linear:   f(t) = b1 t + b2, for comparison; b3-b5 are not used, and it is not truncated
# hill, weibull and gamma are for positive times (logtime) and zero up to t = 0.
# b1-b5: role of each parameter, NA if the model does not use it (it is then always fixed)
#   b1: "asymptote", the response approached at late times; "slope" of the linear model
#   b2: "inflection", the time of the maximum growth rate, or "halfway", the time of half the
#       asymptote; "intercept" of the linear model
#   b3: "scale" on the time axis, positive; "slope" of the hill curve on the
#       log time axis; "shape" of the weibull and gamma curves, larger than 1
#   b4: not used by the curves themselves, see truncpars
#   b5: "asymmetry", positive, with b5 = 1 the symmetric curves above; fixed at 1 by default. The
#       curve is b1 G(z) with (powS, pow, shift):
#       logistic (Richards): G = (1 + b5 exp(-z))^(-1/b5) = F(z - log b5)^(1/b5)
#       hill, weibull, normal, gamma: G = F(z)^b5
#       gompertz: G = 1 - S(z)^b5, a power of the survival function (a power of F is only a shift)
# b4 is the threshold of the truncated curves (truncated=TRUE in gcTMB), for all models but the linear
# one: see truncpars.
# z(t, b2, b3), zinv(z, b2, b3): the transformed time and its inverse
# b2fromz(t, z, b3): b2 for which z(t, b2, b3) = z, for the profiles of the times
# line: the linearisation for the starting values, Finv(y / b1) = a + c h(t): h(t), the slope c of
#   b3 (cslope), and b2 and b3 from a and c (fromline); the gamma model, whose F depends on b3, is
#   linearised for a grid of shapes (shapeF)
# dens(z), F(z), Finv(q): the density and distribution function of the base distribution and the
#   inverse of F; logS(z): the log of its survival function 1 - F(z), Sinvlog(ls): the inverse, z with
#   logS(z) = ls; for the truncated curves, accurate for thresholds far beyond b2. The functions of the
#   gamma base take the shape as their second argument, the others ignore it.
basedists <- list(
  logistic=list(dens=function(z, ...) dlogis(z),
                F=function(z, ...) plogis(z),
                Finv=function(q, ...) qlogis(q),
                logF=function(z, ...) plogis(z, log.p=TRUE),
                Finvlog=function(lf, ...) qlogis(lf, log.p=TRUE),
                logS=function(z, ...) plogis(z, lower.tail=FALSE, log.p=TRUE),
                Sinvlog=function(ls, ...) qlogis(ls, lower.tail=FALSE, log.p=TRUE)),
  gompertz=list(dens=function(z, ...) exp(-z - exp(-z)),
                F=function(z, ...) exp(-exp(-z)),
                Finv=function(q, ...) -log(-log(q)),
                logF=function(z, ...) -exp(-z),
                Finvlog=function(lf, ...) -log(-lf),
                logS=function(z, ...) log(-expm1(-exp(-z))),
                Sinvlog=function(ls, ...) -log(-log1p(-exp(ls)))),
  cloglog=list(dens=function(z, ...) exp(z - exp(z)),
               F=function(z, ...) -expm1(-exp(z)),
               Finv=function(q, ...) log(-log1p(-q)),
               logF=function(z, ...) log(-expm1(-exp(z))),
               Finvlog=function(lf, ...) log(-log1p(-exp(lf))),
               logS=function(z, ...) -exp(z),
               Sinvlog=function(ls, ...) log(-ls)),
  normal=list(dens=function(z, ...) dnorm(z),
              F=function(z, ...) pnorm(z),
              Finv=function(q, ...) qnorm(q),
              logF=function(z, ...) pnorm(z, log.p=TRUE),
              Finvlog=function(lf, ...) qnorm(lf, log.p=TRUE),
              logS=function(z, ...) pnorm(z, lower.tail=FALSE, log.p=TRUE),
              Sinvlog=function(ls, ...) qnorm(ls, lower.tail=FALSE, log.p=TRUE)),
  gamma=list(dens=function(z, sh) dgamma(z, sh),
             F=function(z, sh) pgamma(z, sh),
             Finv=function(q, sh) qgamma(q, sh),
             logF=function(z, sh) pgamma(z, sh, log.p=TRUE),
             Finvlog=function(lf, sh) qgamma(lf, sh, log.p=TRUE),
             logS=function(z, sh) pgamma(z, sh, lower.tail=FALSE, log.p=TRUE),
             Sinvlog=function(ls, sh) qgamma(ls, sh, lower.tail=FALSE, log.p=TRUE))
)

# time axis of the location-scale models (logistic, gompertz, normal)
lstime <- list(logtime=FALSE,
               z=function(t, b2, b3) (t - b2) / b3,
               zinv=function(z, b2, b3) b2 + b3 * z,
               b2fromz=function(t, z, b3) t - b3 * z,
               h=function(t) t,
               cslope=function(b3) 1 / b3,
               fromline=function(a, c) c(-a / c, 1 / c))

# asymmetry by a power of F (see gcmodels)
powF <- list(powS=FALSE, pow=function(a) a, shift=function(a) 0)

gcmodels <- list(
  logistic=c(list(b1="asymptote", b2="inflection", b3="scale", b4=NA, b5="asymmetry", powS=FALSE,
                    pow=function(a) 1 / a, shift=function(a) log(a)), lstime, basedists$logistic),
  gompertz=c(list(b1="asymptote", b2="inflection", b3="scale", b4=NA, b5="asymmetry", powS=TRUE,
                    pow=function(a) a, shift=function(a) 0), lstime, basedists$gompertz),
  hill=c(list(b1="asymptote", b2="halfway", b3="slope", b4=NA, b5="asymmetry", logtime=TRUE,
              z=function(t, b2, b3) b3 * (log(t) - log(b2)),
              zinv=function(z, b2, b3) b2 * exp(z / b3),
              b2fromz=function(t, z, b3) t * exp(-z / b3),
              h=function(t) log(t),
              cslope=function(b3) b3,
              fromline=function(a, c) c(exp(-a / c), c)),
         powF, basedists$logistic),
  weibull=c(list(b1="asymptote", b2="inflection", b3="shape", b4=NA, b5="asymmetry", logtime=TRUE,
                 z=function(t, b2, b3) b3 * (log(t) - log(b2)) + log((b3 - 1) / b3),
                 zinv=function(z, b2, b3) b2 * exp((z - log((b3 - 1) / b3)) / b3),
                 b2fromz=function(t, z, b3) t * exp(-(z - log((b3 - 1) / b3)) / b3),
                 h=function(t) log(t),
                 cslope=function(b3) b3,
                 fromline=function(a, c) c(exp((log((c - 1) / c) - a) / c), c)),
            powF, basedists$cloglog),
  normal=c(list(b1="asymptote", b2="inflection", b3="scale", b4=NA, b5="asymmetry"), powF, lstime, basedists$normal),
  gamma=c(list(b1="asymptote", b2="inflection", b3="shape", b4=NA, b5="asymmetry", logtime=TRUE, shapeF=TRUE,
               z=function(t, b2, b3) t * (b3 - 1) / b2,
               zinv=function(z, b2, b3) z * b2 / (b3 - 1),
               b2fromz=function(t, z, b3) t * (b3 - 1) / z,
               h=function(t) t),
          powF, basedists$gamma),
  linear=list(b1="slope", b2="intercept", b3=NA, b4=NA, b5=NA, linear=TRUE, logtime=FALSE)
)

# links of the parameters b1-b5 of a model m (see gcmodels): the possible links of each parameter, the
# first one is the default. The scale b3 may have the inverse link, so that its linear predictor is
# the rate 1/b3. The log link keeps a parameter positive: the asymptote b1, the times b2 and b4 (the
# threshold of a truncated curve), and b3; the logit link keeps the asymptote b1 in (0, 1), for
# proportions.
modellinks <- function(m){
  list(b1=if (m$b1 %in% "asymptote") c("identity", "log", "logit") else if (m$b1 %in% "slope") c("identity", "log") else "identity",
       b2=if (m$b2 %in% c("inflection", "halfway")) c("identity", "log") else "identity",
       b3=if (m$b3 %in% "scale") c("identity", "inverse", "log") else if (m$b3 %in% c("slope", "shape")) c("identity", "log") else "identity",
       b4=c("identity", "log"),
       b5=c("identity", "log"))
}

# names of the links, in the order of their index in the TMB templates
linknames <- c("identity", "inverse", "log", "logit")

# index of each link, as in the TMB templates: 1 = identity, 2 = inverse, 3 = log, 4 = logit
linkindex <- function(link){
  match(link, linknames)
}

# parameters b1-b5 from their linear predictors eta (a vector, or a matrix with a column for each
# parameter), by the inverse links lnk (see linkindex)
linkpars <- function(lnk, eta){
  if (is.matrix(eta)) return(do.call(cbind, lapply(1:5, function(k) linkfun(lnk[k], eta[, k], inverse=TRUE))))
  sapply(1:5, function(k) linkfun(lnk[k], eta[k], inverse=TRUE))
}

# link function (inverse=FALSE) or inverse link of a parameter, by its index (see linkindex)
linkfun <- function(lnk, x, inverse=FALSE){
  if (inverse) switch(lnk, x, 1/x, exp(x), plogis(x)) else switch(lnk, x, 1/x, log(x), qlogis(x))
}

# derivative of the inverse link lnk at the linear predictor eta
dlinkinv <- function(lnk, eta){
  switch(lnk, rep(1, length(eta)), -1/eta^2, exp(eta), dlogis(eta))
}

# properties of a model, by name or index
gcmodel <- function(mod){
  gcmodels[[mod]]
}

# index of the model, as in the TMB template
modelindex <- function(model){
  which(model == names(gcmodels))
}

# TRUE for the linear model, which is not a growth curve with an asymptote
islinear <- function(mod){
  isTRUE(gcmodel(mod)$linear)
}

# roles of the parameters b1-b5 of a model (see gcmodels); truncated curves use b4 as their threshold
truncpars <- function(mod, trunc=FALSE){
  m <- gcmodel(mod)
  roles <- c(m$b1, m$b2, m$b3, m$b4, m$b5)
  if (trunc && !islinear(mod)) roles[4] <- "threshold"
  roles
}

# parameters not used by the model, which cannot be estimated
unusedpars <- function(mod, trunc=FALSE){
  which(is.na(truncpars(mod, trunc)))
}

# truncation of a fit: TRUE for a curve truncated at the threshold b4
fittrunc <- function(object){
  isTRUE(object$truncated)
}

# b1 of a fit as the maximum growth rate (b1type="rate") instead of the asymptote
fitrate <- function(object){
  identical(object$b1type, "rate")
}

# maximum growth rate of curves with asymptote 1 and parameters f2-f5 (vectors, recycled), as gcrate1
# in the TMB template: the asymptote of a curve whose b1 is the maximum growth rate is b1 / rate1
rate1 <- function(mod, f2, f3, f4, f5, trunc=FALSE){
  if (islinear(mod)) return(rep(1, max(length(f2), 1)))
  p <- cbind(f2, f3, f4, f5)
  key <- do.call(paste, c(as.data.frame(p), sep=""))
  first <- !duplicated(key)
  r <- vapply(which(first), function(i) maxrateR(mod, c(1, p[i, ]), trunc)[["rate"]], numeric(1))
  r[match(key, key[first])]
}

# parameters b1-b5 on their original scale from the linear predictors lpv of a curve of a fit, with b1
# the asymptote also for b1type="rate"
curvepars <- function(fit, lpv){
  b <- linkpars(linkindex(fit$link), lpv)
  if (fitrate(fit)) b[1] <- b[1] / rate1(modelindex(fit$model), b[2], b[3], b[4], b[5], fittrunc(fit))
  b
}

# asymmetric curve G of model m with the asymmetry a = b5 (see gcmodels), as in the TMB template:
# log G(z), log(1 - G(z)), and their inverses; sh = b3 is the shape of the gamma base
logG <- function(m, z, a, sh=NULL){
  if (m$powS) log(-expm1(a * m$logS(z, sh))) else m$pow(a) * m$logF(z - m$shift(a), sh)
}
logSG <- function(m, z, a, sh=NULL){
  if (m$powS) a * m$logS(z, sh) else log(-expm1(m$pow(a) * m$logF(z - m$shift(a), sh)))
}
SGinvlog <- function(m, ls, a, sh=NULL){
  if (m$powS) m$Sinvlog(ls / a, sh) else m$shift(a) + m$Finvlog(log(-expm1(ls)) / m$pow(a), sh)
}
Ginv <- function(m, q, a, sh=NULL){
  if (m$powS) m$Sinvlog(log1p(-q) / a, sh) else m$shift(a) + m$Finvlog(log(q) / m$pow(a), sh)
}

# log(1 - G) at the threshold f4 of a truncated curve of model m; for the models on positive times, a
# threshold up to zero does not truncate the curve
logS4 <- function(m, f2, f3, f4, f5=1){
  if (!m$logtime) return(logSG(m, m$z(f4, f2, f3), f5, f3))
  ifelse(f4 > 0, logSG(m, m$z(ifelse(f4 > 0, f4, 1), f2, f3), f5, f3), 0)
}

# growth curve, same as gcmean in the TMB template: f1 G(z(t)) with the asymmetry f5, and truncated
# at the threshold f4 (trunc=TRUE), f1 (G(z) - G(z4)) / (1 - G(z4)) for t > f4 and zero otherwise,
# z4 = z(f4). The models on positive times are zero up to t = 0. The linear model is f1 t + f2.
gcmeanR <- function(mod, t, f1, f2, f3, f4, f5=1, trunc=FALSE, b1rate=FALSE){
  m <- gcmodel(mod)
  if (isTRUE(m$linear)) return(f1 * t + f2)
  # b1 as the maximum growth rate: the asymptote
  if (b1rate) f1 <- f1 / rate1(mod, f2, f3, f4, f5, trunc)
  pos <- if (m$logtime) t > 0 else rep(TRUE, length(t))
  z <- m$z(ifelse(pos, t, 1), f2, f3)
  Ft <- if (!trunc) exp(logG(m, z, f5, f3)) else ifelse(t > f4, -expm1(logSG(m, z, f5, f3) - logS4(m, f2, f3, f4, f5)), 0)
  ifelse(pos, Ft, 0) * f1
}

# time at which a curve with parameters b2, b3 and asymmetry b5 reaches the fraction q of its
# asymptote, as gctime in the TMB template (residual zero): zinv(Ginv(q)); truncated at the threshold
# b4 (trunc=TRUE), zinv(SGinv((1 - q) S(z(b4)))). For the linear model, the time at which it reaches
# the response q b1.
gctimeR <- function(mod, q, b2, b3, b4=NA, b5=1, trunc=FALSE, b1=1){
  m <- gcmodel(mod)
  if (isTRUE(m$linear)) return(q - b2 / b1)
  z <- if (!trunc) Ginv(m, q, b5, b3) else SGinvlog(m, log1p(-q) + logS4(m, b2, b3, b4, b5), b5, b3)
  m$zinv(z, b2, b3)
}
