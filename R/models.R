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
#   chapmanrichards, wykoff: b4 + b1 G(t), with the intercept b4 (startcurve, wykoff entry); never
#             truncated
#   beta:     zero up to the start b4 (beta entry), so that a truncation at b4 changes nothing
#   preecebaines: model 1 of Preece and Baines (1978), b1 - 2 (b1 - b3) / (exp(b4 (t - b2)) +
#             exp(b5 (t - b2))), never truncated
#   These four have the curve G(t, b2, b3, b4) relative to b1 (or the mean function mean for
#   preecebaines), the inverse tq and the maximum growth rate mr instead of a base distribution and a
#   transformed time (start=TRUE). The position of a model in gcmodels is its index in the TMB
#   templates; modelnames gives the order in which the models are listed to the user.
# hill, weibull and gamma are for positive times (logtime) and zero up to t = 0.
# b1-b5: role of each parameter, NA if the model does not use it (it is then always fixed)
#   b1: "asymptote", the response approached at late times; "slope" of the linear model
#   b2: "inflection", the time of the maximum growth rate, or "halfway", the time of half the
#       asymptote; "intercept" of the linear model
#   b3: "scale" on the time axis, positive; "slope" of the hill curve on the
#       log time axis; "shape" of the weibull and gamma curves, larger than 1; "end" of growth of beta
#   b4: not used by the curves themselves, see truncpars; "intercept" of chapmanrichards and wykoff,
#       b4 + b1 G(t), e.g. the breast height of height-diameter models, fixed at 0 by default
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

# the chapmanrichards curve G(t) = (1 - exp(-t/b3))^c, c = exp(b2/b3), zero up to t = 0 (see
# gcmodels): the fraction G of the asymptote at t, the time tq at which it reaches the fraction q, and
# the time t, fraction G and rate r (asymptote 1) of the maximum growth rate, at the inflection point
# (at t = 0 for curves without an inflection point at a positive time)
startcurve <- function(){
  cfun <- function(b2, b3) exp(b2 / b3)
  list(start=TRUE, logtime=FALSE, cfun=cfun,
       G=function(t, b2, b3, b4){
         z <- pmax(t / b3, 0)
         ifelse(t > 0, (-expm1(-z))^cfun(b2, b3), 0)
       },
       tq=function(q, b2, b3, b4, ...) -b3 * log1p(-q^(1 / cfun(b2, b3))),
       mr=function(b2, b3, b4, b5, b1=1){
         cc <- cfun(b2, b3)
         if (!all(is.finite(c(cc, b3))) || b3 <= 0) return(c(t=NA, G=NA, r=NA))
         if (cc > 1) return(c(t=b3 * log(cc), G=(1 - 1 / cc)^cc, r=b1 * (1 - 1 / cc)^(cc - 1) / b3))
         c(t=0, G=0, r=if (cc == 1) b1 / b3 else Inf)
       })
}

# the beta growth function of Yin et al. (2003), eqn (8): the maximum b1 reached at the end of growth
# te = b3, the time b2 = tm of the maximum growth rate and the start b4 = tb, zero up to b4. With
# u = (t - b4) / D, D = b3 - b4, and d = (b3 - b2) / D, the curve relative to b1 is
# G(u) = (1 + (1 - u) / d) u^(1 / d), which rises from 0 at u = 0 to 1 at u = 1 (te). Capped, it stays
# at 1 from te on (the beta sigmoid growth function of Yin et al.); otherwise it falls back to zero at
# u = 1 + d, the time 2 b3 - b2, and is zero after it. The levels of timeto and the maximum growth rate
# are those of the rising part, up to te.
betamodel <- function(capped){
  list(b1="asymptote", b2="inflection", b3="end", b4=NA, b5=NA, start=TRUE, logtime=FALSE, capped=capped,
       G=function(t, b2, b3, b4){
         D <- b3 - b4
         d <- (b3 - b2) / D
         u <- pmin(pmax((t - b4) / D, 0), if (capped) 1 else 1 + d)
         ifelse(t > b4 & (capped | u < 1 + d), (1 + (1 - u) / d) * u^(1 / d), 0)
       },
       tq=function(q, b2, b3, b4, ...){
         D <- b3 - b4
         d <- (b3 - b2) / D
         u <- vapply(q, function(qi) if (!is.finite(qi) || qi <= 0 || qi >= 1 || !is.finite(d) || d <= 0 || d >= 1) NA_real_ else
           uniroot(function(u) (1 + (1 - u) / d) * u^(1 / d) - qi, c(0, 1), tol=1e-14)$root, numeric(1))
         b4 + D * u
       },
       mr=function(b2, b3, b4, b5, b1=1){
         D <- b3 - b4
         d <- (b3 - b2) / D
         if (!all(is.finite(c(D, d))) || d <= 0 || d >= 1) return(c(t=NA, G=NA, r=NA))
         c(t=b2, G=2 * (1 - d)^(1 / d), r=b1 * (1 + 1 / d) * (1 - d)^(1 / d - 1) / D)
       })
}

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
  linear=list(b1="slope", b2="intercept", b3=NA, b4=NA, b5=NA, linear=TRUE, logtime=FALSE),
  chapmanrichards=c(list(b1="asymptote", b2="inflection", b3="scale", b4="intercept", b5=NA), startcurve()),
  # Yin et al. (2003), capped at the maximum b1 from the end of growth b3 on, see betamodel
  beta_capped=betamodel(capped=TRUE),
  # Wykoff et al. (1982), exp(a + b / (t + 1)): b1 exp(-s / (t + 1)) with s = 2 (b2 + 1), which puts
  # the inflection point at b2, so that b1 = exp(a) and s = -b; positive from t = -1 on
  wykoff=list(b1="asymptote", b2="inflection", b3=NA, b4="intercept", b5=NA, start=TRUE, logtime=FALSE,
              G=function(t, b2, b3, b4){
                u <- t + 1
                ifelse(u > 0, exp(-2 * (b2 + 1) / pmax(u, 1e-300)), 0)
              },
              tq=function(q, b2, b3, b4, ...) -1 - 2 * (b2 + 1) / log(q),
              mr=function(b2, b3, b4, b5, b1=1){
                s <- 2 * (b2 + 1)
                if (!is.finite(s) || s <= 0) return(c(t=NA, G=NA, r=NA))
                c(t=b2, G=exp(-2), r=b1 * 4 * exp(-2) / s)
              }),
  # Preece and Baines (1978), model 1: the adult height b1, the age b2 of the middle of puberty, the
  # height b3 at b2, and the rates b4 < b5 before and during puberty; estimated by default (free)
  preecebaines=list(b1="asymptote", b2="inflection", b3="level", b4="rate", b5="rate", start=TRUE, logtime=FALSE,
                    trunc=FALSE, free=4:5,
                    mean=function(t, b1, b2, b3, b4, b5){
                      u <- t - b2
                      b1 - 2 * (b1 - b3) * exp(-pblse(b4 * u, b5 * u))
                    },
                    # the response q b1: log(exp(b4 u) + exp(b5 u)) = log(2 (b1 - b3) / (b1 (1 - q))),
                    # increasing in u
                    tq=function(q, b2, b3, b4, b5, b1=1){
                      L <- log(2 * (b1 - b3) / (b1 * (1 - q)))
                      u <- vapply(L, function(l) if (!is.finite(l) || !all(is.finite(c(b4, b5))) || b4 <= 0 || b5 <= 0) NA_real_ else
                        uniroot(function(u) pblse(b4 * u, b5 * u) - l, c(-1, 1), extendInt="upX", tol=1e-14)$root, numeric(1))
                      b2 + u
                    },
                    # the pubertal peak of the growth rate: with a = min(b4, b5), b = max(b4, b5) and
                    # w = exp((b - a) u), the larger root of b^2 w^2 - (a^2 - 4ab + b^2) w + a^2 = 0
                    # (the smaller one is the minimum before it, the take-off); none for b/a <= 3 + 2 sqrt(2)
                    mr=function(b2, b3, b4, b5, b1=1){
                      a <- min(b4, b5)
                      b <- max(b4, b5)
                      k <- a^2 - 4 * a * b + b^2
                      disc <- k^2 - 4 * a^2 * b^2
                      if (!all(is.finite(c(a, b, k, disc))) || a <= 0 || k <= 0 || disc <= 0) return(c(t=NA, G=NA, r=NA))
                      w <- (k + sqrt(disc)) / (2 * b^2)
                      u <- log(w) / (b - a)
                      f <- b1 - 2 * (b1 - b3) * exp(-pblse(b4 * u, b5 * u))
                      r <- 2 * (b1 - b3) * (b4 * exp(b4 * u) + b5 * exp(b5 * u)) / (exp(b4 * u) + exp(b5 * u))^2
                      c(t=b2 + u, G=f / b1, r=r)
                    }),
  # Yin et al. (2003), eqn (8) beyond the end of growth b3, back to zero at 2 b3 - b2, see betamodel
  beta=betamodel(capped=FALSE)
)

# the models in the order in which they are listed to the user
modelnames <- c("logistic", "gompertz", "hill", "weibull", "beta", "beta_capped", "normal", "gamma", "chapmanrichards",
                "wykoff", "preecebaines", "linear")

# log(exp(x) + exp(y)), for the Preece-Baines curve
pblse <- function(x, y){
  m <- pmax(x, y)
  m + log(exp(x - m) + exp(y - m))
}

# links of the parameters b1-b5 of a model m (see gcmodels): the possible links of each parameter, the
# first one is the default. The scale b3 may have the inverse link, so that its linear predictor is
# the rate 1/b3. The log link keeps a parameter positive: the asymptote b1, the times b2 and b4 (the
# threshold of a truncated curve), and b3; the logit link keeps the asymptote b1 in (0, 1), for
# proportions.
modellinks <- function(m){
  list(b1=if (m$b1 %in% "asymptote") c("identity", "log", "logit") else if (m$b1 %in% "slope") c("identity", "log") else "identity",
       b2=if (m$b2 %in% c("inflection", "halfway")) c("identity", "log") else "identity",
       b3=if (m$b3 %in% "scale") c("identity", "inverse", "log") else if (m$b3 %in% c("slope", "shape", "end", "level")) c("identity", "log") else "identity",
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

# TRUE for chapmanrichards, beta and wykoff, with their own curves G (see gcmodels), which are never
# truncated: beta is zero up to its start b4, and b4 is the intercept of the other two
isstart <- function(mod){
  isTRUE(gcmodel(mod)$start)
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
  key <- do.call(paste, c(as.data.frame(p), sep=" "))
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

# intercept of the curves of model mod: b4 for chapmanrichards and wykoff, zero otherwise
intercept <- function(mod, b4){
  if (identical(gcmodel(mod)$b4, "intercept")) b4 else 0
}

# TRUE for the models with the intercept b4 (chapmanrichards, wykoff), which are never truncated
hasintercept <- function(mod){
  identical(gcmodel(mod)$b4, "intercept")
}

# growth curve, same as gcmean in the TMB template: f1 G(z(t)) with the asymmetry f5, and truncated
# at the threshold f4 (trunc=TRUE), f1 (G(z) - G(z4)) / (1 - G(z4)) for t > f4 and zero otherwise,
# z4 = z(f4). The models on positive times are zero up to t = 0. The linear model is f1 t + f2.
# chapmanrichards and wykoff are f4 + f1 G(t), with the intercept f4; beta is zero up to its start f4.
gcmeanR <- function(mod, t, f1, f2, f3, f4, f5=1, trunc=FALSE, b1rate=FALSE){
  m <- gcmodel(mod)
  if (isTRUE(m$linear)) return(f1 * t + f2)
  # b1 as the maximum growth rate: the asymptote
  if (b1rate) f1 <- f1 / rate1(mod, f2, f3, f4, f5, trunc)
  # chapmanrichards, beta and wykoff: never truncated (beta is zero up to b4 anyway)
  if (!is.null(m$mean)) return(m$mean(t, f1, f2, f3, f4, f5))
  if (isTRUE(m$start)) return(intercept(mod, f4) + f1 * m$G(t, f2, f3, f4))
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
  if (isTRUE(m$start)) return(m$tq(q, b2, b3, b4, b5, b1=b1))
  z <- if (!trunc) Ginv(m, q, b5, b3) else SGinvlog(m, log1p(-q) + logS4(m, b2, b3, b4, b5), b5, b3)
  m$zinv(z, b2, b3)
}
