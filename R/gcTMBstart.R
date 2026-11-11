#' Starting values for a gcTMB model
#'
#' Computes starting values for the fixed-effect parameters of a growth curve,
#' which can be inspected, modified, and passed to the \code{start} argument of \code{\link{gcTMB}}.
#'
#' Method "lm" linearises the growth curve: for a given asymptote b1, the mean response at each
#' time t, relative to b1, gives Finv(y/b1) = (t - b2)/b3, with the inverse Finv of the function F
#' of the model (qlogis for logistic, -log(-log(q)) for gompertz). A linear regression of Finv(y/b1)
#' on the time gives b3 as the inverse of the slope and b2 as the time at which the line crosses zero;
#' with b2 or b3 given in \code{start}, only the other one is estimated. The observed responses
#' often do not reach the asymptote yet, so b1 is profiled: for asymptotes from just above the
#' largest mean response to four times it, the linearised fit gives b2 and b3, b1 is refitted by least
#' squares for this curve, and the asymptote with the smallest residual sum of squares is used.
#' Mean responses below 1\% or above 99\% of the asymptote carry little information about the
#' linearised curve and are left out. With covariates in \code{fform}, the pooled starting value of
#' each parameter is projected onto its design matrix, giving the intercept or each group the pooled
#' value and contrasts a value of zero. The reserved parameters b4 and b5 start at 0.
#'
#' Method "grid" evaluates the negative log-likelihood (without random effects) at all
#' combinations of a grid of parameter values and returns the best grid point; parameters
#' not included in the grid start at the values of method "lm".
#' Without bounds in \code{gsl} and \code{gsu}, default ranges are used for b2 (the time range
#' extended by a quarter at both ends) and b3 (the sign of the "lm" starting value, with 1/4 to 4
#' times its value, log-spaced; on the scale of its link). Default ranges are only used for parameters
#' with a single coefficient in \code{fform}; the asymptote b1 is only included in the grid with
#' user-specified bounds.
#'
#' @param form Formula, with the names of the response and the time variable, e.g. \code{weight ~ day}
#' @param fform Formula, list of formulas, or NULL. Fixed effects for the model parameters: a formula of the form b1 + b2 ~ Treatment, or a list of such formulas, e.g. list(b1 ~ Treatment, b2 ~ Treatment + Site). Parameters not listed have a single intercept (~ 1). A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~ 1), is also accepted.
#' @param data Data frame. Observations with missing values in the variables of the model or in the weights are left out.
#' @param start Numeric vector or NULL. Starting values for each model parameter. Only NA entries are replaced.
#' @param fix Logical vector or NULL. A vector with the same length as start. If TRUE, a model coefficient is fixed to the starting value. The reserved parameters b4 and b5, which the models do not use, are always fixed, at 0 unless given in start.
#' @param weights Numeric vector or NULL. Weights, one for each row of data.
#' @param family Character string. The distribution of the response, see \code{\link{gcTMB}}.
#' @param model Character string. The growth model: "logistic" or "gompertz".
#' @param link Named character vector or NULL. Links of the parameters b1-b5, see \code{\link{gcTMB}}.
#' @param method Character string. One of "lm" or "grid".
#' @param gsl Numeric vector or NULL. Lower bounds for the grid search, same length as start, on the scale of the coefficients (the rate 1/b3 for b3 with the inverse link); NA (or NULL) for the default range.
#' @param gsu Numeric vector or NULL. Upper bounds for the grid search, same length as start; NA (or NULL) for the default range.
#' @param gsn Integer, a single number or a vector, or NULL. Number of grid points for each parameter, a single number or a vector with the same length as start (default 10). A parameter with gsn = 1 stays at its "lm" starting value.
#'
#' @return A named numeric vector of starting values, one for each model coefficient.
#' @examples
#' gcTMBstart(height ~ DAS, data=FabaHeights25, model="gompertz")
#' # with a given asymptote, from which b2 and b3 are linearised
#' gcTMBstart(height ~ DAS, data=FabaHeights25, start=c(b1=120))
#' # best point of a grid of b2 and b3 values
#' gcTMBstart(height ~ DAS, data=FabaHeights25, method="grid", gsn=c(1, 10, 10, 1, 1))
#' @export
#'
gcTMBstart <- function(form, fform=NULL, data, start=NULL, fix=NULL, weights=NULL, family="gaussian", model="logistic", link=NULL, method="lm", gsl=NULL, gsu=NULL, gsn=NULL){
  dd <- gcTMBdata(form, fform=fform, data=data, weights=weights, family=family, model=model, link=link)
  start <- parvector(start, dd$npi, NA, "start")
  fix <- parvector(fix, dd$npi, FALSE, "fix")
  fix[dd$npi %in% unusedpars(dd$mod)] <- TRUE
  sv <- gcTMBstartvals(dd, start, fix, method=method, gsl=gsl, gsu=gsu, gsn=gsn)
  names(sv) <- coefnames(dd$Xs)
  sv
}


# model data shared by gcTMB and gcTMBstart
gcTMBdata <- function(form, fform=NULL, rform=NULL, data, weights=NULL, family="gaussian", model="logistic", link=NULL){
  fform <- fforms(fform)
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (length(weights) != nrow(data)) stop("weights must have one entry for each row of data.")

  # observations with missing values in the variables of the model or in the weights are left out
  rfl <- if (inherits(rform, "formula")) list(rform) else rform
  vars <- c(all.vars(form), unlist(lapply(fform, all.vars)), unlist(lapply(rfl, function(f) all.vars(f[[3]]))))
  keep <- complete.cases(data[, intersect(vars, names(data)), drop=FALSE]) & !is.na(weights)
  data <- data[keep, , drop=FALSE]
  weights <- weights[keep]

  # response and time
  mf <- model.frame(form, data=data)
  if (ncol(mf) != 2 || !is.numeric(mf[, 2])) stop("form must have the response on the left and a single numeric time variable on the right, e.g. weight ~ day.")
  # fixed effects; terms, factor levels and contrasts are kept for the design matrices of new data
  fterms <- lapply(fform, function(f){
    fmf <- model.frame(f, data=data)
    tt <- terms(fmf)
    list(terms=tt, xlev=.getXlevels(tt, fmf), contrasts=attr(model.matrix(tt, fmf), "contrasts"))
  })
  Xs <- lapply(fterms, function(ft) designmatrix(ft, data))
  npars <- sapply(Xs, ncol)

  # model
  mod <- modelindex(model)
  fam <- which(family == c("gaussian", "poisson", "binomial", "beta", "nbinom2", "nbinom1", "betabinomial", "lognormal"))
  if (length(mod) != 1) stop("Unknown model: ", model, "; one of ", paste(names(gcmodels), collapse=", "), ".")
  if (length(fam) != 1) stop("Unknown family: ", family, "; one of gaussian, lognormal, poisson, binomial, betabinomial, beta, nbinom2, nbinom1.")
  m <- gcmodel(mod)
  link <- parlinks(link, m, model)
  lnk <- linkindex(link)

  # response
  y <- mr <- model.response(mf, "any")
  if (family %in% c("binomial", "betabinomial") & ncol(as.matrix(mr)) == 2){
    bn <- mr[,1] + mr[,2]
    y <- mr[,1]
  } else if (family == "betabinomial"){
    stop("the precision of the beta-binomial family is not identifiable for binary responses; give the response as cbind(successes, failures).")
  } else if (family == "binomial"){
    # binary response: one trial per observation
    bn <- rep(1, length(y))
  } else {
    bn <- rep(0, length(y))
  }
  if (family == "lognormal" && any(!is.na(y) & y <= 0)) stop("the lognormal family needs a positive response.")
  if (family == "beta" && any(!is.na(y) & (y <= 0 | y >= 1))) stop("the beta family needs a response strictly between 0 and 1.")

  # random effects
  re <- reterms(rform, data)
  rpars <- unlist(lapply(re$terms, function(tm) tm$pars))
  unused <- intersect(unusedpars(mod), rpars)
  if (length(unused) > 0) stop("The ", model, " model does not use ", paste("b", unused, sep="", collapse=", "), "; remove it from rform.")

  list(y=y, x=mf[,2], w=weights, bn=bn, fam=fam, mod=mod, lnk=lnk, re=re,
       Xs=Xs, fform=fform, fterms=fterms, npars=npars, npi=rep(1:5, npars),
       family=family, model=model, link=link, data=data)
}

# design matrix of a parameter for data, with the terms, factor levels and contrasts of the fit
designmatrix <- function(ft, data){
  # xlev drops contrasts set on a factor of data, with a warning; contrasts.arg restores them
  mf <- withCallingHandlers(model.frame(ft$terms, data=data, xlev=ft$xlev, na.action=na.pass),
                            warning=function(w) if (grepl("^contrasts dropped", conditionMessage(w))) invokeRestart("muffleWarning"))
  model.matrix(ft$terms, mf, contrasts.arg=ft$contrasts)
}


# starting values for all fixed-effect parameters; only NA entries of start are replaced
gcTMBstartvals <- function(dd, start, fix, method="lm", gsl=NULL, gsu=NULL, gsn=NULL){
  method <- match.arg(method, c("lm", "grid"))
  sv <- lmstart(dd, start, fix)
  if (method == "grid") sv <- gridstart(dd, sv, fix, gsl, gsu, gsn)
  sv
}


# starting values from the linearised growth curve, Finv(y / b1) = (t - b2) / b3, profiled over the
# asymptote b1
lmstart <- function(dd, start, fix){
  npi <- dd$npi
  m <- gcmodel(dd$mod)
  # parameters with a given starting value, and their values on the original scale
  kn <- sapply(1:3, function(i) fix[npi == i][1] | !is.na(start[npi == i][1]))
  kv <- sapply(1:3, function(i) if (kn[i]) linkfun(dd$lnk[i], start[npi == i][1], inverse=TRUE) else NA)

  # response on the mean scale, aggregated by time; geometric means for the medians of lognormal responses
  if (dd$family %in% c("binomial", "betabinomial") & any(dd$bn > 0)){
    yd <- tapply(dd$w * dd$y, dd$x, sum) / tapply(dd$w * dd$bn, dd$x, sum)
  } else if (identical(dd$family, "lognormal")){
    yd <- exp(tapply(dd$w * log(dd$y), dd$x, sum) / tapply(dd$w, dd$x, sum))
  } else {
    yd <- tapply(dd$w * dd$y, dd$x, sum) / tapply(dd$w, dd$x, sum)
  }
  td <- sort(unique(dd$x))
  yd <- as.vector(yd)
  ok <- is.finite(yd)
  td <- td[ok]
  yd <- yd[ok]
  if (length(td) == 0) stop("no finite responses for the starting values.")
  trange <- if (length(td) > 1) diff(range(td)) else 1
  increasing <- length(td) < 2 || yd[length(yd)] >= yd[1]

  # b2 and b3 from the linearised curve for a given asymptote b1; fallbacks for too few points
  linfit <- function(b1){
    q <- yd / b1
    use <- is.finite(q) & q > 0.01 & q < 0.99
    t <- td[use]
    z <- m$Finv(q[use])
    b2 <- kv[2]
    b3 <- kv[3]
    if (!kn[2] && kn[3]){
      if (length(t) > 0) b2 <- mean(t - b3 * z)
    } else if (kn[2] && !kn[3]){
      # z = (t - b2) / b3: slope 1/b3 through the origin
      if (length(t) > 0 && any(t != b2)) b3 <- 1 / (sum((t - b2) * z) / sum((t - b2)^2))
    } else if (!kn[2] && !kn[3] && length(unique(t)) >= 2){
      cf <- coefficients(lm(z ~ t))
      b3 <- 1 / cf[[2]]
      b2 <- -cf[[1]] * b3
    }
    # fallbacks: a scale of an eighth of the time range in the direction of the data, and the time
    # at which the response is closest to half the asymptote
    if (!is.finite(b3) || b3 == 0) b3 <- (if (increasing) 1 else -1) * trange / 8
    if (!is.finite(b2)) b2 <- td[which.min(abs(yd / b1 - 0.5))]
    c(b2, b3)
  }

  # candidate asymptotes, from just above the response furthest from zero to four times it
  ymax <- yd[which.max(abs(yd))]
  if (ymax == 0) ymax <- 1
  b1grid <- if (kn[1]) kv[1] else ymax * exp(seq(log(1.01), log(4), length.out=40))
  b1grid <- vapply(b1grid, function(b) linkrange(b, dd$family), numeric(1))
  cand <- lapply(unique(b1grid), function(b1){
    p <- linfit(b1)
    Fz <- m$F((td - p[1]) / p[2])
    # the asymptote by least squares for this curve, unless given
    if (!kn[1] && sum(Fz^2) > 0) b1 <- linkrange(sum(yd * Fz) / sum(Fz^2), dd$family)
    list(b=c(b1, p), rss=sum((yd - b1 * Fz)^2))
  })
  rss <- vapply(cand, function(cd) cd$rss, numeric(1))
  best <- cand[[if (all(!is.finite(rss))) 1 else which.min(rss)]]$b

  pooled <- c(linkfun(dd$lnk[1], best[1]), linkfun(dd$lnk[2], best[2]), linkfun(dd$lnk[3], best[3]), 0, 0)
  projectstart(dd, start, fix, pooled)
}

# starting values from pooled values of b1-b5: projected onto each design matrix, for the NA
# entries of start that are not fixed
projectstart <- function(dd, start, fix, pooled){
  npi <- dd$npi
  for (i in 1:5){
    id <- which(npi == i)
    X <- dd$Xs[[i]]
    cf <- qr.coef(qr(X), rep(pooled[i], nrow(X)))
    cf[is.na(cf) | abs(cf) < 1e-8 * max(abs(pooled[i]), 1)] <- 0
    fill <- is.na(start[id]) & !fix[id]
    start[id[fill]] <- cf[fill]
  }
  start[is.na(start)] <- 0
  as.vector(start)
}


# best point of a grid of parameter values, evaluated without random effects
gridstart <- function(dd, start, fix, gsl, gsu, gsn){
  np <- length(start)
  if (is.null(gsl)) gsl <- rep(NA, np)
  if (is.null(gsu)) gsu <- rep(NA, np)
  if (is.null(gsn)) gsn <- 10
  if (length(gsl) != np | length(gsu) != np) stop("gsl and gsu must have the same length as start, with NA for default ranges.")
  if (length(gsn) == 1) gsn <- rep(gsn, np)
  if (length(gsn) != np) stop("gsn must be a single number or have the same length as start.")
  gsn[is.na(gsn)] <- 10

  # grid values for each parameter: user bounds, otherwise default ranges
  dflt <- gridrange(dd, start)
  gvals <- lapply(1:np, function(i){
    if (fix[i] | gsn[i] < 2) return(NULL)
    if (!is.na(gsl[i]) & !is.na(gsu[i])) return(seq(gsl[i], gsu[i], length.out=gsn[i]))
    if (is.null(dflt[[i]])) return(NULL)
    r <- dflt[[i]]
    if (!is.na(gsl[i])) r$lower <- gsl[i]
    if (!is.na(gsu[i])) r$upper <- gsu[i]
    if (r$log) sign(r$lower) * exp(seq(log(abs(r$lower)), log(abs(r$upper)), length.out=gsn[i]))
    else seq(r$lower, r$upper, length.out=gsn[i])
  })
  gid <- which(!sapply(gvals, is.null))
  if (length(gid) == 0) return(start)

  obj <- gcTMBobj(dd, start, fix, rform=NULL)
  # position of the grid parameters in obj$par (fixed parameters are mapped out)
  ppos <- cumsum(!fix)[gid]
  grid <- expand.grid(gvals[gid])
  objvals <- apply(grid, 1, function(g){
    par <- obj$par
    par[ppos] <- g
    obj$fn(par)
  })
  start[gid] <- as.numeric(grid[which.min(objvals), ])
  start
}
# default grid ranges for b2 and b3 from the data and the "lm" starting values (on the scale of
# their links); only for parameters with a single coefficient, and not for the asymptote b1
gridrange <- function(dd, start){
  npi <- dd$npi
  x <- dd$x
  rng <- vector("list", length(start))
  single <- function(i) if (dd$npars[i] == 1) which(npi == i) else NULL

  # b2: time of the maximum growth rate within the time range, extended at both ends
  i2 <- single(2)
  if (!is.null(i2)){
    ext <- 0.25*diff(range(x))
    if (ext == 0) ext <- 1
    rng[[i2]] <- list(lower=min(x) - ext, upper=max(x) + ext, log=FALSE)
  }
  # b3: keep the sign of the starting value, between 1/4 and 4 times its value
  i3 <- single(3)
  if (!is.null(i3)){
    b3 <- start[i3]
    if (!is.finite(b3) || b3 == 0) b3 <- 1
    rng[[i3]] <- list(lower=b3/4, upper=b3*4, log=TRUE)
  }
  rng
}


# TMB objective function for given starting values
gcTMBobj <- function(dd, start, fix, rform=NULL){
  npi <- dd$npi
  sb <- lapply(1:5, function(i) start[npi == i])
  re <- if (is.null(rform)) reterms(NULL, length(dd$y)) else dd$re
  tnp <- vapply(re$terms, function(tm) length(tm$pars), integer(1))
  tnl <- vapply(re$terms, function(tm) length(tm$levels), integer(1))
  # log standard deviations and transformed correlations for each term, from fits of each group;
  # random effects at zero
  vs <- if (length(re$terms) > 0) varstart(dd, start, re) else list(theta=numeric(0))
  stheta <- vs$theta
  up <- rep(0, sum(tnp*tnl))
  dlist <- list(model="gcTMB", y=dd$y, x=dd$x, w=dd$w, bn=dd$bn, fam=dd$fam, mod=dd$mod, lnk=as.integer(dd$lnk),
                tnp=tnp, tnl=tnl, tpar=as.integer(unlist(lapply(re$terms, function(tm) tm$pars)) - 1),
                tzb=if (length(re$terms) > 0) tzbvec(re$terms) else integer(0),
                simre=1L, Z=re$Z)
  # no predictions: see predictse
  Xp <- setNames(lapply(dd$Xs, function(X) X[0, , drop=FALSE]), paste("Xp", 1:5, sep=""))
  datalist <- c(dlist, dd$Xs, list(xp=numeric(0)), Xp, list(Zp=re$Z[0, , drop=FALSE]))

  # log_sigma: residual sd (gaussian, of log(y) for lognormal), precision (beta, betabinomial), size theta (nbinom2) or dispersion phi
  # (nbinom1); fixed for poisson and binomial
  if (dd$family == "gaussian") slsigma <- if (is.null(rform)) 1 else vs$log_sigma
  # lognormal: sd of log(y), an upper bound of the residual sd
  if (dd$family == "lognormal") slsigma <- if (!is.null(rform)) vs$log_sigma else log(max(sqrt(sum(dd$w * (log(dd$y) - sum(dd$w * log(dd$y))/sum(dd$w))^2)/sum(dd$w)), 1e-3))
  if (dd$family == "beta") slsigma <- 10
  if (dd$family %in% c("nbinom2", "nbinom1")) slsigma <- 0
  if (dd$family == "betabinomial") slsigma <- 2
  if (dd$family %in% c("poisson", "binomial")) slsigma <- 0
  parameters <- list(b1=sb[[1]], b2=sb[[2]], b3=sb[[3]], b4=sb[[4]], b5=sb[[5]], u=up, theta=as.numeric(stheta), log_sigma=slsigma)

  # fix parameters to starting values
  lmap <- list()
  if (any(fix)){
    idp <- 1:length(start)
    idp[fix] <- NA
    for (i in 1:5){
      lmap[[paste("b", i, sep="")]] <- as.factor(idp[npi == i])
    }
  }
  if (dd$family %in% c("poisson", "binomial")) lmap$log_sigma <- factor(NA)
  if (!(dd$family %in% c("poisson", "binomial"))) parameters$log_sigma <- sigmastart(datalist, parameters, lmap, rform)

  MakeADFun(datalist, parameters, map=lmap, random=if (!is.null(rform)) "u", DLL="gcTMB_TMBExports", silent=TRUE)
}

# starting value of log_sigma (residual sd, precision, size or dispersion), on the scale of the data:
# the maximum likelihood estimate with the coefficients and the variances of the random effects at
# their starting values, starting from the family default; with random effects, the random effects are
# integrated out (Laplace approximation), so that the variation between groups is left to them. The
# family default is kept if the fit fails or is not better.
sigmastart <- function(datalist, parameters, lmap, rform){
  ls0 <- parameters$log_sigma
  map1 <- lmap
  for (i in 1:5){
    b <- paste("b", i, sep="")
    map1[[b]] <- factor(rep(NA, length(parameters[[b]])))
  }
  map1$theta <- factor(rep(NA, length(parameters$theta)))
  if (is.null(rform)) map1$u <- factor(rep(NA, length(parameters$u)))
  res <- tryCatch({
    obj1 <- MakeADFun(datalist, parameters, map=map1, random=if (!is.null(rform)) "u", DLL="gcTMB_TMBExports", silent=TRUE)
    f0 <- obj1$fn(ls0)
    o <- suppressWarnings(nlminb(ls0, obj1$fn, obj1$gr, lower=ls0 - 20, upper=ls0 + 20))
    if (is.finite(o$objective) && (!is.finite(f0) || o$objective < f0)) o$par else ls0
  }, error=function(e) ls0)
  unname(res)
}


# starting values of the variance parameters of the random effect terms re, on the scale of the data:
# with the curve at the starting values, the deviations of each group are estimated by one Gauss-Newton
# step, a weighted least squares fit of the residuals of the group on the derivatives of the mean
# w.r.t. the linear predictors of the parameters with random coefficients (times the covariate of
# a random slope); exact for random asymptotes with the identity link. The standard deviations and
# correlations of the deviations start theta, the residual standard deviation of these fits log_sigma
# (gaussian, and lognormal on the log scale). Starting values that are too small leave the optimiser
# in the flat likelihood near zero variance, e.g. with standard deviations of 1 for a response of
# thousands; the deviations include their estimation error, which errs on the safe side.
varstart <- function(dd, start, re){
  npi <- dd$npi
  eta <- sapply(1:5, function(k) as.vector(dd$Xs[[k]] %*% start[npi == k]))
  eta <- matrix(eta, nrow=length(dd$y))
  meanfun <- function(eta) gcmeanR(dd$mod, dd$x, linkfun(dd$lnk[1], eta[, 1], inverse=TRUE), linkfun(dd$lnk[2], eta[, 2], inverse=TRUE),
                                    linkfun(dd$lnk[3], eta[, 3], inverse=TRUE), linkfun(dd$lnk[4], eta[, 4], inverse=TRUE),
                                    linkfun(dd$lnk[5], eta[, 5], inverse=TRUE))
  mu <- meanfun(eta)
  # response on the scale of the mean (proportions for binomial counts), on the log scale for lognormal
  binom <- dd$family %in% c("binomial", "betabinomial")
  lognorm <- dd$family == "lognormal"
  y <- if (binom) dd$y/pmax(dd$bn, 1) else dd$y
  w <- dd$w * (if (binom) dd$bn else 1)
  r <- if (lognorm) log(y) - log(mu) else y - mu
  # derivatives of the mean (log mean for lognormal) w.r.t. the linear predictors of b1-b5
  rpars <- sort(unique(unlist(lapply(re$terms, function(tm) tm$pars))))
  D <- matrix(0, length(y), 5)
  for (k in rpars){
    h <- 1e-5 * pmax(1, abs(eta[, k]))
    ep <- em <- eta
    ep[, k] <- eta[, k] + h
    em[, k] <- eta[, k] - h
    D[, k] <- (if (lognorm) log(meanfun(ep)) - log(meanfun(em)) else meanfun(ep) - meanfun(em)) / (2*h)
  }
  ok <- is.finite(r) & is.finite(w) & w > 0 & apply(is.finite(D), 1, all)
  # fallback scale of the random coefficients of each parameter, on the scale of the link: 10% of the
  # asymptote b1, 10% of the scale b3 for the time b2, and 10% of the linear predictor for b3
  b <- vapply(1:5, function(k) linkfun(dd$lnk[k], start[npi == k][1], inverse=TRUE), numeric(1))
  fscale <- vapply(1:5, function(k){
    s <- switch(k, abs(b[1]), abs(b[3]), abs(linkfun(dd$lnk[3], b[3])), abs(b[4]), abs(b[5]))
    if (!is.finite(s) || s == 0) 1 else 0.1*s
  }, numeric(1))

  zoff <- 0
  rss <- 0
  rdf <- 0
  theta <- unlist(lapply(re$terms, function(tm){
    np <- length(tm$pars)
    nl <- length(tm$levels)
    nc <- length(tm$cnms)
    Zt <- lapply(seq_len(nc), function(c) re$Z[, zoff + (c - 1)*nl + seq_len(nl), drop=FALSE])
    zoff <<- zoff + nc*nl
    # fallback for slopes: per unit of the covariate
    cov <- vapply(seq_len(np), function(j){
      z <- Zt[[tm$zcol[j]]]@x
      if (tm$cnms[tm$zcol[j]] == "(Intercept)" || length(z) < 2 || sd(z) == 0) 1 else sd(z)
    }, numeric(1))
    fb <- fscale[tm$pars] / cov
    est <- matrix(NA, nl, np)
    for (l in seq_len(nl)){
      M <- sapply(seq_len(np), function(j) as.vector(Zt[[tm$zcol[j]]][, l]) * D[, tm$pars[j]])
      M <- matrix(M, ncol=np)
      i <- which(ok & rowSums(M != 0) > 0)
      if (length(i) <= np) next
      fit <- lm.wfit(M[i, , drop=FALSE], r[i], w[i])
      if (fit$rank < np) next
      est[l, ] <- fit$coefficients
      rss <<- rss + sum(w[i] * fit$residuals^2)
      rdf <<- rdf + length(i) - np
    }
    use <- rowSums(is.na(est)) == 0
    sds <- if (sum(use) >= 2) apply(est[use, , drop=FALSE], 2, sd) else rep(NA, np)
    # at least 10% of the fallback scale, to start away from zero variance
    sds <- ifelse(is.finite(sds), pmax(sds, 0.1*fb), fb)
    cpar <- numeric(np*(np - 1)/2)
    if (np > 1 && sum(use) > np){
      R <- suppressWarnings(cor(est[use, , drop=FALSE]))
      if (all(is.finite(R))){
        # shrink towards zero, away from singular correlation matrices
        R <- 0.9*R + 0.1*diag(np)
        C <- t(chol(R))
        L <- C / diag(C)
        cpar <- t(L)[upper.tri(L)]
      }
    }
    c(log(sds), cpar)
  }))
  ls <- if (rdf > 0 && rss > 0) log(sqrt(rss/rdf)) else 0
  list(theta=unname(theta), log_sigma=ls)
}


# a random effect standard deviation that ends close to zero may be stuck in the flat likelihood near
# zero variance rather than estimated as zero: refit once with these standard deviations at their
# starting values, and keep the better fit
retryvariances <- function(obj, opt, terms, lower, upper, ctrl){
  it <- which(names(obj$par) == "theta")
  isd <- it[unlist(lapply(terms, function(tm){
    np <- length(tm$pars)
    c(rep(TRUE, np), rep(FALSE, np*(np - 1)/2))
  }))]
  small <- isd[opt$par[isd] < obj$par[isd] + log(1e-3)]
  if (length(small) == 0 || !is.finite(opt$objective)) return(opt)
  p <- opt$par
  p[small] <- obj$par[small]
  opt2 <- tryCatch(suppressWarnings(nlminb(start=p, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl)),
                   error=function(e) NULL)
  if (!is.null(opt2) && is.finite(opt2$objective) && opt2$objective < opt$objective - 1e-6) opt2 else opt
}


# right-hand side formulas X1-X5 for b1-b5 from fform: a formula b1 + b2 ~ Treatment, a list of
# such formulas, or a list of five right-hand side formulas
fforms <- function(fform){
  res <- rep(list(~ 1), 5)
  names(res) <- paste("X", 1:5, sep="")
  if (is.null(fform)) return(res)
  if (inherits(fform, "formula")) fform <- list(fform)
  if (!is.list(fform) || !all(sapply(fform, inherits, "formula"))) stop("fform must be a formula or a list of formulas.")
  twosided <- sapply(fform, length) == 3
  if (!any(twosided)){
    if (length(fform) != 5) stop("A list of right-hand side formulas in fform needs one formula for each of b1 to b5.")
    res[] <- fform
    return(res)
  }
  if (!all(twosided)) stop("Each formula in fform needs parameters on the left-hand side, e.g. b1 + b2 ~ Treatment.")
  done <- integer(0)
  for (f in fform){
    pnames <- all.vars(f[[2]])
    if (!all(pnames %in% paste("b", 1:5, sep=""))) stop("The left-hand side of a formula in fform can only contain b1 to b5.")
    pars <- match(pnames, paste("b", 1:5, sep=""))
    if (any(pars %in% done)) stop("Each parameter can only appear once in fform.")
    done <- c(done, pars)
    rhs <- f
    rhs[2] <- NULL
    res[pars] <- list(rhs)
  }
  res
}

# coefficient names of the form b2:Stressdrought
coefnames <- function(Xs){
  unlist(lapply(1:5, function(i) paste("b", i, ":", colnames(Xs[[i]]), sep="")))
}

# random effect terms from rform, a formula or a list of formulas of the form b1 + b2 ~ 1 | group
# (random intercepts) or b2 ~ x | group (random intercepts and slopes in x); nested grouping factors
# (1 | a/b) give one term per level of nesting. A term has a random coefficient for each parameter
# and each column of its design (intercept, slopes), correlated within the term: pars, the parameter
# of each coefficient, zcol, its column of the design, and labels, e.g. b2 and b2:x; xterms, the
# terms, factor levels and contrasts of the design, for the design matrices of new data. Z has a
# block for each column of each term, with a column for each level.
reterms <- function(rform, data){
  if (is.null(rform)){
    n <- if (is.numeric(data)) data else nrow(data)
    return(list(terms=list(), Z=Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(n, 0))))
  }
  if (inherits(rform, "formula")) rform <- list(rform)
  terms <- list()
  Zl <- list()
  for (f in rform){
    if (length(f) != 3) stop("Each random effect formula needs parameters on the left-hand side, e.g. b1 + b2 ~ 1 | group.")
    pnames <- all.vars(f[[2]])
    if (!all(pnames %in% paste("b", 1:5, sep=""))) stop("The left-hand side of a random effect formula can only contain b1 to b5.")
    pars <- sort(match(pnames, paste("b", 1:5, sep="")))
    frh <- f
    frh[2] <- NULL
    rt <- mkReTrms(findbars(frh), model.frame(subbars(frh), data=data))
    for (j in seq_along(rt$Ztlist)){
      cn <- rt$cnms[[j]]
      zb <- reblocks(rt$Ztlist[[j]], length(cn))
      # design of the intercept and slopes, from the left-hand side of the bar of the term
      xmf <- model.frame(as.formula(call("~", str2lang(names(rt$Ztlist)[j])[[2]])), data=data)
      xt <- terms(xmf)
      xterms <- list(terms=xt, xlev=.getXlevels(xt, xmf), contrasts=attr(model.matrix(xt, xmf), "contrasts"))
      if (!identical(colnames(designmatrix(xterms, data)), cn)) stop("the design of the random effect term ", names(rt$Ztlist)[j], " is not supported.")
      # each parameter with each column of the design
      cf <- expand.grid(zcol=seq_along(cn), par=pars)
      # name and formula identify the term again in new data (predict)
      terms[[length(terms) + 1]] <- list(pars=cf$par, zcol=cf$zcol, cnms=cn, xterms=xterms,
                                         labels=paste("b", cf$par, ifelse(cn[cf$zcol] == "(Intercept)", "", paste(":", cn[cf$zcol], sep="")), sep=""),
                                         name=names(rt$Ztlist)[j], group=sub("^.*\\| ", "", names(rt$Ztlist)[j]),
                                         levels=zb$levels, formula=frh)
      Zl <- c(Zl, zb$blocks)
    }
  }
  list(terms=terms, Z=do.call(cbind, Zl))
}

# design blocks of a transposed random effect design Zt of lme4 with nc columns per level (rows by
# level, then column): for each column a matrix with a row for each observation and a column for each
# level, and the levels
reblocks <- function(Zt, nc){
  blocks <- lapply(seq_len(nc), function(c) t(Zt[seq(c, nrow(Zt), by=nc), , drop=FALSE]))
  list(blocks=blocks, levels=rownames(Zt)[seq(1, nrow(Zt), by=nc)])
}

# first column in Z (zero-based) of the design block of each random coefficient of all terms
tzbvec <- function(terms){
  zoff <- 0
  unlist(lapply(terms, function(tm){
    nl <- length(tm$levels)
    zb <- zoff + (tm$zcol - 1)*nl
    zoff <<- zoff + length(tm$cnms)*nl
    as.integer(zb)
  }))
}

# labels of the random coefficients of a term, e.g. b2 and b2:x; parameters for fits without them
relabels <- function(tm){
  if (is.null(tm$labels)) paste("b", tm$pars, sep="") else tm$labels
}

# links of the parameters b1-b5 of model m (see gcmodels), named: the defaults of the model (see
# modellinks), replaced by link, a character vector of links named by parameter
parlinks <- function(link, m, model){
  allowed <- modellinks(m)
  out <- vapply(allowed, function(a) a[1], character(1))
  if (!is.null(link)){
    if (!is.character(link)) stop("link must be a character vector.")
    if (is.null(names(link))) stop("link must be named by parameter, e.g. c(b3=\"inverse\").")
    if (any(!(names(link) %in% names(out))) || any(duplicated(names(link)))) stop("the names of link must be distinct parameter names b1-b5.")
    for (p in names(link)){
      if (is.na(linkindex(link[[p]]))) stop("Unknown link: ", link[[p]], "; one of ", paste(linknames, collapse=", "), ".")
    }
    out[names(link)] <- link
  }
  for (p in names(out)){
    if (!(out[[p]] %in% allowed[[p]])) stop("the ", out[[p]], " link is not available for ", p, " of the ", model, " model; ",
                                             if (length(allowed[[p]]) == 1) "it has the " else "one of ", paste(allowed[[p]], collapse=", "),
                                             if (length(allowed[[p]]) == 1) " link." else ".")
  }
  out
}

# start or fix as a vector with one entry per coefficient: an unnamed vector is used as it is,
# a vector named by parameter (e.g. c(b1=0, b5=0)) gives the value for all coefficients of that parameter
parvector <- function(v, npi, default, arg){
  out <- rep(default, length(npi))
  if (is.null(v)) return(out)
  if (is.null(names(v))){
    if (length(v) != length(npi)) stop(arg, " must have one entry for each of the ", length(npi), " model coefficients, or be named by parameter, e.g. c(b5=0).")
    return(v)
  }
  pn <- match(names(v), paste("b", 1:5, sep=""))
  if (any(is.na(pn)) | any(duplicated(pn))) stop("the names of ", arg, " must be distinct parameter names b1-b5.")
  for (i in seq_along(v)) out[npi == pn[i]] <- v[[i]]
  out
}

# keep an asymptote inside the range of the mean of the family: inside (0, 1) for proportions, and
# positive for counts and lognormal responses
linkrange <- function(f, family){
  if (family %in% c("binomial", "betabinomial", "beta")) return(min(max(f, 1e-4), 1 - 1e-4))
  if (family %in% c("poisson", "nbinom2", "nbinom1", "lognormal")) return(max(f, 1e-4))
  f
}
