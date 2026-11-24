# autocorrelation of the residuals within groups: corAR1 and corCAR1

d <- FabaHeights25
d$Plant_ID <- factor(d$Plant_ID)

test_that("the constructors match those of nlme", {
  a <- corAR1(0.3, form=~ DAS | Plant_ID)
  expect_s3_class(a, c("corAR1", "corStruct"))
  expect_equal(as.numeric(a), log(1.3/0.7))
  expect_equal(attr(a, "formula"), ~ DAS | Plant_ID)
  expect_false(attr(a, "fixed"))
  expect_equal(as.numeric(corCAR1(0.3)), qlogis(0.3))
  expect_error(corAR1(1), "between -1 and 1")
  expect_error(corCAR1(0), "between 0 and 1")
  skip_if_not_installed("nlme")
  expect_equal(unclass(corAR1(-0.4, form=~ DAS | Plant_ID, fixed=TRUE)), unclass(nlme::corAR1(-0.4, form=~ DAS | Plant_ID, fixed=TRUE)),
               ignore_attr="formula")
  expect_equal(unclass(corCAR1(0.6, form=~ DAS | Plant_ID)), unclass(nlme::corCAR1(0.6, form=~ DAS | Plant_ID)), ignore_attr="formula")
})

test_that("the likelihood is the multivariate normal density of the residuals of each group", {
  # unequally spaced times with a gap, unsorted, and a group with a single observation
  set.seed(3)
  dd <- data.frame(t=c(5, 1, 2, 4.5, 9, 3, 1, 2, 4, 7), g=c(1, 1, 1, 1, 1, 2, 3, 3, 3, 3))
  dd$y <- 50/(1 + exp(-(dd$t - 4)/1.5)) + rnorm(10, 0, 2)
  for (cs in list(corCAR1(0.5, form=~ t | g), corAR1(-0.3, form=~ 1 | g))){
    fit <- gcTMB(y ~ t, data=dd, truncated=FALSE, correlation=cs, start=c(b1=50, b2=4, b3=1.5))
    mu <- predict(fit)
    s <- sigma(fit)
    phi <- corpars(fit)[["rho"]]
    ll <- sum(sapply(split(seq_len(nrow(dd)), dd$g), function(i){
      tt <- if (inherits(cs, "corCAR1")) dd$t[i] else seq_along(i)
      V <- s^2 * outer(tt, tt, function(a, b) phi^abs(a - b))
      e <- dd$y[i] - mu[i]
      -0.5*(length(e)*log(2*pi) + determinant(V)$modulus + sum(e*solve(V, e)))
    }))
    expect_equal(as.numeric(logLik(fit)), ll, tolerance=1e-8)
  }
})

test_that("the fits without random effects agree with nlme::gnls", {
  skip_if_not_installed("nlme")
  for (k in 1:2){
    mk <- function(v) if (k == 1) corCAR1(v, form=~ DAS | Plant_ID) else corAR1(v, form=~ 1 | Plant_ID)
    fit <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, correlation=mk(0.3))
    expect_equal(fit$opt$convergence, 0)
    cf <- coef(fit)
    g <- suppressWarnings(nlme::gnls(Height ~ b1/(1 + exp(-(DAS - b2)/b3)), data=d, start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]),
                                     correlation=mk(corpars(fit)[["rho"]]), control=nlme::gnlsControl(returnObject=TRUE)))
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(g)), tolerance=1e-7)
    expect_equal(unname(cf[1:3]), unname(coef(g)), tolerance=1e-5)
    expect_equal(corpars(fit)[["rho"]], unname(coef(g$modelStruct$corStruct, unconstrained=FALSE)), tolerance=1e-5)
    # gnls reports sigma with n - p degrees of freedom
    expect_equal(sigma(fit) * sqrt(nobs(fit)/(nobs(fit) - 3)), g$sigma, tolerance=1e-5)
  }
})

test_that("precision weights and the lognormal family agree with nlme::gnls", {
  skip_if_not_installed("nlme")
  d$w <- rep(c(1, 2, 0.5, 4), length.out=nrow(d))
  cs <- corCAR1(0.9, form=~ DAS | Plant_ID)
  fit <- gcTMB(Height ~ DAS, data=d, weights=d$w, truncated=FALSE, correlation=cs)
  cf <- coef(fit)
  g <- suppressWarnings(nlme::gnls(Height ~ b1/(1 + exp(-(DAS - b2)/b3)), data=d, start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]),
                                   correlation=cs, weights=nlme::varFixed(~ 1/w), control=nlme::gnlsControl(returnObject=TRUE)))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(g)), tolerance=1e-7)
  # lognormal: gaussian log(y) with mean log f, and the Jacobian of the log
  fl <- gcTMB(Height ~ DAS, data=d, family="lognormal", truncated=FALSE, correlation=cs)
  cf <- coef(fl)
  gl <- suppressWarnings(nlme::gnls(log(Height) ~ log(b1) - log(1 + exp(-(DAS - b2)/b3)), data=d, start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]),
                                    correlation=cs, control=nlme::gnlsControl(returnObject=TRUE)))
  expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(gl)) - sum(log(d$Height)), tolerance=1e-7)
  expect_equal(unname(cf[1:3]), unname(coef(gl)), tolerance=1e-4)
})

test_that("with random effects, the likelihood is exact for a random asymptote and at least that of nlme", {
  fr <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, truncated=FALSE, correlation=corCAR1(form=~ DAS | Plant_ID))
  cf <- coef(fr)
  s <- sigma(fr)
  phi <- corpars(fr)[["rho"]]
  sdu <- sqrt(VarCorr(fr)[[1]][1])
  # the logistic curve is linear in the random asymptote: the marginal distribution of each plant is
  # multivariate normal, and the Laplace approximation exact
  ll <- sum(sapply(split(d, d$Plant_ID), function(g){
    Fz <- 1/(1 + exp(-(g$DAS - cf[[2]])/cf[[3]]))
    V <- sdu^2 * tcrossprod(Fz) + s^2 * outer(g$DAS, g$DAS, function(a, b) phi^abs(a - b))
    e <- g$Height - cf[[1]]*Fz
    -0.5*(length(e)*log(2*pi) + determinant(V)$modulus + sum(e*solve(V, e)))
  }))
  expect_equal(as.numeric(logLik(fr)), ll, tolerance=1e-8)
  skip_if_not_installed("nlme")
  nm <- nlme::nlme(Height ~ b1/(1 + exp(-(DAS - b2)/b3)), data=d, fixed=b1 + b2 + b3 ~ 1, random=nlme::pdDiag(b1 ~ 1), groups=~ Plant_ID,
                   start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]), correlation=corCAR1(phi, form=~ DAS | Plant_ID), method="ML",
                   control=nlme::nlmeControl(returnObject=TRUE, msMaxIter=500))
  # nlme's alternating algorithm stops close to the optimum
  expect_true(as.numeric(logLik(fr)) >= as.numeric(logLik(nm)) - 1e-6)
  expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(nm)), tolerance=1e-3)
  expect_equal(corpars(fr)[["rho"]], unname(coef(nm$modelStruct$corStruct, unconstrained=FALSE)), tolerance=1e-2)
})

test_that("the grouping defaults to the innermost grouping factor of the random effects", {
  f1 <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corCAR1(form=~ DAS))
  f2 <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance=1e-8)
  expect_equal(f1$corgroup, "Plant_ID")
  expect_output(print(f1), "continuous-time AR\\(1\\) in DAS within Plant_ID")
  # without random effects and without a group: a single series, which has duplicated times here
  expect_error(gcTMB(Height ~ DAS, data=d, correlation=corCAR1(form=~ DAS)), "unique values within groups")
})

test_that("checks of the correlation structures", {
  expect_error(gcTMB(Height ~ DAS, data=d, family="poisson", correlation=corAR1(form=~ 1 | Plant_ID)), "gaussian and lognormal")
  d$t <- d$DAS + 0.5*(d$Plant_ID == "6")
  expect_error(gcTMB(Height ~ DAS, data=d, correlation=corAR1(form=~ t/7 | Plant_ID)), "integer valued")
  expect_error(gcTMB(Height ~ DAS, data=d, correlation=list(1)), "corAR1\\(\\) or corCAR1\\(\\)")
  expect_error(gcTMB(Height ~ DAS, data=d, correlation=corAR1(form=~ DAS | nope)), "nope")
  # F-tests need independent residuals
  f0 <- gcTMB(Height ~ DAS, data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  f1 <- gcTMB(Height ~ DAS, fform=b1 ~ factor(Treatment), data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  expect_error(anova(f0, f1, test="F"), "without autocorrelation")
  expect_null(f0$rdf)
})

test_that("a fixed rho, and the parameter in coef, confint and summary", {
  ff <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corCAR1(0.5, form=~ DAS | Plant_ID, fixed=TRUE))
  expect_equal(corpars(ff), c(rho=0.5))
  expect_false(any(grepl("cor_par", names(coef(ff)))))
  expect_output(print(ff), "rho: 0.5 \\(fixed\\)")
  fe <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corAR1(form=~ 1 | Plant_ID))
  expect_equal(tail(names(coef(fe)), 1), "cor_par:AR1")
  expect_equal(corpars(fe)[["rho"]], tanh(coef(fe)[["cor_par:AR1"]]/2))
  ci <- confint(fe, parm="cor_par:AR1", method="wald")
  expect_equal(rownames(ci), "rho")
  expect_true(ci[1] < corpars(fe) && corpars(fe) < ci[2])
  expect_output(print(summary(fe)), "Correlation structure: AR\\(1\\) in the order of the observations within Plant_ID")
  expect_length(corpars(gcTMB(Height ~ DAS, data=d)), 0)
  # nlme's constructors give the same fit
  skip_if_not_installed("nlme")
  fn <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=nlme::corAR1(form=~ 1 | Plant_ID))
  expect_equal(as.numeric(logLik(fn)), as.numeric(logLik(fe)), tolerance=1e-8)
})

test_that("simulated residuals have the correlation of the fit, and normalized residuals are independent", {
  fit <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  phi <- corpars(fit)[["rho"]]
  s <- sigma(fit)
  dat <- fit$obj$env$data
  sims <- simulate(fit, nsim=400, seed=1, random=FALSE)
  # residuals of the simulations around the curves of the fitted plants, standardised
  R <- (as.matrix(sims) - predict(fit))/s
  expect_equal(mean(R^2), 1, tolerance=0.03)
  has <- dat$corprev >= 0
  # pairs of consecutive observations 7 days apart (most of the design)
  i <- which(has & dat$cordt == 7)
  expect_gt(length(i), 20)
  r7 <- mean(R[i, ] * R[dat$corprev[i] + 1, ])
  expect_equal(r7, phi^7, tolerance=0.05)
  # the fit to the observed data: normalized residuals uncorrelated within plants
  rn <- residuals(fit, type="normalized")
  expect_equal(length(rn), nobs(fit))
  expect_lt(abs(cor(rn[has], rn[dat$corprev[has] + 1])), 0.2)
  rp <- residuals(fit, type="pearson")
  expect_gt(cor(rp[has], rp[dat$corprev[has] + 1]), 0.3)
  # without autocorrelation, the normalized residuals are the pearson residuals
  f0 <- gcTMB(Height ~ DAS, data=d)
  expect_equal(residuals(f0, type="normalized"), residuals(f0, type="pearson"))
})

test_that("simulate with newdata draws correlated residuals within groups", {
  fit <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, correlation=corCAR1(form=~ DAS | Plant_ID))
  phi <- corpars(fit)[["rho"]]
  nd <- expand.grid(DAS=c(10, 11, 13), Plant_ID=c("a", "b"))
  sims <- simulate(fit, nsim=4000, seed=2, newdata=nd)
  R <- (as.matrix(sims) - attr(sims, "mean")[, 1])/sigma(fit)
  C <- cor(t(R))
  expect_equal(C[1, 2], phi, tolerance=0.02)
  expect_equal(C[1, 3], phi^3, tolerance=0.02)
  expect_lt(abs(C[1, 4]), 0.05)
})

test_that("rho is recovered from data simulated with known parameters", {
  set.seed(11)
  # 30 curves, 12 times each, CAR(1) residuals with rho 0.6 per unit of time and a random asymptote
  tt <- c(1, 2, 3.5, 4, 5, 6.5, 7, 8, 9.5, 10, 11, 12)
  sim1 <- function(){
    do.call(rbind, lapply(1:30, function(g){
      b1 <- 50 + rnorm(1, 0, 5)
      e <- as.vector(t(chol(outer(tt, tt, function(a, b) 0.6^abs(a - b)))) %*% rnorm(length(tt)))
      data.frame(g=g, t=tt, y=b1/(1 + exp(-(tt - 6)/1.5)) + 2*e)
    }))
  }
  est <- t(replicate(10, {
    fit <- gcTMB(y ~ t, rform=b1 ~ 1 | g, data=sim1(), truncated=FALSE, correlation=corCAR1(form=~ t | g))
    c(corpars(fit), sigma=sigma(fit), b2=coef(fit)[[2]])
  }))
  expect_equal(mean(est[, "rho"]), 0.6, tolerance=0.05)
  expect_equal(mean(est[, "sigma"]), 2, tolerance=0.05)
  expect_equal(mean(est[, "b2"]), 6, tolerance=0.02)
})

test_that("the methods work with autocorrelation", {
  fit <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  nd <- data.frame(DAS=c(20, 40), Plant_ID=d$Plant_ID[1])
  p <- predict(fit, newdata=nd, se.fit=TRUE, interval="prediction")
  expect_true(all(is.finite(as.matrix(p))))
  bt <- bootstrap(fit, nsim=5, seed=1, progress=FALSE)
  expect_true(all(is.finite(bt$estimates[!is.na(bt$estimates[, 1]), "cor_par:CAR1"])))
  fit0 <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d)
  expect_gt(anova(fit0, fit)$LR[2], 0)
  # zero weights: the observation is left out of the chain
  w <- rep(1, nrow(d))
  w[d$DAS == 28] <- 0
  fw <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, weights=w, correlation=corCAR1(form=~ DAS | Plant_ID))
  fs <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d[w > 0, ], correlation=corCAR1(form=~ DAS | Plant_ID))
  expect_equal(as.numeric(logLik(fw)), as.numeric(logLik(fs)), tolerance=1e-6)
  expect_s3_class(update(fit, sigma=2), "gcTMB")
  expect_equal(corpars(update(fit, sigma=2)), corpars(fit))
})

test_that("update sets rho, keeping the other parameters", {
  fit <- gcTMB(Height ~ DAS, data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  fu <- update(fit, rho=0.7)
  expect_equal(corpars(fu), c(rho=0.7))
  expect_equal(fixef(fu), fixef(fit))
  expect_equal(sigma(fu), sigma(fit))
  expect_true(is.na(converged(fu)))
  expect_output(print(fu), "rho: 0.7")
  # the log-likelihood is the multivariate normal density at the new rho
  mu <- predict(fit)
  ll <- sum(vapply(split(seq_len(nrow(d)), d$Plant_ID), function(i){
    S <- sigma(fit)^2 * 0.7^abs(outer(d$DAS[i], d$DAS[i], "-"))
    L <- t(chol(S))
    z <- forwardsolve(L, d$Height[i] - mu[i])
    -0.5*sum(z^2) - sum(log(diag(L))) - 0.5*length(i)*log(2*pi)
  }, numeric(1)))
  expect_equal(as.numeric(logLik(fu)), ll, tolerance=1e-8)
  # simulations use the new rho
  sims <- simulate(fu, nsim=1, seed=1)
  expect_equal(length(sims[[1]]), nrow(d))
  # AR(1) allows negative values; errors
  fa <- gcTMB(Height ~ DAS, data=d, correlation=corAR1(form=~ 1 | Plant_ID))
  expect_equal(corpars(update(fa, rho=-0.3)), c(rho=-0.3))
  expect_error(update(fit, rho=-0.3), "between 0 and 1")
  expect_error(update(fa, rho=1), "between -1 and 1")
  expect_error(update(gcTMB(Height ~ DAS, data=d), rho=0.5), "no autocorrelation")
  expect_error(update(fit, rho=0.5, model="gompertz"), "not both")
})
