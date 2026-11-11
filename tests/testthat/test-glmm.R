# Some growth models are generalised linear (mixed) models in another parameterisation. Both gcTMB
# and lme4/glmmTMB use the Laplace approximation, so the marginal log-likelihoods agree at the optimum.
#
#   binomial or beta, logistic, b1 = 1, inverse link of b3 (rate r):
#     logit(p) = r t - r b2, random b2: sd = sd(u) r
# Gaussian models with random asymptotes b1 are linear in the random effects, so the Laplace
# approximation of gcTMB and the linearisation of nlme are both exact.

skip_if_not_installed("lme4")

asym1 <- list(start=c(b1=1), fix=c(b1=TRUE), link=c(b3="inverse"))

simbinom <- function(seed=2){
  d <- simgroups(ng=12, seed=seed)
  u <- rnorm(12, 0, 1.5)
  d$n <- 20
  d$y <- rbinom(nrow(d), d$n, plogis(0.5*(d$day - (10 + u[d$id]))))
  d
}

test_that("binomial logistic without random effects agrees with glm", {
  d <- simbinom()
  f <- gcTMB(cbind(y, n - y) ~ day, data=d, family="binomial", start=asym1$start, fix=asym1$fix, link=asym1$link)
  g <- glm(cbind(y, n - y) ~ day, data=d, family=binomial)
  b <- coef(g)
  expect_equal(f$opt$objective, -as.numeric(logLik(g)), tolerance=1e-7)
  expect_equal(unname(bcoef(f)), unname(c(-b[1]/b[2], b[2])), tolerance=1e-5)
  expect_equal(residuals(f), unname(residuals(g, type="pearson")), tolerance=1e-5)
})

test_that("binomial logistic with random b2 agrees with glmer", {
  d <- simbinom()
  f <- gcTMB(cbind(y, n - y) ~ day, rform=b2 ~ 1 | id, data=d, family="binomial",
             start=asym1$start, fix=asym1$fix, link=asym1$link)
  g <- lme4::glmer(cbind(y, n - y) ~ day + (1 | id), data=d, family=binomial)
  b <- lme4::fixef(g)
  sdg <- sqrt(lme4::VarCorr(g)$id[1])
  expect_equal(f$opt$objective, -as.numeric(logLik(g)), tolerance=1e-6)
  expect_equal(unname(bcoef(f)), unname(c(-b[1]/b[2], b[2])), tolerance=1e-4)
  expect_equal(gcTMB:::recov(f)[[1]]$sd, unname(sdg/b[2]), tolerance=1e-3)
})

test_that("beta logistic agrees with glmmTMB, with and without random effects", {
  skip_if_not_installed("glmmTMB")
  d <- simgroups(times=seq(0, 20, 2), ng=12, seed=2)
  u <- rnorm(12, 0, 0.5)
  mu <- plogis(0.3*(d$day - 10) + u[d$id])
  d$y <- rbeta(nrow(d), mu*30, (1 - mu)*30)
  f1 <- gcTMB(y ~ day, data=d, family="beta", start=asym1$start, fix=asym1$fix, link=asym1$link)
  g1 <- glmmTMB::glmmTMB(y ~ day, data=d, family=glmmTMB::beta_family())
  b <- glmmTMB::fixef(g1)$cond
  expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-6)
  expect_equal(unname(bcoef(f1)), unname(c(-b[1]/b[2], b[2])), tolerance=1e-4)
  expect_equal(unname(exp(coef(f1)["log_phi"])), sigma(g1), tolerance=1e-4)
  f2 <- gcTMB(y ~ day, rform=b2 ~ 1 | id, data=d, family="beta", start=asym1$start, fix=asym1$fix, link=asym1$link)
  g2 <- glmmTMB::glmmTMB(y ~ day + (1 | id), data=d, family=glmmTMB::beta_family())
  expect_equal(f2$opt$objective, -as.numeric(logLik(g2)), tolerance=1e-6)
})

test_that("random asymptotes: gcTMB and nlme maximise the same exact likelihood", {
  skip_if_not_installed("nlme")
  for (m in models){
    d <- simre(m, sd2=0)
    f <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=d, model=m)
    n <- nlme::nlme(formulas[[m]], fixed=b1 + b2 + b3 ~ 1, random=b1 ~ 1, groups=~id, data=d,
                    start=setNames(bcoef(f), c("b1", "b2", "b3")), method="ML")
    np <- c(nlme::fixef(n), log(as.numeric(nlme::VarCorr(n)[1, 2])), log(n$sigma))
    expect_equal(as.numeric(f$obj$fn(np)), -as.numeric(logLik(n)), tolerance=1e-6, label=m)
    expect_lte(f$opt$objective, -as.numeric(logLik(n)) + 1e-6)
    expect_equal(unname(coef(f)), unname(np), tolerance=1e-3, label=m)
  }
})

test_that("random times of the maximum growth rate agree approximately with nlme", {
  skip_if_not_installed("nlme")
  d <- simre(sd1=0, sd2=1.5)
  f <- gcTMB(y ~ day, rform=b2 ~ 1 | id, data=d)
  n <- nlme::nlme(formulas$logistic, fixed=b1 + b2 + b3 ~ 1, random=b2 ~ 1, groups=~id, data=d,
                  start=setNames(bcoef(f), c("b1", "b2", "b3")), method="ML")
  se <- sqrt(diag(vcov(f)))[1:3]
  expect_true(all(abs(bcoef(f) - nlme::fixef(n)) < 0.5*se))
  expect_equal(gcTMB:::recov(f)[[1]]$sd, as.numeric(nlme::VarCorr(n)[1, 2]), tolerance=0.05)
})

test_that("random slopes of the asymptote agree with nlme", {
  skip_if_not_installed("nlme")
  # two series per individual, under the conditions x = -1 and x = 1
  d <- simgroups(ng=12, nrep=2, seed=4)
  d$x <- rep(c(-1, 1), length.out=nrow(d))
  u <- cbind(rnorm(12, 0, 4), rnorm(12, 0, 3))
  d$y <- (50 + 5*d$x + u[d$id, 1] + u[d$id, 2]*d$x)/(1 + exp(-(d$day - 10)/2.5)) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, fform=b1 ~ x, rform=b1 ~ x | id, data=d)
  n <- nlme::nlme(y ~ b1/(1 + exp(-(day - b2)/b3)), fixed=list(b1 ~ x, b2 + b3 ~ 1), random=b1 ~ x, groups=~id, data=d,
                  start=unname(bcoef(f)), method="ML")
  expect_lte(f$opt$objective, -as.numeric(logLik(n)) + 1e-5)
  expect_equal(f$opt$objective, -as.numeric(logLik(n)), tolerance=1e-4)
  expect_equal(unname(bcoef(f)), unname(nlme::fixef(n)), tolerance=1e-3)
  # random slopes need their covariates for predictions
  expect_error(predict(f, newdata=data.frame(day=1, id="1")))
  expect_length(predict(f, newdata=data.frame(day=1, id="1", x=0.5)), 1)
})
