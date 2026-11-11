# Bootstrap, likelihood ratio tests, simulation and diagnostic plots.

test_that("bootstrap refits equal gcTMB fits to the simulated responses", {
  d <- simgauss("gompertz")
  f <- gcTMB(y ~ day, data=d, model="gompertz")
  bt <- bootstrap(f, nsim=5, seed=1, progress=FALSE)
  expect_s3_class(bt, "bootstrap.gcTMB")
  sims <- simulate(f, nsim=5, seed=1)
  f1 <- gcTMB(y ~ day, data=transform(d, y=sims[[1]]), model="gompertz")
  expect_equal(unname(bt$estimates[1, ]), unname(coef(f1)), tolerance=1e-5)
  expect_output(print(bt), "5 of 5 fits converged")
})

test_that("bootstrap intervals of coefficients and times", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  bt <- bootstrap(f, nsim=200, seed=2, progress=FALSE)
  cb <- confint(f, method="boot", boot=bt)
  cw <- confint(f, method="wald")
  expect_equal(unname(cb[1:3, ]), unname(cw[1:3, ]), tolerance=0.05)
  tt <- timeto(f, c(10, 90))
  tb <- confint(tt, boot=bt)
  # the times of the bootstrap coefficients
  q <- bt$estimates[, 2] + bt$estimates[, 3]*qlogis(0.9)
  expect_equal(unname(tb[2, ]), unname(quantile(q, c(0.025, 0.975))), tolerance=1e-8)
  expect_error(confint(f, boot=bootstrap(gcTMB(y ~ day, data=simgauss("gompertz"), model="gompertz"), nsim=2, progress=FALSE)), "different model")
})

test_that("anova gives likelihood ratio tests and F-tests", {
  d <- rbind(cbind(simgauss("logistic", seed=1), trt="a"),
             cbind(simgauss("logistic", seed=2, b=c(60, 10, 2.5)), trt="b"))
  d$trt <- factor(d$trt)
  f0 <- gcTMB(y ~ day, data=d)
  f1 <- gcTMB(y ~ day, fform=b1 ~ trt, data=d)
  a <- anova(f0, f1)
  expect_equal(a$Chisq[2], 2*(as.numeric(logLik(f1)) - as.numeric(logLik(f0))))
  expect_equal(a$Df[2], 1)
  expect_lt(a[["Pr(>Chisq)"]][2], 0.001)
  # F-test as anova of nls
  n0 <- nlsfit("logistic", d)
  n1 <- nls(y ~ (a1 + a2*(trt == "b"))/(1 + exp(-(day - b2)/b3)), data=d, start=list(a1=50, a2=10, b2=10, b3=2.5))
  af <- anova(f0, f1, test="F")
  an <- anova(n0, n1)
  expect_equal(af$F[2], an[["F value"]][2], tolerance=1e-6)
  expect_equal(af[["Pr(>F)"]][2], an[["Pr(>F)"]][2], tolerance=1e-6)
  expect_error(anova(f0), "at least two")
  expect_error(anova(f0, gcTMB(y ~ day, data=simgauss("logistic"))), "same data")
})

test_that("anova of models with and without random effects, and bootstrap tests", {
  d <- simre()
  f0 <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=d)
  f1 <- gcTMB(y ~ day, rform=b1 + b2 ~ 1 | id, data=d)
  a <- anova(f0, f1)
  expect_equal(a$Df[2], 2)
  ab <- anova(f0, f1, test="boot", nsim=10, seed=1, progress=FALSE)
  expect_true(all(c("Pr(boot)") %in% names(ab)))
  expect_length(attr(ab, "boot")[[1]], 10)
})

test_that("simulate has the mean and variance of the fitted model", {
  d <- simre()
  f <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=d)
  s <- simulate(f, nsim=2000, seed=1)
  expect_equal(dim(s), c(nrow(d), 2000))
  # new random effects: the marginal variance includes the variance of b1 times F^2
  Fz <- predict(f, random=FALSE)/bcoef(f)[[1]]
  vc <- gcTMB:::recov(f)[[1]]
  v <- residpar(f)^2 + vc$sd^2*Fz^2
  expect_equal(rowMeans(s), predict(f, random=FALSE), tolerance=0.02)
  expect_equal(apply(s, 1, var), v, tolerance=0.15)
  # with the fitted random effects
  s2 <- simulate(f, nsim=2000, seed=1, random=FALSE)
  expect_equal(rowMeans(s2), predict(f), tolerance=0.02)
})

test_that("fixed effects and variance components are recovered", {
  d <- simre(ng=40, sd1=5, sd2=1, seed=11)
  f <- gcTMB(y ~ day, rform=b1 + b2 ~ 1 | id, data=d)
  se <- sqrt(diag(vcov(f)))[1:3]
  expect_true(all(abs(bcoef(f) - truepars) < 3*se))
  expect_equal(unname(gcTMB:::recov(f)[[1]]$sd), c(5, 1), tolerance=0.35)
})

test_that("plot gives diagnostic ggplots", {
  f <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=simre())
  p <- plot(f, which=1:4, draw=FALSE, seed=1)
  expect_equal(names(p), c("fitted", "time", "qq", "ranef"))
  expect_true(all(vapply(p, inherits, logical(1), "ggplot")))
  expect_equal(p$time$labels$x, "Time")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(f))
  expect_error(plot(gcTMB(y ~ day, data=simgauss("logistic")), which=4), "no random effects")
})
