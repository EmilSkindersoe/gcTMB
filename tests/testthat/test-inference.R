# Bootstrap, likelihood ratio tests, simulation and diagnostic plots.

test_that("bootstrap refits equal gcTMB fits to the simulated responses", {
  d <- simgauss("gompertz")
  f <- gcTMB(y ~ day, truncated=FALSE, data=d, model="gompertz")
  bt <- bootstrap(f, nsim=5, seed=1, progress=FALSE)
  expect_s3_class(bt, "bootstrap.gcTMB")
  sims <- simulate(f, nsim=5, seed=1)
  f1 <- gcTMB(y ~ day, truncated=FALSE, data=transform(d, y=sims[[1]]), model="gompertz")
  expect_equal(unname(bt$estimates[1, ]), unname(coef(f1)), tolerance=1e-5)
  expect_output(print(bt), "5 of 5 fits converged")
})

test_that("bootstrap intervals of coefficients and times", {
  f <- gcTMB(y ~ day, truncated=FALSE, data=simgauss("logistic"))
  bt <- bootstrap(f, nsim=200, seed=2, progress=FALSE)
  cb <- confint(f, method="boot", boot=bt)
  cw <- confint(f, method="wald")
  expect_equal(unname(cb[1:3, ]), unname(cw[1:3, ]), tolerance=0.05)
  tt <- timeto(f, c(10, 90))
  tb <- confint(tt, boot=bt)
  # the times of the bootstrap coefficients
  q <- bt$estimates[, 2] + bt$estimates[, 3]*qlogis(0.9)
  expect_equal(unname(tb[2, ]), unname(quantile(q, c(0.025, 0.975))), tolerance=1e-8)
  expect_error(confint(f, boot=bootstrap(gcTMB(y ~ day, truncated=FALSE, data=simgauss("gompertz"), model="gompertz"), nsim=2, progress=FALSE)), "different model")
})

test_that("anova gives likelihood ratio tests and F-tests", {
  d <- rbind(cbind(simgauss("logistic", seed=1), trt="a"),
             cbind(simgauss("logistic", seed=2, b=c(60, 10, 2.5)), trt="b"))
  d$trt <- factor(d$trt)
  f0 <- gcTMB(y ~ day, truncated=FALSE, data=d)
  f1 <- gcTMB(y ~ day, truncated=FALSE, fform=b1 ~ trt, data=d)
  a <- anova(f0, f1)
  expect_equal(a$LR[2], 2*(as.numeric(logLik(f1)) - as.numeric(logLik(f0))))
  expect_equal(a$Df[2], 1)
  expect_lt(a[["Pr(>LR)"]][2], 0.001)
  # F-test as anova of nls
  n0 <- nlsfit("logistic", d)
  n1 <- nls(y ~ (a1 + a2*(trt == "b"))/(1 + exp(-(day - b2)/b3)), data=d, start=list(a1=50, a2=10, b2=10, b3=2.5))
  af <- anova(f0, f1, test="F")
  an <- anova(n0, n1)
  expect_equal(af$F[2], an[["F value"]][2], tolerance=1e-6)
  expect_equal(af[["Pr(>F)"]][2], an[["Pr(>F)"]][2], tolerance=1e-6)
  expect_error(anova(f0), "at least two")
  expect_error(anova(f0, gcTMB(y ~ day, truncated=FALSE, data=simgauss("logistic"))), "same data")
})

test_that("anova of models with and without random effects, and bootstrap tests", {
  d <- simre()
  f0 <- gcTMB(y ~ day, truncated=FALSE, rform=b1 ~ 1 | id, data=d)
  f1 <- gcTMB(y ~ day, truncated=FALSE, rform=b1 + b2 ~ 1 | id, data=d)
  a <- anova(f0, f1)
  expect_equal(a$Df[2], 2)
  ab <- anova(f0, f1, test="boot", nsim=10, seed=1, progress=FALSE)
  expect_true(all(c("Pr(boot)") %in% names(ab)))
  expect_length(attr(ab, "boot")[[1]], 10)
})

test_that("simulate has the mean and variance of the fitted model", {
  d <- simre()
  f <- gcTMB(y ~ day, truncated=FALSE, rform=b1 ~ 1 | id, data=d)
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
  f <- gcTMB(y ~ day, truncated=FALSE, rform=b1 + b2 ~ 1 | id, data=d)
  se <- sqrt(diag(vcov(f)))[1:3]
  expect_true(all(abs(bcoef(f) - truepars) < 3*se))
  expect_equal(unname(gcTMB:::recov(f)[[1]]$sd), c(5, 1), tolerance=0.35)
})

test_that("plot gives diagnostic ggplots", {
  f <- gcTMB(y ~ day, truncated=FALSE, rform=b1 ~ 1 | id, data=simre())
  p <- plot(f, which=1:4, draw=FALSE, seed=1)
  expect_equal(names(p), c("fitted", "time", "qq", "ranef", "curves"))
  expect_true(all(vapply(p, inherits, logical(1), "ggplot")))
  expect_equal(p$time$labels$x, "day")
  expect_equal(p$time$labels$title, "Residuals vs day")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(f))
  expect_error(plot(gcTMB(y ~ day, truncated=FALSE, data=simgauss("logistic")), which=4), "no random effects")
})

test_that("the QQ plot of non-normal residuals warns", {
  fp <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE))
  expect_warning(plot(fp, residuals="pearson", draw=FALSE), "not normally distributed")
  expect_warning(plot(fp, residuals="response", which=3, draw=FALSE), "not normally distributed")
  expect_no_warning(plot(fp, draw=FALSE, seed=1))
  expect_no_warning(plot(fp, residuals="pearson", which=1:2, draw=FALSE))
  f <- gcTMB(y ~ day, truncated=FALSE, data=simgauss("logistic"))
  expect_no_warning(plot(f, residuals="pearson", draw=FALSE))
})

test_that("plot with type='curves' draws the curves of the individuals and the population", {
  d <- simre()
  f <- gcTMB(y ~ day, truncated=FALSE, rform=b1 ~ 1 | id, data=d)
  p <- plot(f, type="curves", draw=FALSE)$curves
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # points, the curve of each individual and one population curve
  expect_equal(nrow(b$data[[1]]), nrow(d))
  expect_equal(length(unique(b$data[[2]]$group)), length(unique(d$id)))
  expect_equal(length(unique(b$data[[3]]$group)), 1)
  # the curve of an individual at its observed times is its fitted curve
  ig <- p$layers[[2]]$data
  id1 <- ig[ig$id == ig$id[1], ]
  expect_equal(id1$y, predict(f, newdata=id1))
  # population curve with the random effects at zero
  pg <- p$layers[[3]]$data
  expect_equal(pg$y, predict(f, newdata=pg, random=FALSE))
  # curves at given times
  pt <- plot(f, type="curves", draw=FALSE, times=seq(0, 5, by=0.5))$curves
  expect_equal(sort(unique(pt$layers[[3]]$data$time)), seq(0, 5, by=0.5))
  expect_error(plot(f, type="curves", times=1), "at least two times")

  # covariates: a population curve for each level, coloured
  d$grp <- factor(rep(c("a", "b"), length.out=length(unique(d$id))))[as.integer(factor(d$id))]
  fg <- gcTMB(y ~ day, fform=b1 ~ grp, truncated=FALSE, rform=b1 ~ 1 | id, data=d)
  pgp <- plot(fg, type="curves", draw=FALSE)$curves
  expect_equal(length(unique(pgp$layers[[3]]$data$curve)), 2)
  expect_equal(pgp$labels$colour, "grp")

  # without random effects: only the population curve
  f0 <- gcTMB(y ~ day, truncated=FALSE, data=simgauss("logistic"))
  p0 <- plot(f0, type="curves", draw=FALSE)$curves
  expect_equal(length(p0$layers), 2)
  # titles: a single curve without random effects, curves with them or with covariates
  expect_equal(p0$labels$title, "Fitted curve")
  expect_equal(p$labels$title, "Fitted curves: individuals and population")
  dg <- simgauss("logistic")
  dg$g <- factor(rep(c("a", "b"), length.out=nrow(dg)))
  pg0 <- plot(gcTMB(y ~ day, fform=b1 ~ g, truncated=FALSE, data=dg), type="curves", draw=FALSE)$curves
  expect_equal(pg0$labels$title, "Fitted curves")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(f, type="curves"))
  # validation plots with the curves last
  pb <- plot(f, type=c("validation", "curves"), which=1:4, draw=FALSE, seed=1)
  expect_equal(names(pb), c("fitted", "time", "qq", "ranef", "curves"))
  expect_silent(plot(f, type=c("curves", "validation")))
  expect_equal(names(plot(f, type=c("curves", "validation"), draw=FALSE)), c("fitted", "time", "qq", "curves"))
  # validation: the curves are added last, as a single graph without a legend
  pv <- plot(fg, draw=FALSE)
  expect_equal(names(pv), c("fitted", "time", "qq", "curves"))
  expect_equal(pv$curves$theme$legend.position, "none")
  expect_s3_class(pv$curves$facet, "FacetNull")
  # curves alone: the individuals with the data and the population curves in two panels, one legend
  expect_s3_class(pgp$facet, "FacetWrap")
  bg <- ggplot2::ggplot_build(pgp)
  expect_equal(unique(bg$data[[1]]$PANEL), factor(1, levels=1:2))
  expect_equal(unique(bg$data[[2]]$PANEL), factor(1, levels=1:2))
  expect_equal(unique(bg$data[[3]]$PANEL), factor(2, levels=1:2))
  expect_false(identical(pgp$theme$legend.position, "none"))
  # without random effects: a single graph
  expect_s3_class(p0$facet, "FacetNull")
})
