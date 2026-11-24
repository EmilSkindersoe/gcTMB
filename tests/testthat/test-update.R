fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)

test_that("update with the estimates reproduces the fit", {
  fu <- update(fit, fixef=fixef(fit), VarCorr=VarCorr(fit), sigma=sigma(fit))
  expect_equal(as.numeric(logLik(fu)), as.numeric(logLik(fit)), tolerance=1e-6)
  expect_equal(coef(fu), coef(fit), tolerance=1e-6)
  expect_equal(VarCorr(fu), VarCorr(fit), tolerance=1e-6)
  expect_true(all(is.na(vcov(fu))))
  expect_true(is.na(converged(fu)))
  # the original fit is unchanged
  expect_true(isTRUE(converged(fit)))
})

test_that("update sets coefficients, covariance matrices and sigma", {
  vc <- VarCorr(fit)$Plant_ID
  vc2 <- vc * c(2, 1)  # not symmetric
  expect_error(update(fit, VarCorr=list(Plant_ID=vc2)), "symmetric")
  vc2 <- diag(c(2, 3)) %*% vc %*% diag(c(2, 3))
  fu <- update(fit, fixef=c(b1=150, b5=4), VarCorr=vc2, sigma=3)
  expect_equal(unname(fixef(fu)[c("b1:(Intercept)", "b5:(Intercept)")]), c(150, 4))
  expect_equal(unname(VarCorr(fu)$Plant_ID), unname(vc2), tolerance=1e-8, ignore_attr=TRUE)
  expect_equal(sigma(fu), 3)
  expect_equal(unname(fixef(fu)["b2:(Intercept)"]), unname(fixef(fit)["b2:(Intercept)"]))
  expect_error(update(fit, fixef=c(b9=1)), "unknown coefficient")
  expect_output(print(fu), "set by update")
  # refit when no parameters are given
  fg <- update(fit, model="gompertz")
  expect_equal(fg$model, "gompertz")
})

test_that("simulate at new times and for new groups", {
  nd <- expand.grid(DAS=c(10, 30, 50), Plant_ID=c("a", "b"))
  s <- simulate(fit, nsim=3, seed=1, newdata=nd)
  expect_equal(dim(s), c(6, 3))
  expect_equal(simulate(fit, nsim=3, seed=1, newdata=nd), s, ignore_attr=TRUE)
  # fitted groups with random=FALSE and almost no residual noise: the fitted curves
  fu <- update(fit, sigma=1e-6)
  nd <- data.frame(DAS=c(20, 40, 60), Plant_ID=FabaHeights25$Plant_ID[1])
  s <- simulate(fu, nsim=1, newdata=nd, random=FALSE)
  expect_equal(s[[1]], unname(predict(fu, newdata=nd)), tolerance=1e-4)
  # new groups: the variance of the asymptote between curves at a late time
  nd <- data.frame(DAS=1000, Plant_ID=as.character(1:4000))
  s <- simulate(fit, nsim=1, seed=2, newdata=nd)[[1]]
  vb1 <- VarCorr(fit)$Plant_ID[1, 1]
  expect_equal(var(s), vb1 + sigma(fit)^2, tolerance=0.1)
  expect_equal(mean(s), unname(fixef(fit)[1]), tolerance=0.01)
})

test_that("simulate returns the simulated random effects and true means", {
  # fitted data: the true means are the curves with the simulated random effects
  s <- simulate(fit, nsim=2, seed=3)
  re <- attr(s, "ranef")$sim_2$Plant_ID
  expect_equal(dim(re), dim(ranef(fit)$Plant_ID))
  fu <- update(fit, sigma=1e-8)
  nd <- FabaHeights25[FabaHeights25$Plant_ID == FabaHeights25$Plant_ID[1], ]
  pars <- fixef(fit)[1:3] + c(unlist(re[as.character(nd$Plant_ID[1]), ]), 0)
  fp <- update(fit, fixef=c(b1=pars[[1]], b2=pars[[2]], b3=pars[[3]]))
  mu <- attr(s, "mean")$sim_2[FabaHeights25$Plant_ID == FabaHeights25$Plant_ID[1]]
  expect_equal(mu, unname(predict(fp, newdata=nd, random=FALSE)), tolerance=1e-8)
  # random=FALSE: the fitted random effects
  s0 <- simulate(fit, nsim=1, seed=3, random=FALSE)
  expect_equal(attr(s0, "ranef")$sim_1, ranef(fit), ignore_attr=TRUE, tolerance=1e-8)
  expect_equal(attr(s0, "mean")$sim_1, unname(fitted(fit)), tolerance=1e-6)
  # newdata: the means are the observations without residual noise
  nd <- expand.grid(DAS=c(10, 40, 80), Plant_ID=c("a", "b"))
  s1 <- simulate(fit, nsim=1, seed=4, newdata=nd)
  s2 <- simulate(update(fit, sigma=1e-8), nsim=1, seed=4, newdata=nd)
  expect_equal(attr(s1, "mean")$sim_1, s2$sim_1, tolerance=1e-6)
  expect_equal(rownames(attr(s1, "ranef")$sim_1$Plant_ID), c("a", "b"))
})

test_that("random effects fits start from the fit without random effects when the Laplace approximation fails", {
  # starting values of the logistic self-starter, at which the inner optimisation fails
  st <- c(b1 = 116.62856, b2 = 32.31784, b3 = 10.97301)
  f1 <- gcTMB(Height ~ DAS, rform = b1 + b2 + b3 ~ 1 | Plant_ID, data = FabaHeights25, start = st)
  f2 <- gcTMB(Height ~ DAS, rform = b1 + b2 + b3 ~ 1 | Plant_ID, data = FabaHeights25)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
})

test_that("control sdreport=FALSE gives the same estimates without standard errors", {
  fn <- gcTMB(Height ~ DAS, rform = b1 + b2 ~ 1 | Plant_ID, data = FabaHeights25, control = list(sdreport = FALSE))
  expect_equal(coef(fn), coef(fit))
  expect_equal(as.numeric(logLik(fn)), as.numeric(logLik(fit)))
  expect_equal(ranef(fn), ranef(fit))
  expect_equal(VarCorr(fn), VarCorr(fit))
  expect_equal(sigma(fn), sigma(fit))
  expect_true(all(is.na(vcov(fn))))
  expect_true(isTRUE(converged(fn)))
  expect_output(print(summary(fn)), "not computed")
  expect_equal(predict(fn, newdata = data.frame(DAS = 30), random = FALSE),
               predict(fit, newdata = data.frame(DAS = 30), random = FALSE))
  s <- simulate(fn, nsim = 1, seed = 1, newdata = data.frame(DAS = 30, Plant_ID = "new"))
  expect_equal(dim(s), c(1, 1))
  f0 <- gcTMB(Height ~ DAS, data = FabaHeights25, control = list(sdreport = FALSE))
  expect_equal(coef(f0), coef(gcTMB(Height ~ DAS, data = FabaHeights25)))
})

test_that("control restart=TRUE restarts a fit that did not converge, and converged accepts a stable optimum", {
  # three plants and three correlated random effects: a flat likelihood, where nlminb reports
  # false convergence at the maximum
  f3 <- gcTMB(Height ~ DAS, rform = b1 + b2 + b3 ~ 1 | Plant_ID, data = FabaHeights25)
  nd <- expand.grid(DAS = seq(10, 80, length.out = 8), Plant_ID = paste0("p", 1:3))
  nd$Height <- simulate(f3, newdata = nd, seed = 1)[[1]]
  ctrl <- list(sdreport = FALSE)
  f0 <- suppressWarnings(gcTMB(Height ~ DAS, rform = b1 + b2 + b3 ~ 1 | Plant_ID, data = nd, control = ctrl))
  f1 <- suppressWarnings(gcTMB(Height ~ DAS, rform = b1 + b2 + b3 ~ 1 | Plant_ID, data = nd,
                               control = c(ctrl, restart = TRUE)))
  expect_null(f0$restart)
  # whether nlminb reports false convergence here depends on the platform's floating-point
  # arithmetic (it does on Windows, not on Linux and macOS); the test below forces a restart
  if (isTRUE(converged(f0))) skip("the fit converges on this platform, so there is nothing to restart")
  expect_false(isTRUE(converged(f0)))
  expect_true(is.list(f1$restart))
  expect_lt(abs(f1$restart$improvement), 1e-6)
  expect_gte(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_true(isTRUE(converged(f1)))
  expect_false(isTRUE(converged(f1, lltol = 0)))
  # a fit that converges is not restarted
  expect_null(gcTMB(Height ~ DAS, rform = b1 + b2 ~ 1 | Plant_ID, data = FabaHeights25,
                    control = list(restart = TRUE))$restart)
})

test_that("control restart=TRUE restarts a fit stopped by the iteration limit", {
  # two iterations are too few on every platform, so the optimiser does not report convergence
  f1 <- suppressWarnings(gcTMB(Height ~ DAS, data = FabaHeights25,
                               control = list(iter.max = 2, sdreport = FALSE, restart = TRUE)))
  expect_true(is.list(f1$restart))
  expect_false(f1$restart$convergence == 0)
  expect_true(is.finite(f1$restart$improvement))
  expect_gte(f1$restart$improvement, 0)
})

test_that("singular fits count as converged unless singular=TRUE", {
  fs <- gcTMB(Height ~ DAS, rform = b1 + b2 ~ 1 | Plant_ID, data = FabaHeights25)
  fs$ssdrl[rownames(fs$ssdrl) == "theta", 1][1] <- log(1e-10)
  expect_false(isTRUE(converged(fs, singular = TRUE)))
  expect_match(attr(converged(fs, singular = TRUE), "reasons"), "singular", all = FALSE)
  expect_false(any(grepl("singular", attr(converged(fs), "reasons"))))
})

test_that("update can set every parameter of the model", {
  # treatment effects, fixed b4 and b5, correlated random effects, sigma and autocorrelation
  dd <- FabaHeights25
  dd$trt <- factor(dd$Treatment)
  f <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, rform=b1 + b2 ~ 1 | Plant_ID, data=dd,
             correlation=corCAR1(form=~ DAS | Plant_ID))
  fe <- fixef(f)
  vc <- VarCorr(f)$Plant_ID
  newfe <- c(fe + 1, b4=2, b5=1.5)
  fu <- update(f, fixef=newfe, VarCorr=list(Plant_ID=2*vc), sigma=2*sigma(f), rho=0.5)
  expect_equal(fixef(fu)[names(fe)], fe + 1)
  expect_equal(unname(fu$obj$env$parList(fu$opt$par)$b4), 2)
  expect_equal(unname(fu$obj$env$parList(fu$opt$par)$b5), 1.5)
  expect_equal(VarCorr(fu)$Plant_ID, 2*vc, ignore_attr=TRUE, tolerance=1e-8)
  expect_equal(sigma(fu), 2*sigma(f))
  expect_equal(corpars(fu), c(rho=0.5))
  # all the entries of coef are changed as expected
  cf <- coef(f)
  ex <- cf
  b <- grepl("^b[1-5]:", names(cf))
  ex[b] <- cf[b] + 1
  sdn <- grepl("^log_sd:", names(cf))
  ex[sdn] <- cf[sdn] + log(2)/2
  ex["log_sigma"] <- cf["log_sigma"] + log(2)
  ex["cor_par:CAR1"] <- 0
  expect_equal(coef(fu), ex, tolerance=1e-8)
  # a fixed autocorrelation and the dispersion of another family can be set too
  ff <- update(f, correlation=corCAR1(0.8, form=~ DAS | Plant_ID, fixed=TRUE))
  expect_equal(corpars(update(ff, rho=0.3)), c(rho=0.3))
  fn <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="nbinom2", fix=c(b4=FALSE))
  expect_equal(sigma(update(fn, sigma=3)), 3)
})
