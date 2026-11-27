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
  # update only sets parameters: other arguments, or none, point to refit
  expect_error(update(fit, model="gompertz"), "has no argument model.*use refit")
  expect_error(update(fit), "use refit")
  fg <- refit(fit, model="gompertz")
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
  fu <- update(f, fixef=newfe, VarCorr=list(Plant_ID=2*vc), sigma=2*sigma(f), corpar=0.5)
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
  ff <- refit(f, correlation=corCAR1(0.8, form=~ DAS | Plant_ID, fixed=TRUE))
  expect_equal(corpars(update(ff, corpar=0.3)), c(rho=0.3))
  fn <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="nbinom2", fix=c(b4=FALSE))
  expect_equal(sigma(update(fn, sigma=3)), 3)
})

test_that("refit starts from the data or from the estimates of the fit", {
  d <- FabaHeights25
  d$trt <- factor(d$Treatment)
  f <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, correlation=corCAR1(form=~ DAS | Plant_ID))
  # from the data: the call as stats::update, without starting values
  fd <- refit(f, fform=b1 + b2 ~ trt)
  expect_null(fd$call$start)
  # from the fit: the intercepts of b1 and b2 and b3 carry over, the new coefficients start from the data
  ff <- refit(f, fform=b1 + b2 ~ trt, startfrom="fit")
  sl <- eval(ff$call$start)
  st <- sl$fixef
  expect_length(st, length(ff$start))
  cn <- ff$coefnames
  keep <- cn %in% f$coefnames
  # the estimates, and the values of the fixed coefficients (b4 = 0, b5 = 1)
  pv <- setNames(f$start, f$coefnames)
  pv[!f$fix] <- coef(f)[seq_len(sum(!f$fix))]
  expect_equal(st[keep], unname(pv[cn[keep]]))
  expect_true(all(is.na(st[!keep])))
  # the variance of the random effects and sigma carry over
  expect_equal(sl$VarCorr$Plant_ID, unname(VarCorr(f)$Plant_ID), ignore_attr=TRUE)
  expect_equal(sl$sigma, sigma(f))
  # the correlation structure starts at the estimate of rho
  expect_equal(eval(ff$call$correlation[["value"]]), unname(corpars(f)))
  # both reach the same optimum
  expect_equal(as.numeric(logLik(ff)), as.numeric(logLik(fd)), tolerance=1e-4)
  # a new correlation structure, or a start given, is used as given
  fg <- refit(f, correlation=corExp(form=~ DAS | Plant_ID), startfrom="fit")
  expect_s3_class(eval(fg$call$correlation), "corExp")
  expect_null(fg$call$correlation$value)
  # a new growth model: no coefficients from the fit
  fm <- refit(f, model="gompertz", startfrom="fit")
  expect_null(eval(fm$call$start)$fixef)
  # a variance function starts at its estimates
  fv <- gcTMB(Height ~ DAS, data=d, variance=varPower())
  fv2 <- refit(fv, fform=b1 ~ trt, startfrom="fit")
  expect_equal(eval(fv2$call$variance[["value"]]), unname(varpars(fv)))
  # a new random effect: the variance of b1 carries over, b2 and the correlation start from the data
  fb <- refit(f, rform=b1 + b2 ~ 1 | Plant_ID, startfrom="fit")
  S <- eval(fb$call$start)$VarCorr$Plant_ID
  expect_equal(S[1, 1], VarCorr(f)$Plant_ID[1, 1])
  expect_true(all(is.na(S[-1])))
  # setting parameters and refitting are separate functions
  expect_error(update(f, fform=b1 ~ trt, fixef=c(b1=150)), "use refit")
  expect_error(update(f, fixef=c(b1=150), startfrom="fit"), "has no argument startfrom")
  expect_error(refit(f, fixef=c(b1=150)), "use update")
  # an unnamed formula updates the formula; NULL removes a component
  fn <- refit(f, . ~ ., correlation=NULL)
  expect_null(fn$correlation)
  expect_equal(as.numeric(logLik(fn)), as.numeric(logLik(gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d))))
})

test_that("start can be a list with the coefficients, the covariance matrices of the random effects and sigma", {
  d <- FabaHeights25
  dd <- gcTMBdata(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=d)
  fix <- defaultfix(NULL, dd$npi, dd$mod)
  st <- gcTMBstartvals(dd, defaultstart(parvector(NULL, dd$npi, NA, "start"), dd), fix)
  S <- matrix(c(100, 6, 6, 4), 2)
  sv <- splitstart(list(fixef=NULL, VarCorr=list(Plant_ID=S), sigma=2), dd, b1 + b2 ~ 1 | Plant_ID)
  obj <- gcTMBobj(dd, st, fix, rform=b1 + b2 ~ 1 | Plant_ID, startvar=sv)
  pl <- obj$env$parList()
  expect_equal(pl$log_sigma, log(2))
  expect_equal(pl$theta, covtheta(S, 2, "Plant_ID"))
  # NA entries start from the data: here the variance of b2 and the correlation
  obj0 <- gcTMBobj(dd, st, fix, rform=b1 + b2 ~ 1 | Plant_ID)
  S2 <- matrix(c(100, NA, NA, NA), 2)
  obj2 <- gcTMBobj(dd, st, fix, rform=b1 + b2 ~ 1 | Plant_ID, startvar=list(VarCorr=S2))
  expect_equal(obj2$env$parList()$theta, c(log(10), obj0$env$parList()$theta[2:3]))
  # a fit from a start list reaches the fit from the data
  f0 <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=d)
  f1 <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=d, start=list(fixef=c(b1=120), VarCorr=S, sigma=2))
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)), tolerance=1e-5)
  # errors
  expect_error(gcTMB(Height ~ DAS, data=d, start=list(coef=1)), "only have the components")
  expect_error(gcTMB(Height ~ DAS, data=d, start=list(VarCorr=S)), "needs random effects")
  expect_error(gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", start=list(sigma=1)), "no residual standard deviation")
  expect_error(gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=d, start=list(VarCorr=S)), "must be 1 x 1")
})
