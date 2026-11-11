# Profile likelihood intervals of the coefficients and of the times of timeto.

test_that("profile intervals of gaussian models are the profile intervals of nls (faba bean heights)", {
  for (m in models){
    f <- gcTMB(height ~ DAS, data=FabaHeights25, model=m)
    n <- fabanls(m)
    ci <- suppressMessages(confint(f, parm=1:3))
    cn <- suppressMessages(confint(n))
    expect_equal(unname(ci), unname(cn), tolerance=1e-3, label=m)
  }
})

test_that("confint gives profile intervals by default and Wald intervals on request", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  expect_message(ci <- confint(f), "profile")
  expect_equal(rownames(ci), c("b1:(Intercept)", "b2:(Intercept)", "b3:(Intercept)", "sigma"))
  cw <- confint(f, method="wald")
  se <- sqrt(diag(vcov(f)))
  expect_equal(unname(cw[1:3, 2]), unname(bcoef(f) + qt(0.975, df.residual(f))*se[1:3]))
  expect_true(all(ci[, 1] < c(bcoef(f), exp(coef(f)[["log_sigma"]])) & ci[, 2] > ci[, 1]))
})

test_that("the profile of the logistic t50 is the profile of b2", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  ct <- suppressMessages(confint(timeto(f, 50)))
  cb <- suppressMessages(confint(f, parm=2))
  expect_equal(unname(ct), unname(cb), tolerance=1e-4)
})

test_that("profile intervals of times are likelihood ratio intervals", {
  for (m in models){
    d <- simgauss(m)
    f <- gcTMB(y ~ day, data=d, model=m)
    q <- 0.9
    ci <- suppressMessages(confint(timeto(f, 100*q)))
    n0 <- nlsfit(m, d)
    rss0 <- deviance(n0)
    cutoff <- nrow(d)*log(1 + qf(0.95, 1, df.residual(f))/df.residual(f))
    # at each limit, the curves through (time, q b1) fitted by nls reach the cutoff of the deviance
    for (tl in ci){
      F <- gcTMB:::gcmodel(match(m, models))$F
      n1 <- nls(y ~ b1*F((day - (tl - b3*Finv[[m]](q)))/b3), data=d, start=list(b1=coef(n0)[[1]], b3=coef(n0)[[3]]))
      expect_equal(nrow(d)*log(deviance(n1)/rss0), cutoff, tolerance=1e-3, label=m)
    }
  }
})

test_that("profiles of absolute times and of times with covariates", {
  d <- rbind(cbind(simgauss("logistic", seed=1), trt="a"),
             cbind(simgauss("logistic", seed=2, b=c(60, 12, 2.5)), trt="b"))
  d$trt <- factor(d$trt)
  f <- gcTMB(y ~ day, fform=b1 + b2 ~ trt, data=d)
  tt <- timeto(f, 30, newdata=data.frame(trt=c("a", "b")), type="absolute")
  ci <- suppressMessages(confint(tt))
  expect_true(all(ci[, 1] < coef(tt) & coef(tt) < ci[, 2]))
  # close to the Wald intervals for this well-identified curve
  expect_equal(unname(ci), unname(confint(tt, method="wald")), tolerance=0.02)
})

test_that("profiles with random effects", {
  d <- simre()
  f <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=d)
  pr <- profile(f, parm=c("b2:(Intercept)", "log_sd:b1|id"), progress=FALSE)
  ci <- confint(pr)
  expect_equal(rownames(ci), c("b2:(Intercept)", "sd:b1|id"))
  expect_true(ci[1, 1] < bcoef(f)[[2]] && bcoef(f)[[2]] < ci[1, 2])
  sdb1 <- gcTMB:::recov(f)[[1]]$sd
  expect_true(ci[2, 1] < sdb1 && sdb1 < ci[2, 2])
  tt <- timeto(f, 90, random=FALSE)
  ct <- suppressMessages(confint(tt, progress=FALSE))
  expect_true(ct[1, 1] < coef(tt) && coef(tt) < ct[1, 2])
})

test_that("intervals beyond the profiled level are missing with a warning", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  pr <- profile(f, parm=2, level=0.5, progress=FALSE)
  expect_warning(ci <- confint(pr, level=0.999), "do not reach")
  expect_true(anyNA(ci))
})

test_that("profiles can be plotted", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  p <- plot(profile(f, progress=FALSE))
  expect_s3_class(p, "ggplot")
  pt <- plot(profile(timeto(f, c(10, 90)), progress=FALSE))
  expect_s3_class(pt, "ggplot")
})
